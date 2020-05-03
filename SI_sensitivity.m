close all
clear all

addpath("functions", "result")

% This file simulates the sensitivity curves.

%% User Input
% number of threads
%parpool(4); 
% number of data points per cluster
N_k = 50; 
% Monte Carlo iterations
MC = 5;

% range of outliers
out_range = [-20 20; -20 20]; %range of outliers
% steps between outliers
step_eps = 10; 


% Select combinations of EM and BIC to be simulated
% 1: Gaussian, 2: t, 3: Huber, 4: Tukey
em_bic = [1 1;
          2 2;
          2 4;
          3 3;
          3 4];
% t:
nu = 3;
% Huber:
qH = 0.8;
%Tukey
cT = 4.685;     
      

%% Data generation
% grid for single outlier
x = out_range(1,1):step_eps:out_range(1,2); 
y = out_range(2,1):step_eps:out_range(2,2);
[X, Y] = meshgrid(x,y);

embic_iter = length(em_bic);
eps_iter = numel(X);

for ii_eps = 1:eps_iter
    for ii_mc = 1:MC
        [data(:,:,ii_eps,ii_mc), labels_true, r, N, K_true, mu_true, S_true] = data_31(N_k, 0);

        % replacement outlier
        N_repl = 1;
        index_repl = randperm(N,N_repl).';
        data(index_repl,:,ii_eps,ii_mc) = [X(ii_eps) Y(ii_eps)];
    end
end
L_max = 2*K_true; % search range


%% density definitions
% Huber:
cH = sqrt(chi2inv(qH, r));
bH = chi2cdf(cH^2,r+2) + cH^2/r*(1-chi2cdf(cH^2,r));
aH = gamma(r/2)/pi^(r/2) / ( (2*bH)^(r/2)*(gamma(r/2) - igamma(r/2, cH^2/(2*bH))) + (2*bH*cH^r*exp(-cH^2/(2*bH)))/(cH^2-bH*r) );

g = {@(t)g_gaus(t, r);
     @(t)g_t(t, r, nu);
     @(t)g_huber(t, r, cH, bH, aH)};

rho = {@(t)rho_gaus(t, r);
       @(t)rho_t(t, r, nu);
       @(t)rho_huber(t, r, cH, bH, aH);
       @(t)rho_tukey(t, r, cT)};

psi = {@(t)psi_gaus(t);
       @(t)psi_t(t, r, nu);
       @(t)psi_huber(t, r, cH, bH);
       @(t)psi_tukey(t, cT)};

eta = {@(t)eta_gaus(t);
       @(t)eta_t(t, r, nu); 
       @(t)eta_huber(t, r, cH, bH);
       @(t)eta_tukey(t, cT)};
   
tic
for ii_eps = 1:eps_iter
    parfor ii_mc = 1:MC
        bic = zeros(L_max, 3, embic_iter);   
        like = zeros(L_max, 3, embic_iter);   
        pen = zeros(L_max, 3, embic_iter); 
        for iEmBic = 1:embic_iter
            for ll = 1:L_max
                %% EM
                [mu_est, S_est, t, R] = EM_RES(data(:,:,ii_eps,ii_mc), ll, g{em_bic(iEmBic,1)}, psi{em_bic(iEmBic,1)});
                mem = (R == max(R,[],2));

                %% BIC
                [bic(ll, 1, iEmBic), like(ll, 1, iEmBic), pen(ll, 1, iEmBic)] = BIC_F(data(:,:,ii_eps,ii_mc), S_est, mu_est, t, mem, rho{em_bic(iEmBic,2)}, psi{em_bic(iEmBic,2)}, eta{em_bic(iEmBic,2)});    
                [bic(ll, 2, iEmBic), like(ll, 2, iEmBic), pen(ll, 2, iEmBic)] = BIC_A(S_est, t, mem, rho{em_bic(iEmBic,2)}, psi{em_bic(iEmBic,2)}, eta{em_bic(iEmBic,2)});                    
                [bic(ll, 3, iEmBic), like(ll, 3, iEmBic), pen(ll, 3, iEmBic)] = BIC_S(S_est, t, mem, rho{em_bic(iEmBic,2)});
            end
        end
        
        bic_final(ii_mc, ii_eps, :, :, :) = bic;
        like_final(ii_mc, ii_eps, :, :, :) = like;
        pen_final(ii_mc, ii_eps, :, :, :) = pen;
    end
    disp(num2str(ii_eps) + "/" + num2str(eps_iter))
    toc
end


%% Evaluation

for iEmBic = 1:embic_iter
    for ii_eps = 1:eps_iter
        for k = 1:size(bic_final, 4)
            BICmax = (permute(bic_final(:,ii_eps,:,k,iEmBic), [1 3 4 2 5]) == max(permute(bic_final(:,ii_eps,:,k,iEmBic), [1 3 4 2 5]),[],2));
% 
%             K_hat = (BICmax .* (1:L_max));
%             K_hat = K_hat(K_hat ~= 0);
%             MAE(k,ii_eps,ii_embic) = sum(abs(K_true - K_hat))/MC;

            K_true_det = repmat([(K_true == 1:max(K_true)) zeros(1, L_max-max(K_true))],MC,1) == 1;
            K_true_under = repmat([~(K_true == 1:max(K_true-1)) zeros(1, L_max-max(K_true-1))],MC,1) == 1;

            p_under(k,ii_eps,iEmBic) = sum(BICmax(K_true_under), 'all')/MC;
            p_det(k,ii_eps,iEmBic) = sum(BICmax(K_true_det))/MC;
            p_over(k,ii_eps,iEmBic) = 1 - p_det(k,ii_eps,iEmBic) - p_under(k,ii_eps,iEmBic);
        end
    end
end

%% plot

g_names = ["Gaus", "t", "Huber", "Tukey"];
names = ["Finite", "Asymptotic", "Schwarz"];
p_det_2 = permute(p_det, [2 1 3]);
[data, labels_true, r, N, K_true, mu_true, S_true] = data_31(N_k, 0);

for iEmBic = 1:embic_iter
    for k_bic = 1:size(bic_final, 4)
        Z = reshape(p_det_2(:,k_bic,iEmBic),size(X));
        fig = figure;
        [M,c] = contour(X,Y,Z);
        c.LineWidth = 1.5;
        hold on
        plot_scatter([labels_true, data], K_true, r)
        title("EM-" + g_names(em_bic(iEmBic,1)) + ", BIC-" + g_names(em_bic(iEmBic,2)) + "-" + names(k_bic))
        colorbar
        caxis([0 1])
  
        % save to .csv
        T = array2table([contour(X,Y,Z).']); 
        writetable(T,"result/sensitivity_EM_" + g_names(em_bic(iEmBic,1)) + "_BIC_" + g_names(em_bic(iEmBic,2))+ "_" + names(k_bic) + "_Nk_" + num2str(N_k) + "_step_" + num2str(step_eps) + "_MC_" + num2str(MC) + "_1.csv", 'Delimiter','tab')
        T = array2table([data, labels_true]); 
        T.Properties.VariableNames = ["x", "y", "label"];
        writetable(T,"result/sensitivity_EM_" + g_names(em_bic(iEmBic,1)) + "_BIC_" + g_names(em_bic(iEmBic,2))+ "_" + names(k_bic) + "_Nk_" + num2str(N_k) + "_step_" + num2str(step_eps) + "_MC_" + num2str(MC) + "_2.csv", 'Delimiter','tab')
    end
end
