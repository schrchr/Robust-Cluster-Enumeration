close all
clear all
addpath("functions", "result");

% This file simulates the probability of detection over the percentage of outliers.
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

%% User Input
% number of threads
%parpool(4); 
% percentage of replacement outliers
epsilon = 0:0.1:0.35;
% number of data points per cluster
N_k = 250; 
% Monte Carlo iterations
MC = 10;
% Select combinations of EM and BIC to be simulated
% 1: Gaussian, 2: t, 3: Huber, 4: Tukey
em_bic = [1 1;
          2 2;
          2 4;
          3 3;
          3 4];
           
% design parameter
% t:
nu = 3;
% Huber:
qH = 0.8;
%Tukey
cT = 4.685;

%% data generation     
embic_iter = length(em_bic);
eps_iter = length(epsilon);

for iEpsilon = 1:eps_iter
    for iMC = 1:MC
        [data(:,:,iEpsilon,iMC), labels_true, r, N, K_true, mu_true, S_true] = data_31(N_k, epsilon(iEpsilon));
    end
end

L_max = 2*K_true; % search range

%% model definitions
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

bic_final = zeros(MC, eps_iter, L_max, 3, embic_iter);
like_final = zeros(MC, eps_iter, L_max, 3, embic_iter);
pen_final = zeros(MC, eps_iter, L_max, 3, embic_iter);
   
%% Cluster Enumeration
tic
for iEpsilon = 1:eps_iter
    parfor iMC = 1:MC
        bic = zeros(L_max, 3, embic_iter);   
        like = zeros(L_max, 3, embic_iter);   
        pen = zeros(L_max, 3, embic_iter); 
        for iEmBic = 1:embic_iter
            for ll = 1:L_max
                %% EM
                [mu_est, S_est, t, R] = EM_RES(data(:,:,iEpsilon,iMC), ll, g{em_bic(iEmBic,1)}, psi{em_bic(iEmBic,1)});
                mem = (R == max(R,[],2));

                %% BIC
                [bic(ll, 1, iEmBic), like(ll, 1, iEmBic), pen(ll, 1, iEmBic)] = BIC_F(data(:,:,iEpsilon,iMC), S_est, mu_est, t, mem, rho{em_bic(iEmBic,2)}, psi{em_bic(iEmBic,2)}, eta{em_bic(iEmBic,2)});    
                [bic(ll, 2, iEmBic), like(ll, 2, iEmBic), pen(ll, 2, iEmBic)] = BIC_A(S_est, t, mem, rho{em_bic(iEmBic,2)}, psi{em_bic(iEmBic,2)}, eta{em_bic(iEmBic,2)});                    
                [bic(ll, 3, iEmBic), like(ll, 3, iEmBic), pen(ll, 3, iEmBic)] = BIC_S(S_est, t, mem, rho{em_bic(iEmBic,2)});
            end
        end
        
        bic_final(iMC, iEpsilon, :, :, :) = bic;
        like_final(iMC, iEpsilon, :, :, :) = like;
        pen_final(iMC, iEpsilon, :, :, :) = pen;
    end
    disp(num2str(epsilon(iEpsilon)))
    toc
end


%% Evaluation
for iEmBic = 1:embic_iter
    for iEpsilon = 1:eps_iter
        for k = 1:size(bic_final, 4)
            BICmax = (permute(bic_final(:,iEpsilon,:,k,iEmBic), [1 3 4 2 5]) == max(permute(bic_final(:,iEpsilon,:,k,iEmBic), [1 3 4 2 5]),[],2));

            K_true_det = repmat([(K_true == 1:max(K_true)) zeros(1, L_max-max(K_true))],MC,1) == 1;
            K_true_under = repmat([~(K_true == 1:max(K_true-1)) zeros(1, L_max-max(K_true-1))],MC,1) == 1;

            p_under(k,iEpsilon,iEmBic) = sum(BICmax(K_true_under), 'all')/MC;
            p_det(k,iEpsilon,iEmBic) = sum(BICmax(K_true_det))/MC;
            p_over(k,iEpsilon,iEmBic) = 1 - p_det(k,iEpsilon,iEmBic) - p_under(k,iEpsilon,iEmBic);
        end
    end
end

%% Plot & Save

marker = {'o','s','d','*','x','^','v','>','<','p','h', '+','o'};
g_names = ["Gaus", "t", "Huber", "Tukey"];
names = ["Finite", "Asymptotic", "Schwarz"];

for iEmBic = 1:embic_iter
    fig = figure;
    h = plot(epsilon, p_det(:,:,iEmBic).', 'LineWidth', 1.5);
    hold on
    grid on
    set(h,{'Marker'}, {marker{1:size(bic_final, 4)}}.')
    xlabel("% of outliers")
    ylabel("probability of detection")
    ylim([0 1])
    legend(names, 'Location', 'southwest')
    title("Nk-" + num2str(N_k) + ", EM-" + g_names(em_bic(iEmBic,1)) + ", BIC-" + g_names(em_bic(iEmBic,2)))

    % save to .csv
    T = array2table([epsilon.', p_det(:,:,iEmBic).']);
    T.Properties.VariableNames = ["x", names];
    writetable(T,"result/outliers_EM_" + g_names(em_bic(iEmBic,1)) + "_BIC_" + g_names(em_bic(iEmBic,2)) + "_MC_" + num2str(MC) + "_Nk_" + num2str(N_k) + ".csv", 'Delimiter','tab')
end


for iEmBic = 1:embic_iter
    names_all(iEmBic,:) = ["EM: " + g_names(em_bic(iEmBic,1)) + ", BIC: " + g_names(em_bic(iEmBic,2)) + "-" + names];
end

fig = figure;
for iEmBic = 1:embic_iter
    h = plot(epsilon, p_det(:,:,iEmBic).', 'LineWidth', 1.5);
    hold on
    set(h,{'Marker'}, {marker{1:size(bic_final, 4)}}.')
end
grid on
xlabel("% of outliers")
ylabel("probability of detection")
ylim([0 1])
legend(names_all.', 'Location', 'southeast')
title("Nk-" + num2str(N_k))

% save to .csv
axis = get(gca,'Children');  
fig_x = axis.XData;
fig_x = fig_x.';
fig_y = flip(cell2mat({axis.YData}.').', 2);
fig_leg = flip(string({axis.DisplayName}));
T = array2table([fig_x, fig_y]);
T.Properties.VariableNames = ["x", fig_leg];
writetable(T,"result/outliers_all_MC_" + num2str(MC) + "_Nk_" + num2str(N_k) + ".csv", 'Delimiter','tab')


for iEmBic = 1:embic_iter
    names_3(iEmBic,:) = ["EM: " + g_names(em_bic(iEmBic,1)) + ", BIC: " + g_names(em_bic(iEmBic,2))];
end

for ii_bic = 1:size(bic_final, 4)% ii_embic = 1:embic_iter
    fig = figure;
    h = plot(epsilon, permute(p_det(ii_bic,:,:),[2 3 1]), 'LineWidth', 1.5);
    hold on
    set(h,{'Marker'}, {marker{1:embic_iter}}.')
    grid on
    xlabel("% of outliers")
    ylabel("probability of detection")
    ylim([0 1])
    legend(names_3, 'Location', 'southwest')
    title("Nk-" + num2str(N_k) + ", BIC-" + names(ii_bic))

    % save to .csv
    T = array2table([epsilon.', permute(p_det(ii_bic,:,:), [2 3 1])]);
    T.Properties.VariableNames = ["x", names_3.'];
    writetable(T,"result/outliers_BIC_" + names(ii_bic) + "_MC_" + num2str(MC) + "_Nk_" + num2str(N_k) + ".csv", 'Delimiter','tab')
end
