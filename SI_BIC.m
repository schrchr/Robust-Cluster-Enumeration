clear all
close all
addpath("functions", "result");

% This file simulates the BIC, likelihood and penalty terms for a given data set.
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

%% User input
% number of Monte Carlo iterations
MC = 5; 
% percentage of replacement outliers
epsilon = 0.04;
% Number of samples per cluster
N_k = 100;
% Select combinations of EM and BIC to be simulated
% 1: Gaussian, 2: t, 3: Huber, 4: Tukey
%        EM BIC
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
     
%% Cluster Enumeration      
tic
embic_iter = length(em_bic);
for iMC = 1:MC
    [data, labels_true, r, N, K_true, mu_true, S_true] = data_31(N_k, epsilon);
    L_max = 2*K_true; % search range

    %% model definitions
    % Huber parameters
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

    for ii_embic = 1:embic_iter
        for ll = 1:L_max
            %% EM
            [mu_est, S_est, t, R] = EM_RES(data, ll, g{em_bic(ii_embic,1)}, psi{em_bic(ii_embic,1)});
            mem = (R == max(R,[],2));

            %% BIC
            [bic(iMC, ll, 1, ii_embic), like(iMC, ll, 1, ii_embic), pen(iMC, ll, 1, ii_embic)] = BIC_F(data, S_est, mu_est, t, mem, rho{em_bic(ii_embic,2)}, psi{em_bic(ii_embic,2)}, eta{em_bic(ii_embic,2)});    
            [bic(iMC, ll, 2, ii_embic), like(iMC, ll, 2, ii_embic), pen(iMC, ll, 2, ii_embic)] = BIC_A(S_est, t, mem, rho{em_bic(ii_embic,2)}, psi{em_bic(ii_embic,2)}, eta{em_bic(ii_embic,2)});
            [bic(iMC, ll, 3, ii_embic), like(iMC, ll, 3, ii_embic), pen(iMC, ll, 3, ii_embic)] = BIC_S(S_est, t, mem, rho{em_bic(ii_embic,2)});                 
        end
    end
    toc
    disp(num2str(iMC))
end

%% Averaging over MC
bic_avg = mean(permute(bic, [2 3 4 1]), 4);
like_avg = mean(permute(like, [2 3 4 1]),4);
pen_avg = mean(permute(pen, [2 3 4 1]), 4);

%% Plots

figure
plot_scatter([labels_true data], K_true, r)

marker = {'o','s','d','*','x','^','v','>','<','p','h', '+','o'};
names = ["Finite", "Asymptotic", "Schwarz"];
g_names = ["Gaus", "t", "Huber", "Tukey"];

% BICs
for ii_embic = 1:embic_iter
    figure
    h =  plot(bic_avg(:,:,ii_embic), 'LineWidth', 1.5);
    grid on
    xlabel("number of clusters")
    ylabel("BIC")
    set(h,{'Marker'}, {marker{1:size(bic, 3)}}.')
    legend(names, 'Location', 'southwest')
    title("Nk: " + num2str(N_k) + ", eps: " + num2str(epsilon) + ", EM-" + g_names(em_bic(ii_embic,1)) + ", BIC-" + g_names(em_bic(ii_embic,2)))

    % save to .csv
    T = array2table([[1:L_max].', bic_avg(:,:,ii_embic)]);
    T.Properties.VariableNames = ["x", names];
    writetable(T,"result/bic_Nk_" + num2str(N_k)  + "_eps_" + num2str(round(epsilon*100, 0)) + "_EM_" + g_names(em_bic(ii_embic,1)) + "_BIC_" + g_names(em_bic(ii_embic,2))+ ".csv", 'Delimiter','tab')
end

% Likelihood terms
fig = figure;
for ii_embic = 1:embic_iter
    h =  plot(like_avg(:,1,ii_embic), 'LineWidth', 1.5);
    hold on
    set(h,{'Marker'}, {marker{ii_embic}}.')
    leg_names(ii_embic,:) = [ "EM: " + g_names(em_bic(ii_embic,1)) + ", BIC: " + g_names(em_bic(ii_embic,2))];
end
grid on
xlabel("number of clusters")
ylabel("Likelihood")
legend(leg_names, 'Location', 'southeast')
title("Nk: " + num2str(N_k) + ", eps: " + num2str(epsilon))

% save to .csv
T = array2table([[1:L_max].', permute(like_avg(:,1,:), [1 3 2])]);
T.Properties.VariableNames = ["x", leg_names.'];
writetable(T,"result/like_Nk_" + num2str(N_k) + "_eps_" + num2str(round(epsilon*100, 0)) + ".csv", 'Delimiter','tab')


% Penalty terms
for ii_embic = 1:embic_iter
    fig = figure;
    h =  plot(pen_avg(:,:,ii_embic), 'LineWidth', 1.5);
    grid on
    xlabel("number of clusters")
    ylabel("Penalty")
    set(h,{'Marker'}, {marker{1:size(pen, 3)}}.')
    legend(names, 'Location', 'southwest')
    title("Nk: " + num2str(N_k) + ", eps: " + num2str(epsilon) + ", EM-" + g_names(em_bic(ii_embic,1)) + ", BIC-" + g_names(em_bic(ii_embic,2)))

    % save to .csv
    T = array2table([[1:L_max].', pen_avg(:,:,ii_embic)]);
    T.Properties.VariableNames = ["x", names];
    writetable(T,"result/pen_Nk_" + num2str(N_k) + "_eps_" + num2str(round(epsilon*100, 0)) + "_EM_" + g_names(em_bic(ii_embic,1)) + "_BIC_" + g_names(em_bic(ii_embic,2)) + ".csv", 'Delimiter','tab')
end