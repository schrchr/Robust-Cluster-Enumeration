clear all
close all
addpath("functions", "result");

% This file simulates the BIC, likelihood and penalty terms for a given data set.
%
% created by Christian A. Schroth, 3. May 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

%% User input
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
% Tukey
cT = 4.685;
     
%% Your data here
epsilon = 0.15; % percentage of replacement outliers
N_k = 250; % Number of samples per cluster
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

      
embic_iter = length(em_bic);
S_est = cell(L_max,embic_iter);
mu_est = cell(L_max,embic_iter);
for ii_embic = 1:embic_iter
    for ll = 1:L_max
        %% EM
        [mu_est{ll,ii_embic}, S_est{ll,ii_embic}, t, R] = EM_RES(data, ll, g{em_bic(ii_embic,1)}, psi{em_bic(ii_embic,1)});
        mem = (R == max(R,[],2));

        %% BIC
        [bic(ll, 1, ii_embic), like(ll, 1, ii_embic), pen(ll, 1, ii_embic)] = BIC_F(data, S_est{ll,ii_embic}, mu_est{ll,ii_embic}, t, mem, rho{em_bic(ii_embic,2)}, psi{em_bic(ii_embic,2)}, eta{em_bic(ii_embic,2)});    
        [bic(ll, 2, ii_embic), like(ll, 2, ii_embic), pen(ll, 2, ii_embic)] = BIC_A(S_est{ll,ii_embic}, t, mem, rho{em_bic(ii_embic,2)}, psi{em_bic(ii_embic,2)}, eta{em_bic(ii_embic,2)});
        [bic(ll, 3, ii_embic), like(ll, 3, ii_embic), pen(ll, 3, ii_embic)] = BIC_S(S_est{ll,ii_embic}, t, mem, rho{em_bic(ii_embic,2)});                 
    end
end


%% Plots
x = -20:0.1:20;
y = -20:0.1:20;
[X, Y] = meshgrid(x,y);
g_names = ["Gaussian", "t", "Huber"];

marker = {'o','s','d','*','x','^','v','>','<','p','h', '+','o'};
names = ["Finite", "Asymptotic", "Schwarz"];
g_names = ["Gaus", "t", "Huber", "Tukey"];

% BICs
for ii_embic = 1:embic_iter
    figure
    subplot(1,2,1)
    plot_scatter([labels_true data], K_true, r)
    for m = 1:K_true
        Z = mvnpdf([X(:) Y(:)],mu_est{K_true,ii_embic}(:,m).',S_est{K_true,ii_embic}(:,:,m));
        Z = reshape(Z,size(X));
        contour(X,Y,Z)
    end
    title("EM: " + g_names(em_bic(ii_embic,1)) + " at K = " + num2str(K_true))
    xlabel('Feature 1')
    ylabel('Feature 2')
    
    subplot(1,2,2)
    h =  plot(bic(:,:,ii_embic), 'LineWidth', 1.5);
    grid on
    xlabel("number of clusters")
    ylabel("BIC")
    set(h,{'Marker'}, {marker{1:size(bic, 2)}}.')
    legend(names, 'Location', 'southwest')
    title("Nk: " + num2str(N_k) + ", eps: " + num2str(epsilon) + ", EM-" + g_names(em_bic(ii_embic,1)) + ", BIC-" + g_names(em_bic(ii_embic,2)))
end
