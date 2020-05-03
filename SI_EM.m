close all
clear all
addpath("functions", "result");

% Example for the EM algorithm. Compares the results of a Gaussian, t and
% Huber based EM over different percentages of outliers.

%% User input
% percentage of replacement outliers
epsilon = 0.15;
% Number of samples per cluster
N_k = 50;
% Select degree of freedom for t distribution:
nu = 3;
% Select tuning paramter for Huber distribution:
qH = 0.8;

%% create data   
[data, labels_true, r, N, K_true, mu_true, S_true] = data_31(N_k, epsilon);

%% model definitions
% Huber parameters
cH = sqrt(chi2inv(qH, r));
bH = chi2cdf(cH^2,r+2) + cH^2/r*(1-chi2cdf(cH^2,r));
aH = gamma(r/2)/pi^(r/2) / ( (2*bH)^(r/2)*(gamma(r/2) - igamma(r/2, cH^2/(2*bH))) + (2*bH*cH^r*exp(-cH^2/(2*bH)))/(cH^2-bH*r) );

g = {@(t)g_gaus(t, r);
     @(t)g_t(t, r, nu);
     @(t)g_huber(t, r, cH, bH, aH)};

psi = {@(t)psi_gaus(t);
       @(t)psi_t(t, r, nu);
       @(t)psi_huber(t, r, cH, bH)};

% needed for plots    
x = -20:0.1:20;
y = -20:0.1:20;
[X, Y] = meshgrid(x,y);
g_names = ["Gaussian", "t", "Huber"];

figure
plot_scatter([labels_true data], K_true, r)
for m = 1:K_true
    Z = mvnpdf([X(:) Y(:)],mu_true(:,m).',S_true(:,:,m));
    Z = reshape(Z,size(X));
    contour(X,Y,Z)
end
title("Model: True")
xlabel('Feature 1')
ylabel('Feature 2')

for iModel = 1:length(g)
    % perform EM algorithm
    [mu_est, S_est, t, R] = EM_RES(data, K_true, g{iModel}, psi{iModel});

    figure
    plot_scatter([labels_true data], K_true, r)
    for m = 1:K_true
        Z = mvnpdf([X(:) Y(:)],mu_est(:,m).',S_est(:,:,m));
        Z = reshape(Z,size(X));
        contour(X,Y,Z)
    end
    title("Model: " + g_names(iModel))
    xlabel('Feature 1')
    ylabel('Feature 2')
end

