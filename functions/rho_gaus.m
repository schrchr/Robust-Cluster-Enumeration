function rho = rho_gaus(t, r)
% Computes rho(t) of the Gaussian distribution
%
% Inputs:
%       t  - (N, 1) squared Mahalanobis distance
%       r  - (1, 1) dimension
%
% Outputs:
%       rho - (N, 1) rho(t) of Gaussian distribution
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    rho = r./2.*log(2.*pi) + t./2;
    %rho = -log((2.*pi).^(-r./2) .* exp(-t./2)); 
end