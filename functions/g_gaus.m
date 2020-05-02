function g = g_gaus(t, r)
% Computes g(t) of the Gaussian distribution
%
% Inputs:
%       t  - (N, 1) squared Mahalanobis distance
%       r  - (1, 1) dimension
%
% Outputs:
%       g - (N, 1) g(t) of Gaussian distribution
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    clip = 400; % clipping to avoid zero

    g = (2.*pi).^(-r/2) .* exp(-t./2);
    g(t >= clip) = (2.*pi).^(-r/2) .* exp(-clip./2);
end