function eta = eta_t(t, r, nu)
% Computes eta(t) of the t distribution
%
% Inputs:
%       t  - (N, 1) squared Mahalanobis distance
%       r  - (1, 1) dimension
%       nu - (1, 1) degree of freedom
%
% Outputs:
%       eta - (N, 1) eta(t) of t distribution
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing
    
    eta = - 0.5 .* (nu + r)./(nu + t).^(2);
end