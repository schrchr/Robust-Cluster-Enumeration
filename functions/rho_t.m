function rho = rho_t(t, r, nu)
% Computes rho(t) of the t distribution
%
% Inputs:
%       t  - (N, 1) squared Mahalanobis distance
%       r  - (1, 1) dimension
%       nu - (1, 1) degree of freedom
%
% Outputs:
%       rho - (N, 1) rho(t) of t distribution
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    rho = -log(gamma((nu + r)/2)) + log(gamma(nu/2)) + 0.5.*r.*log(pi.*nu) + (nu+r)./2.*log(1 + t./nu);
    %rho = -log( gamma((nu + r)/2)/(gamma(nu/2)*(pi*nu).^(r/2)) * (1 + t/nu).^(-(nu+r)/2) );
end