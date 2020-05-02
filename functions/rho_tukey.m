function rho = rho_tukey(t, r, cT)
% Computes rho(t) of the Tukey loss function
%
% Inputs:
%       t  - (N, 1) squared Mahalanobis distance
%       r  - (1, 1) dimension
%       cT  - (1, 1) tuning parameter
%
% Outputs:
%       rho - (N, 1) rho(t) of Tukey loss function     
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    rho = zeros(length(t),1);
    rho(t <= cT^2) = r./2.*log(2.*pi) + t(t <= cT^2).^3 ./(6.*cT.^4) - t(t <= cT.^2).^2 ./ (2.*cT.^2) + t(t <= cT.^2)./2;
    rho(t > cT^2) = r./2.*log(2.*pi) + cT.^2./6;
end
