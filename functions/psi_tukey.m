function psi = psi_tukey(t, cT)
% Computes psi(t) of the Tukey loss function
%
% Inputs:
%       t  - (N, 1) squared Mahalanobis distance
%       cT  - (1, 1) tuning parameter
%
% Outputs:
%       psi - (N, 1) psi(t) of Tukey loss function
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    psi = zeros(length(t),1);
    psi(t <= cT^2) = 0.5.*(1 - t(t <= cT^2)./cT^2).^2;
    %psi(t > cT^2) = 0;
end