function eta = eta_tukey(t, cT)
% Computes eta(t) of the Tukey loss function
%
% Inputs:
%       t  - (N, 1) squared Mahalanobis distance
%       cT  - (1, 1) tuning parameter
%
% Outputs:
%       eta - (N, 1) eta(t) of Tukey loss function
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    eta = zeros(length(t),1);
    eta(t <= cT^2) = t(t <= cT.^2)./cT.^4  - 1./cT.^2;
    %eta(t > cT^2) = 0;
end