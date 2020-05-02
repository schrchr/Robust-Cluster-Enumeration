function t = mahalanobisDistance(x, mu, S)
% Computes the squared Mahalanobis distance.
%
% Inputs:
%       x - (N, r) data, r - dimension
%       mu - (r, 1) cluster centroid
%       S - (r, r) cluster scatter matrix
%
% Outputs:
%       t - (N, 1) squared Mahalanobis distances
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    % make sure, that S is not an ill conditioned matrix
    t = dot(((x.' - mu).' / S).',(x.' - mu), 1).';

%% Alternative 1, slow
%     S_inv = inv(S);
%     N = length(x);
%     t = zeros(N,1);
%     for n = 1:N
%         t(n,1) = (x(n,:).' - mu).' * S_inv * (x(n,:).' - mu);
%     end
 
%% Alternative 2, fast, but not precise
%     S_inv = inv(S);
%     x_centered = (x - mu.').';
%     [L, ~] = chol(S_inv,'lower');
%     t = dot(L*x_centered, L*x_centered, 1).';  
end