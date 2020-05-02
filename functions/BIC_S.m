function [bic, like, pen] = BIC_S(S_est, t, mem, rho)
% computes the BIC of a RES distribution with Schwarz Penalty Term
%
% Inputs:
%        S_est - (r, r, ll) estimated Scatter matrix of cluster m
%        t - (N, ll) squared Mahalanobis distances of data points in cluster m
%        mem - (N, ll) cluster memberships
%        rho - eta of density generator g
%
% Outputs: 
%        bic - (1, 1) bic
%        pen - (1, 1) penalty term
%        like - (1, 1) likelihood term
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    N_m = sum(mem);
    r = size(S_est,1);
    ll = size(S_est,3);
    q = 1/2*r*(r+3);
    N = size(t, 1);
    
    if (N == 0)
       N = 1; 
    end
    
    temp_rho = zeros(1,ll);
    logdetS = zeros(1,ll);

    for m = 1:ll
        temp_rho(m) = sum(rho(t(mem(:,m), m)));
        logdetS(m) = log(det(S_est(:,:,m)));
    end

    like = - sum(temp_rho) + sum(N_m(N_m > 0) .* log(N_m(N_m > 0))) - sum(N_m(N_m > 0) .* logdetS(N_m > 0))/2;
    pen =  - q*ll/2*log(N);

    bic = like + pen;    
end

