function [bic, like, pen] = BIC_A(S_est, t, mem, rho, psi, eta)
% computes the BIC of a RES distribution with the asymptotic penalty term
%
% Inputs:
%        S_est - (r, r, ll) estimated Scatter matrix of cluster m
%        t - (N, ll) squared Mahalanobis distances of data points in cluster m
%        mem - (N, ll) cluster memberships
%        rho - rho of density generator g
%        psi - psi of density generator g
%        eta - eta of density generator g
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
    
    temp_rho = zeros(1,ll);
    temp_psi = zeros(1,ll);
    temp_eta = zeros(1,ll);
    logdetS = zeros(1,ll);
    epsilon = zeros(1,ll);
    
    for m = 1:ll
        temp_rho(m) = sum(rho(t(mem(:,m), m)));
        temp_psi(m) = sum(psi(t(mem(:,m), m)));
        temp_eta(m) = sum(eta(t(mem(:,m), m)));
        
        epsilon(m) = max([abs(temp_psi(m)), abs(temp_eta(m)), N_m(m)]);
        
        logdetS(m) = log(det(S_est(:,:,m)));
    end

    like = - sum(temp_rho(temp_rho > 0)) + sum(N_m(N_m > 0) .* log(N_m(N_m > 0))) - sum(N_m .* logdetS)/2;
    pen =  - 1/2 * q *sum(log(epsilon(epsilon > 0)));

    bic = like + pen;
end