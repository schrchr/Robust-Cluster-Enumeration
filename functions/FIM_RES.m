function [J] = FIM_RES(x_hat_m, t_m, S_est, psi, eta, D)
% computes FIM of one cluster for a given RES distribution
%
% Inputs:
%        x_hat_m - (r, N_m) data matrix
%        t_m - (N_m, 1) squared Mahalanobis distances of data points in cluster m
%        S_est - (r, r) estimated Scatter matrix of cluster m
%        psi - psi of density generator g
%        eta - eta of density generator g
%        D - (r^2, 1/2*r*(r+1)) duplication amtrix
%
% Outputs: 
%        J - (q,q) FIM, q = 1/2*r*(r+3)
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    r = size(S_est,1);
    N_m = length(t_m);

    %F_mumu
    temp_eta = zeros(r, r, N_m);
    for n = 1:N_m
        temp_eta(:,:,n) = eta(t_m(n)) .* x_hat_m(:,n) * x_hat_m(:,n).';
    end
    F_mumu = -4.*(S_est\sum(temp_eta, 3))/S_est - S_est \ eye(r,r) .* sum(psi(t_m)) .* 2;

    %F_muS
    temp_eta = zeros(r, r^2, N_m);
    for n = 1:N_m
          temp_eta(:,:,n) = eta(t_m(n)) .* kron((S_est\x_hat_m(:,n) * x_hat_m(:,n).')/S_est, x_hat_m(:,n).'/S_est);
    end
    F_muS = -2 .* sum(temp_eta, 3) * D;

    % F_Smu
    F_Smu = F_muS.';

    temp_eta = zeros(r^2, r^2, N_m);
    for n = 1:N_m
        temp_eta(:,:,n) = eta(t_m(n)) .* kron((S_est\x_hat_m(:,n) * x_hat_m(:,n).')/S_est, (S_est\x_hat_m(:,n) * x_hat_m(:,n).')/S_est);
    end 
    F_SS = - D.' * sum(temp_eta, 3) * D - N_m./2 * D.' * (kron(S_est, S_est) \ eye(r^2,r^2)) * D;
         
    J = [-F_mumu -F_muS; 
         -F_Smu -F_SS];
    
    % there are cases in which the FIM is NOT positive semidefinite. To
    % catch these cases, we use the function nearestSPD(), which calculates
    % the nearest positive semidefinite matrix
    if(det(J) < 0)
        % warning("J not positive definite, use nearestSPD()")         
        J = nearestSPD(J);
    end
end