function [bic, like, pen] = BIC_F(data, S_est, mu_est, t, mem, rho, psi, eta)
% computes the BIC of a RES distribution based on a finite sample penalty term
%
% Inputs:
%        data  - (N, r) data samples
%        S_est - (r, r, ll) estimated Scatter matrix of all clusters
%        mu_est - (r, ll) estimated mean values of all clusters
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
    D = duplicationMatrix(r);
    q = 1/2*r*(r+3);
    
    temp_rho = zeros(1,ll);
    logdetS = zeros(1,ll);
    detJ = zeros(1,ll);

    for m = 1:ll
        x_hat_m = data(mem(:,m), 2:end).' - mu_est(:,m);
        t_m = t(mem(:,m), m);
        J = FIM_RES(x_hat_m, t_m, S_est(:,:,m), psi, eta, D);
        detJ(m) = det(J);
        temp_rho(m) = sum(rho(t(mem(:,m), m)));
        logdetS(m) = log(det(S_est(:,:,m)));
        
        if(detJ(m) < 0)
            warning("negative determinant, J still not positive definite") 
            detJ(m) = detJ(m) + 10^(-10);
            if(detJ(m) < 0)
                detJ(m) = abs(detJ(m));
            end
        elseif(detJ(m) == 0 && N_m(m) == 0)
            warning("cluster without data point, zero determinant")  
            detJ(m) = 1;
        elseif(detJ(m) == 0)
            warning("zero determinant")  
            detJ(m) = detJ(m) + 10^(-10);
            if(detJ(m) < 0)
                detJ(m) = 1;
            end
        end
    end

    like = - sum(temp_rho) + sum(N_m(N_m > 0) .* log(N_m(N_m > 0))) - sum(N_m(N_m > 0) .* logdetS(N_m > 0))/2;
    pen =  - 1/2 * sum(log(detJ)) + ll*q/2*log(2*pi) - ll*log(ll) ;

    bic = like + pen;  
end