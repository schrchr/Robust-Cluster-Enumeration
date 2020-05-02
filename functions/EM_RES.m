function [mu_hat, S_hat, t, R] = EM_RES(data, ll, g, psi)
% EM algorithm for mixture of RES distributions defined by g, psi
%
% Inputs:
%        data - (N, r) data matrix (without lables!)
%        ll - (1,1) number of clusters
%        g - anonymous function of density generator, e.g. @(t)g_gaus(t, r)
%        psi - anonymous function of psi, e.g. @(t)psi_gaus(t)
%
% Outputs: 
%        mu_hat - (r, ll) final estimate of cluster centroids
%        S_hat - (r, r, ll) final estimate of cluster scatter matrices
%        t - (N, ll) Mahalanobis distances to all clusters
%        R - (N, ll) estimates of the posterior probabilities per cluster
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing
% 
% with parts from:
%
% [1] F. K. Teklehaymanot, M. Muma, and A. M. Zoubir, "Bayesian Cluster Enumeration Criterion for Unsupervised Learning",
%     IEEE Trans. Signal Process. (accepted),
%     [Online-Edition: https://arxiv.org/abs/1710.07954v2], 2018.
%
% [2] F. K. Teklehaymanot, M. Muma, and A. M. Zoubir, "Novel Bayesian Cluster
%     Enumeration Criterion for Cluster Analysis With Finite Sample Penalty Term", 
%     in Proc. 43rd IEEE Int. conf. on Acoustics, Speech and Signal Process. (ICASSP), pp. 4274-4278, 2018, 
%     [Online-Edition: https://www.researchgate.net/publication/322918028]

limit = 1e-6; % a value that determines when the EM algorithm should terminate
em_max_iter = 200; % maximum number of iterations of the EM algorithm
reg_value = 1e-6; % regularization value used to regularize the covariance matrix in the EM algorithm

%% variable initializations
r = size(data, 2);
N = size(data, 1);

v = zeros(N,ll);
v_diff = zeros(N,ll);
tau = zeros(1,ll);
S_hat = zeros(r,r,ll);
t = zeros(N,ll);
log_likelihood = zeros(1, em_max_iter);

%% initialization using K-means++ 
warning('off', "stats:kmeans:FailedToConvergeRep")
%[clu_memb_kmeans, mu_Kmeans] = kmeans(data, ll, 'MaxIter', 10, 'Replicates', 5);
% cityblock distance equals a kmedian
[clu_memb_kmeans, mu_Kmeans] = kmeans(data, ll, 'Distance', 'cityblock', 'MaxIter', 10, 'Replicates', 5);
mu_hat = mu_Kmeans.';
for m = 1:ll
    tau(m) = sum(clu_memb_kmeans == m)/N;
end

for m = 1:ll
    x_hat = data(clu_memb_kmeans == m).' - mu_hat(:,m);
    N_m = sum(clu_memb_kmeans == m);
    S_hat(:,:,m) = (x_hat * x_hat.') ./ N_m;

    % Check if the sample covariance matrix is positive definite
    [~, indicator] = chol(S_hat(:,:,m));
    if (indicator ~= 0 || cond(S_hat(:,:,m)) > 30)
        S_hat(:,:,m) = 1/(r*N_m)*sum(diag(x_hat'*x_hat))*eye(r); % diagonal covariance matrix whose diagonal entries are identical
        [~,indicator] = chol(S_hat(:,:,m));
        %disp("Initial Covariance not PSD")
        if indicator ~= 0
            S_hat(:,:,m) = eye(r); % if the estimated covariance matrix is still singular, then set it to identity
        end
    end
    t(:,m) = mahalanobisDistance(data, mu_hat(:,m), S_hat(:,:,m));
end

%% EM algorithm
for ii = 1:em_max_iter
    % E-step
    v_lower = zeros(N,ll);
    for j = 1:ll
        v_lower(:,j) = tau(j) * det(S_hat(:,:,j))^(-1/2) * g(t(:,j));
    end
    for m = 1:ll
        v(:,m) = tau(m) .* det(S_hat(:,:,m))^(-1/2) .* g(t(:,m)) ./ sum(v_lower,2);
        v_diff(:,m) = v(:,m) .* psi(t(:,m));
    end

    % M-step
    for m = 1:ll
        mu_hat(:,m) = sum(v_diff(:,m) .* data) / sum(v_diff(:,m));
        S_hat(:,:,m) = 2 .* (v_diff(:,m).' .* (data.' - mu_hat(:,m)) * (data.' - mu_hat(:,m)).') ./ sum(v(:,m)) + reg_value*eye(r);
        tau(m) = sum(v(:,m))/N;   
        t(:,m) = mahalanobisDistance(data, mu_hat(:,m), S_hat(:,:,m));
    end

    % convergence check
    v_conv = zeros(N, ll);
    for m = 1:ll
        v_conv(:,m) = tau(m) .* det(S_hat(:,:,m)).^(-1/2) .* g(t(:,m));
    end
    log_likelihood(ii) = sum(log(sum(v_conv.')));

    if(ii > 1)
        if(abs(log_likelihood(ii)-log_likelihood(ii-1)) < limit)
            % break if difference between two consecutive steps is small
            break;
        end   
    end
end
% calculate posterior probabilities
R = v_conv./sum(v_conv, 2);

%% diagonal loading
% If the estimated "Matrix is close to singular or badly scaled" it cannot be inverted. 
% https://math.stackexchange.com/questions/261295/to-invert-a-matrix-condition-number-should-be-less-than-what
% If S_hat has a large condition number, a small number in comparison to the
% matrix entries is added. This step should be subject for further tweaking.
for m = 1:ll
    cond_S = cond(S_hat(:,:,m));
    if(cond_S > 30)
        %warning("S with large condition number")
        S_hat(:,:,m) = S_hat(:,:,m) + 0.01 * 10^floor(log10(trace(S_hat(:,:,m)))) * log10(cond_S) * eye(r);
    end
end

%% plots log-likelihood over iterations for debug 
% figure
% plot(log_likelihood(log_likelihood ~= 0))
% grid on
% xlabel("number of iterations")
% ylabel("log-likelihood")