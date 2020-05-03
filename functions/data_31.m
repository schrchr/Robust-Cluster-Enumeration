% close all
% clear all
% 
% epsilon = 0.2;
% N_k = 200;

function [data, labels, r, N, K_true, mu_true, scatter_true] = data_31(N_k, epsilon)
% Creates three Gaussian clusters
% based on http://arxiv.org/pdf/1811.12337v1
%
% Inputs: 
%        N_k - (1,1) number of data vectors per cluster
%        epsilon - (1,1) percentage of replacement outlieres
%
% Outputs:
%        data  - (N, r+1) data
%        labels - (N,1) true labels
%        r - (1,1) number of features/dimensions in the generated data set
%        N - (1,1) total number of samples in the data set
%        K_true - (1,1) true number of clusters in the data set
%        mu_true - (r, K_true) true mean values
%        scatter_true - (r, r, K_true) true scatter matrices
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    out_range = [-20 20; -20 20];
    K_true = 3; % number of clusters
    r = 2; % number of features/dimension

    mu_true = zeros(r, K_true);
    mu_true(:,1) = [0; 5];
    mu_true(:,2) = [5; 0];
    mu_true(:,3) = [-5; 0];

    scatter_true = zeros(r, r, K_true);
    scatter_true(:,:,1) = [2 0.5; 0.5 0.5];
    scatter_true(:,:,2) = [1 0; 0 0.1];
    scatter_true(:,:,3) = [2 -0.5; -0.5 0.5];


    N = K_true*N_k; % total number of data points

    data = [];
    for k = 1:K_true
        data = [data ;[ones(N_k,1)*k, mvnrnd(mu_true(:,k), scatter_true(:,:,k),N_k)]];
    end

    % randomly permute data
    data = data(randperm(N), :);

    % replacement outlier
    N_repl = round(N * epsilon);
    index_repl = randperm(N,N_repl).';

    data_rpl = [];
    for ir = 1:r
        data_rpl = [data_rpl, rand(N_repl,1)*(out_range(ir,2) - out_range(ir,1)) + out_range(ir,1)];
    end
    data(index_repl,:) = [ones(N_repl,1)*(K_true+1), data_rpl];
    
    labels = data(:,1);
    data = data(:,2:end);
end