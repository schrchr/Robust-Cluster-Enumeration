function plot_scatter(data, K_true, r)
% Plots data with given labels for r = 1,2,3
%
% Inputs:
%       data  - (N, r+1) data(:,1) includes labels from 1 to K_true+1,
%               where K_true+1 are the outliers, data(:,2:r+1) includes the actual data points
%       K_true  - (1, 1) true number of clusters
%       r  - (1, 1) dimension
%
% created by Christian A. Schroth, 30. April 2020
%
% "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
% Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
% submitted to IEEE Transactions on Signal Processing

    for k = 1:K_true+1
         if(r == 1)
            if(k == K_true+1)
                scatter(data(data(:,1)==k, 2),zeros(sum(data(:,1)==k), 1),'k', 'o')
            else
                scatter(data(data(:,1)==k, 2),zeros(sum(data(:,1)==k), 1),'filled')
                hold on
            end
         elseif(r == 2)
            if(k == K_true+1)
                scatter(data(data(:,1)==k, 2),data(data(:,1)==k, 3),15,'k', 'o')
            else
                scatter(data(data(:,1)==k, 2),data(data(:,1)==k, 3),15,'filled')
                hold on
            end
        elseif(r == 3)
            if(k == K_true+1)
                scatter3(data(data(:,1)==k, 2),data(data(:,1)==k, 3),data(data(:,1)==k, 4),'k', 'o')
            else
                scatter3(data(data(:,1)==k, 2),data(data(:,1)==k, 3),data(data(:,1)==k, 4),'filled')
                hold on
            end
        end
    end
    xlabel('Feature 1')
    ylabel('Feature 2')
    zlabel('Feature 3')
end