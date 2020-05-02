%% The duplication matrix

% For details, see: 
%
% [1] F. K. Teklehaymanot, M. Muma, and A. M. Zoubir, "Bayesian Cluster Enumeration Criterion for Unsupervised Learning",
%     IEEE Trans. Signal Process. (accepted),
%     [Online-Edition: https://arxiv.org/abs/1710.07954v2], 2018.
%
% [2] F. K. Teklehaymanot, M. Muma, and A. M. Zoubir, "Novel Bayesian Cluster
%     Enumeration Criterion for Cluster Analysis With Finite Sample Penalty Term", 
%     in Proc. 43rd IEEE Int. conf. on Acoustics, Speech and Signal Process. (ICASSP), pp. 4274-4278, 2018, 
%     [Online-Edition: https://www.researchgate.net/publication/322918028]


% Copyright (c) 2018 Freweyni K. Teklehaymanot. All rights reserved.


% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, version 3 of the License.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Input:
%        num_features - number of features in the data set


% Output:
%        dup_mat - the duplication matrix


function dup_mat = duplicationMatrix(num_features)

dup_mat_transpose = zeros(.5*num_features*(num_features+1),num_features^2); % transpose of the duplication matrix

for j = 1:num_features
    
    for i = 1:num_features
        
        u = zeros(.5*num_features*(num_features+1),1);
        idx = (j-1)*num_features+i-.5*j*(j-1);
        u(idx,1) = 1;
        
        Y = zeros(num_features,num_features);
        Y(i,j) = 1;
        Y(j,i) = 1;
        d = u*reshape(Y,[],1)';
        
        if i>=j
            dup_mat_transpose = dup_mat_transpose + d;
        end
        
    end
    
end

dup_mat = dup_mat_transpose'; 

end
