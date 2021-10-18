function [data_final, label_int, r, N, K_true] = data_radar()

    load Features_PCA_ABCD.mat;
    
    labels = string(PCA_Strides.Persons);
    N = length(labels);
    label_int = zeros(N,1);
    label_int(labels == "A") = 1;
    label_int(labels == "B") = 2;
    label_int(labels == "C") = 3;
    label_int(labels == "D") = 4;
    
    data = PCA_Strides.Features; 
    
    %data = data ./ max(abs(data));  
    data = data ./ mean(data);
    
    data_final = [label_int data]; 

    K_true = 4; 
    r = 5;
end