%Define a function which computes the Vandermonde matrix, outputted as M
%x is a vector of data points
%n is the degree to which we wish to fit our polynomial to
function [M] = myVandermonde(n,x)
    
    %Assign L to be the number of data points we have
    L = length(x);
    
    %Make a matrix of zeros so the appropriate Vandermonde matrix can be
    %constructed, with the correct size depending on n and L
    M = zeros(L,n+1);

    %Set the first column of M to be a column of 1s
    M(:,1) = 1;
      
    %For the j-th column of M, set it to be x_i^(j-1) using our previous
    %column and multiplying element-wise by the vector x
    for j=2:n+1
        M(:,j) = M(:,j-1).*x;
    end
    
end