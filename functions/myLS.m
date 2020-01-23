%Define a function which computes the LS solution for each value of n, and
%returns a vector of the residual sum of squares for each value of n
%nValues is a vector consisting of integer values of n to find the LS
%solutions for
function[sumResidualSquares] = myLS(nValues)

    %Load the data for the crater so it can be used in the code
    load('h.mat')
    load('x.mat')   

    %Initialise a vector which will store the residual sum of squares for
    %each polynomial fit, to then be later updated and outputted by the
    %function
    sumResidualSquares = zeros(1,length(nValues));

    %Go through each value of n
    for i=1:length(nValues)

        %Assign n to be the degree for which we are fitting our polynomial
        n = nValues(i);

        %Use my function myVandermonde to obtain the Vandermonde matrix for
        %this degree n, and call this matrix A
        A = myVandermonde(n,x);

        %Use the provided function house with A as its input to get the 
        %matrix R in the QR decomposition, as well as the matrix H which 
        %stores the vectors to be used for the Householder reflections
        [H,R] = house(A);
        
        %Set N to be the number of data points that we have in our LS
        %problem, and initialise b to be h, from which we will then 
        %override the values of b
        N = size(A,1);    
        b = h;
        
        %Override the values of b according to the appropriate Householder
        %reflection vectors which will then be used to solve our least 
        %squares problem and find the best polynomial fit of degree n
        for k=1:(n+1)
            b(k:N) = b(k:N) - (2.0)*(H(k:N,k)'*b(k:N))*H(k:N,k);        
        end

        %Since we want the reduced QR factorisation of R, set Rhat to be 
        %a square matirx by taking only the top part of R by omitting the
        %chunk of zeros beneath Rhat in R, so that Rhat is and (n+1)x(n+1)
        %upper triangular matrix
        Rhat = R(1:n+1,:);

        %Initialise a vector coeff of n+1 zeros, which will store the
        %coefficients of the fitted polynomial of degree n in ascending
        %order of the powers of x
        coeff = zeros(n+1,1);

        %Use backwards substitution to find the coefficients of our 
        %polynomial, using the fact that Rhat is an (n+1)x(n+1) upper
        %triangular matrix
        coeff(n+1) = b(n+1)/Rhat(n+1,n+1);
        for j=n:-1:1
            coeff(j) = (b(j) - dot(Rhat(j,:),coeff') )/Rhat(j,j);
        end

        %Compute the residual sum of squares by getting the norm of the 
        %vector whose entries are the difference of the fitted polynomial
        %and the original data for each data point
        sumResidualSquares(i) = norm(A*coeff - h);

        figure() 
        %Plot the original data
        plot(x,h)
        hold on
        %Plot the fitted polynomial using our Vandermonde matrix and the
        %vector of coefficients found
        plot(x,A*coeff)
        hold off
        title(sprintf('Graph Of h_{%s} Against Original Data h_{i}',num2str(n)))
        legend('Original data','Fitted polynomial')
        xlabel('Point on diameter relative to centre')
        ylabel('Height of crater')
        
        %Display the 6 monomial coefficients of lowest degree for each LS
        %solution
        disp('Monomial coefficients of 6 lowest degree terms are:')
        disp(coeff(1:6))
        
    end
    
end