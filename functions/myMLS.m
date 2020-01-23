%Define a function which computes the MLS solution for each value of n and 
%each value of eta stored in nValues and etaValues respectively, and
%returns a matrix of the residual sum of squares for each combination of
%values of n and eta
%epsilon is a small real number that is used to apply the accelerated MLS
%nValues is a vector consisting of integer values of n 
%etaValues is a vector consisting of real values of the parameter eta
%xDistn is a number that determines in what way we should sample the 
%original data points - note that if xDistn = 1 then we use all the data 
function[sumResidualSquares] = myMLS(epsilon,etaValues,nValues,xDistn)

    %Load the data for the crater so it can be used in the code
    load('h.mat')
    load('x.mat')  
    
    %Make copies of the original data so they are not overriden if we take
    %samples of it (depending on the value of xDistn, we may or may not
    %take a sample)
    xData = x;
    hData = h;
    
    %If xDistn > 1 then we only take every xDistn-th point to reduce the
    %number of data points used, reducing the cost of the MLS method
    if xDistn > 1
        x = x(1:xDistn:end);
        h = h(1:xDistn:end);
    end
    
    %If xDistn < 1, we take a random sample of this percentage of points
    %from the original data, again to reduce the cost of the MLS method
    if xDistn < 1
        sampleSize = round(xDistn*length(x));
        indices = randperm(length(x));
        indices = indices(1:sampleSize);
        x = x(indices);
        h = h(indices);
    end
    
    %Initialise a matrix which will store the residual sum of squares for
    %each case, where each row represents the same value of n and
    %each column represents the same value of eta
    sumResidualSquares = zeros(length(nValues),length(etaValues));

    %Go through each value of n
    for p = 1:length(nValues)
        %Go through each value of eta
        for q = 1:length(etaValues)

            %Initialise a vector of points called xValues, and at each of
            %these points we will compute the MLS solution at that point 
            %and store the value of the best-fitted polynomial in hValues
            %Note that I made xValues match the same points as in the
            %original data so comparisons are point-wise easier to make
            xValues =  -0.5:0.001:0.5;
            hValues = zeros(length(xValues),1);

            %Set n and eta so they can be used in the MLS method
            n = nValues(p);
            eta = etaValues(q);

            %Initialise our weight function so it can be called upon below
            weightFunction = @(r) exp(-(r^2)/eta);

            %Compute the appropriate Vandermonde matrix
            A = myVandermonde(n,x);

            %Apply the MLS method at each point in xValues
            for i=1:length(xValues)

                %Compute the weight values as per the instructions in the
                %algorithm
                weightValues = zeros(length(x),1);
                xValue = xValues(i);                
                for j=1:length(x)
                    weightValues(j) = weightFunction(abs(x(j) - xValue));

                    %If epsilon ~= 0, then set any weight values that are
                    %less than epsilon to 0
                    if epsilon ~= 0
                        if weightValues(j) < epsilon
                            weightValues(j) = 0;
                        end
                    end
                end

                %To implement MLS, multiply each row of A and h by the
                %respective weight values as described in Question 3
                AMLS = weightValues .* A;
                HMLS = weightValues .* h;

                %To use the accelerated MLS method, if we have any all-zero 
                %rows in AMLS and HMLS then remove them so we have a
                %smaller system of linear equations
                if epsilon ~= 0
                    AMLS(all(~AMLS,2),:) = [];
                    HMLS(all(~HMLS,2),:) = [];
                end

                %Use the provided function house with A as its input to get
                %the matrix R in the QR decomposition, as well as the 
                %matrix H which stores the vectors to be used in the 
                %Householder reflections
                [H,R] = house(AMLS);    

                %Set N to be the number of data points that we have in our 
                %MLS problem, and initialise b to be HLMS, from 
                %which we will then override the values of b
                N = size(AMLS,1);    
                b = HMLS;

                %Override the values of b according to the appropriate 
                %Householder reflection vectors which will then be used to
                %solve our MLS problem and find the best polynomial fit of
                %degree n
                for k=1:(n+1)
                    b(k:N) = b(k:N) - (2.0)*(H(k:N,k)'*b(k:N))*H(k:N,k);        
                end

                %Since we want the reduced QR factorisation of R, set Rhat 
                %to be the square matirx by taking only the top part of R 
                %by omitting the chunk of zeros beneath Rhat in R, so that 
                %Rhat is and (n+1)x(n+1) upper triangular matrix
                Rhat = R(1:n+1,:);

                %Initialise a vector coeff of n+1 zeros, which will store 
                %the coefficients of the fitted polynomial of degree n
                coeff = zeros(n+1,1);

                %Use backwards substitution to find the coefficients of our
                %polynomial,using the fact that Rhat is an nxn upper 
                %triangular matrix
                coeff(n+1) = b(n+1)/Rhat(n+1,n+1);
                for j=n:-1:1
                    coeff(j) = (b(j) - dot(Rhat(j,:),coeff') )/Rhat(j,j);
                end

                %Using the least squares polynomial found, evaluate the
                %ploynomial at the point xValue to obtain the MLS solution 
                %at this point, and store it in hValues
                xPowers = zeros(1,n+1);
                for j=1:n+1
                    xPowers(j) = xValue^(j-1);
                end                
                hValues(i) = xPowers * coeff;

            end

            %Compute the residual sum of squares by getting the norm of the 
            %vector whose entries are the difference of the fitted
            %polynomial and the original data at each data point
            sumResidualSquares(p,q) = norm(hValues - hData);

            figure();
            %Plot the original data
            plot(xData,hData);
            hold on;
            %Plot the solutiuon given by MLS
            plot(xValues,hValues);
            hold off;
            title(sprintf('Graph Of MLS Solution h Against Original Data h_{i} With n = %s And eta = %s',num2str(n),num2str(eta)))
            legend('Original data','MLS solution')
            xlabel('Point on diameter relative to centre')
            ylabel('Height of crater')
        
        end        
    end
    
end