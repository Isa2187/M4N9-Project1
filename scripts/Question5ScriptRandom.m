%Note that in each application of the accelerated MLS method below, I have
%fixed epsilon = 10^(-6), n = 3 and eta = 10^(-3)

%Trial 1 of applying the MLS method by randomly taking 5% of the original
%data, storing the residual sum of squares for this solution in 
%sumResidualSquares5a
sumResidualSquares5a = myMLS(10^(-6),10^(-3),3,0.05);

%Trial 2 of applying the MLS method by randomly taking 5% of the original
%data, storing the residual sum of squares for this solution in 
%sumResidualSquares5b
sumResidualSquares5b = myMLS(10^(-6),10^(-3),3,0.05);

%Trial 3 of applying the MLS method by randomly taking 5% of the original
%data, storing the residual sum of squares for this solution in 
%sumResidualSquares5c
sumResidualSquares5c = myMLS(10^(-6),10^(-3),3,0.05);

