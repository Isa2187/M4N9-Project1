%Note that in each application of the accelerated MLS method below, I have
%fixed epsilon = 10^(-6), n = 3 and eta = 10^(-3)

%Reduce the cost of the MLS method by taking 1 in every 2 points, storing 
%the residual sum of squares for this solution in sumResidualSquares5i
sumResidualSquares5i = myMLS(10^(-6),10^(-3),3,2);

%Reduce the cost of the MLS method by taking 1 in every 3 points, storing 
%the residual sum of squares for this solution in sumResidualSquares5ii
sumResidualSquares5ii = myMLS(10^(-6),10^(-3),3,3);

%Reduce the cost of the MLS method by taking 1 in every 4 points, storing 
%the residual sum of squares for this solution in sumResidualSquares5iii
sumResidualSquares5iii = myMLS(10^(-6),10^(-3),3,4);

%Reduce the cost of the MLS method by taking 1 in every 5 points, storing 
%the residual sum of squares for this solution in sumResidualSquares5iv
sumResidualSquares5iv = myMLS(10^(-6),10^(-3),3,5);

%Reduce the cost of the MLS method by taking 1 in every 10 points, storing 
%the residual sum of squares for this solution in sumResidualSquares5v
sumResidualSquares5v = myMLS(10^(-6),10^(-3),3,10);

%Reduce the cost of the MLS method by taking 1 in every 20 points, storing 
%the residual sum of squares for this solution in sumResidualSquares5vi
sumResidualSquares5vi = myMLS(10^(-6),10^(-3),3,20);

%Reduce the cost of the MLS method by taking 1 in every 30 points, storing 
%the residual sum of squares for this solution in sumResidualSquares5vii
sumResidualSquares5vii = myMLS(10^(-6),10^(-3),3,30);