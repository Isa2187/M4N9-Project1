%Find the MLS solution using the accelerated MLS method by setting epsilon
%to be 10^(-6), fixing eta = 10^(-3) and letting n = 1, 3 and 10
sumResidualSquares4b = myMLS(10^(-6),10^(-3),[1 3 10],1)