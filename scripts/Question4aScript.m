%Find the MLS solution using the accelerated MLS method by setting epsilon
%to be 10^(-6), fixing n = 3 and letting eta = 10^(-2), 10^(-3) and 10^(-4)
sumResidualSquares4a = myMLS(10^(-6),[10^(-2),10^(-3),10^(-4)],3,1)