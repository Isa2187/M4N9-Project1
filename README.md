# M4N9-Project1
First project of Computational Linear Algebra (M4N9) module taken in 4th year. (Grade = 100%)

All code was done in MATLAB, and further details of the task can be found in the folder Project_Files.

Using a set of data points, this project revolved around generating a continuous function to represent a surface, in this case the depth of a crater on the moon. Using the least squares algorithm, I first fitted polynomials of different degrees to the data points, and saw how this affected the fit. Following this, the moving least squares algorithm was implemented, which is similar in nature to the standard least squares algorithm but considers a weight function that increases the significance of local points to ensure a more accurate fit.

Whilst more accurate, the moving least squares algorithm is more computationally expensive, so further methods were examined to try and increase the accuracy of the fit at a low compuational cost. These included taking a random selection of data points and modelling on this smaller data set, as well as using an 'accelerated' version of moving least squares, which neglected data points that fell into a negligible range.
