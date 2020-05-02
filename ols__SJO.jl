
"""
OLS estimation procedure adapted from Bradley Setzler
https://juliaeconomics.com/2014/06/15/introductory-example-ordinary-least-squares/
2020-05-03
"""

using Distributions, LinearAlgebra, Random

Random.seed!(2);

#PRELIMINARIES-----------------------------------------------------------------*
N = 1000;           # size of sample
K = 3;              # number of variables

#GENERATE INDEPENDENT  AND DEPENDENT VARIABLES---------------------------------*
genX = MvNormal(Matrix(1.0I,K,K))
X = rand(genX,N)'

X_noconstant = X;
constant = ones(N);
X = hcat(constant,X_noconstant);

genEpsilon = Normal(0,1);
epsilon = rand(genEpsilon,N);

trueParams = [0.1, 0.5, -0.3, 0.0];

Y = X*trueParams + epsilon

#RECOVER THE TRUE COEFFICIENTS-------------------------------------------------*

# explore the linear form of the equation and the exogeneity assumpton
#to obtain coefficents
# the inverse of matrix multiplied by the matrix gives you the identity matrix
#Y = X*b + ϵ

## start by multiplying through with an inverse
b = X\Y;                 #(X\X) is an identiy
b = (X\X)\(X\Y);         # full specifiction with inclusion of identity

## start by multiplying through with a transpose
b = (X'X)\X'*Y;
b = inv(X'X)*X'*Y

#END---------------------------------------------------------------------------*
