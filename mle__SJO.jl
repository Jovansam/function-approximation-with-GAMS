
"""
MLE estimation procedure adapted from Bradley Setzler
https://juliaeconomics.com/2014/06/16/numerical-maximum-likelihood-the-ols-example/
2020-05-03
"""

using Distributions, Random, Optim, LinearAlgebra;

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

Y = X*trueParams + epsilon;

#RECOVER THE TRUE COEFFICIENTS BY MLE------------------------------------------*

function normal_loglikelihood(ρ)
    β = ρ[1:4]      # four regressions parameters to estimate including constant
    σ = ρ[5]        # but we also need to estimate the standard deviation
    ϵ = Y - X*β     # the residuals from the estimation

    dist_type = Normal(0,sqrt(σ^2))    # specify properties of distribution

    # return the logarithm of a Normal probability density at ϵ[i]
    # you have to provide the type of distribution. That way you do not have
    # to code the likelihood as for instance you would in GAMS.
    contributions = logpdf.(dist_type,ϵ)      # equivalent to: logpdf() = log(pdf())

    loglikelihood = sum(contributions)  # sum over the logged residuals

    return -loglikelihood       # we return a negative because the optimizing function minimizes
end

# We shall use Optim for optimization.
# You could also use JuMP, Convex.jl, or even NLopt
params0 = [.1,.2,.3,.4,.5]
optimum = optimize(normal_loglikelihood,params0,ConjugateGradient())
loglik_MLE = optimum.minimum;
params_MLE = optimum.minimizer;

println("The loglikelihood is: $loglik_MLE")
println("The estimated parameters are: $params_MLE")
println("The true parameters are: $trueParams")

#END---------------------------------------------------------------------------*
