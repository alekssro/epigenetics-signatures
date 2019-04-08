######################################################
##	MAJORIZE MINIMIZE ALGORITHM FOR THE NQP PROBLEM ##
######################################################
##
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
##
## The algorithm minimizes the quadratic form 
## f(x) = 1/2x'Ax + b'x subject to x >= 0 
## for a symmetric positive definite matrix A and 
## a negative vector b. 
## 
## The solution to the 
## non-negative quadratic programming problem (NQP) is
## based on the MM algorithm (Lee and Seung, 1999)
## 
## Name: MajorizeMinimizeNQP
## 
## Input: 
## A: positive definite K by K matrix with positive entries.
## b: negative vector of length K. 
## init: initial value of x, vector of length K. Optional. 
## maxIter: maximum number of iterations.
## tol: tolerance. Algorithm breaks when the objective function
## decreases tol.
##
## Output:
## x: minimizer of objective function f.
## objfct: final value of objective function f(x).
## xList: values of x in each iteration.
## objfctList: values of f(x) in each iteration.
## time: accumulated time after each iteration. 
##
###############################################################
MajorizeMinimizeNQP <- function(A,b,init,maxIter=1000,tol=0.001){
    ## Objective function
    f = function(x) as.numeric(1/2*crossprod(x, A)%*%x+crossprod(b, x))
    ## Initial values of xList, objFctList and time
    K = ncol(A)
    if (missing(init)) {
        x = matrix(runif(K), nrow=K, byrow=T)
    } else {
        x = init
    }
    xList = x
    objfctList = c(f(x))
    time = c(0)
    start = Sys.time()
    ##-------------------------------------------
    ## Iterative procedure for solving the NQP
    ##-------------------------------------------
    ## Minimization 
    for (i in 1:maxIter){
        x = -b*x/(A %*% x) 
        objfct = f(x)
        # if (abs(objfct - objfctList[i]) <= tol) break
        xList = cbind(xList, x)
        objfctList = c(objfctList, objfct)
        time = c(time, Sys.time() - start)
    } 
    return (list(x=x, objfct=objfct, xList = xList, objfctList = objfctList, time=time))
}