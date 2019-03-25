################################################################################
# Simple implementation of EM algorithm for latent-class model
# You can change:
# - The dimension of the observable variables
# - The distributions of the observable variables conditional on the latent
#   variable
# Number of observable variables and dimension of latent variable are hard coded
# Download R for free at r-project.org
################################################################################


# setup ========================================================================

# initialize random no. generator
set.seed(12345)

# no. of vars
k <- 4

# dimension of vars
d <- 4

# dimension of latent variable
p <- 2

# no. of obs
n <- 500


# probability distributions  ===================================================

# pz is the distribution of z
pz <- c(.7, .3)

# pxab is the distribution for x_a when z=b
px11 <- c(.1, .2, .3, .4) 
px12 <- c(.2, .2, .2, .4)
px21 <- c(.5, .25, .2, .05)
px22 <- c(.3, .4, .2, .1)
px31 <- c(.3, .2, .1, .4)
px32 <- c(.5, .2, .2, .1)
px41 <- c(.1, .8, .05, .05)
px42 <- c(.4, .2, .2, .2)

# turn these into a matrix for each value of z
px1 <- cbind(px11, px21, px31, px41)
px2 <- cbind(px12, px22, px32, px42)


# generate data ================================================================

# draw realization of z
# draw z-specific realization of each x
# observed x is the draw corresponding to the draw of z
z <- sample(1:2, n, replace=T, prob=pz)
x1 <- matrix(0, n, k)
x2 <- matrix(0, n, k)
x <- matrix(0, n, k)
for (i in 1:k) {
  x1[,i] <- sample(1:d, n, replace=T, prob=px1[,i])
  x2[,i] <- sample(1:d, n, replace=T, prob=px2[,i])
  x[z==1,i] <- x1[z==1,i]
  x[z==2,i] <- x2[z==2,i]
}


# estimation ===================================================================

# initializations
diff <- 1
tol <- 10^(-5)
iter <- 0

# use perturbed initial parameters to avoid a local max where
# all probabiliy distributions are equal
# P(X_k=x_j | z) calculated as e(b_j)/sum_k e(b_k) to constrain probabilities
# to unit interval
b <- rnorm(2*k*d, 0, 2)

# pi = P(Z=z) = E(q)
pi <- c(.5, .5)

# update q
update.q <-function(b) {
  b1 <- b[1:(d*k)]
  b2 <- b[(d*k+1):(2*d*k)]
  prob1 <- exp(matrix(b1, d, k)) / matrix(rep(apply(exp(matrix(b1, d, k)), 2, sum), d), d, k, byrow=T)
  prob2 <- exp(matrix(b2, d, k)) / matrix(rep(apply(exp(matrix(b2, d, k)), 2, sum), d), d, k, byrow=T)
  l1 <- prob1[x[,1], 1]*prob1[x[,2], 2]*prob1[x[,3], 3]*prob1[x[,4], 4]
  l2 <- prob2[x[,1], 1]*prob2[x[,2], 2]*prob2[x[,3], 3]*prob2[x[,4], 4]
  q <- cbind(pi[1]*l1/(pi[1]*l1+pi[2]*l2), pi[2]*l2/(pi[1]*l1+pi[2]*l2))
  return(q)
}

# expected log likelihood
logl <- function(b){
  b1 <- b[1:(d*k)]
  b2 <- b[(d*k+1):(2*d*k)]
  prob1 <- exp(matrix(b1, d, k)) / matrix(rep(apply(exp(matrix(b1, d, k)), 2, sum), d), d, k, byrow=T)
  prob2 <- exp(matrix(b2, d, k)) / matrix(rep(apply(exp(matrix(b2, d, k)), 2, sum), d), d, k, byrow=T)
  l1 <- prob1[x[,1], 1]*prob1[x[,2], 2]*prob1[x[,3], 3]*prob1[x[,4], 4]
  l2 <- prob2[x[,1], 1]*prob2[x[,2], 2]*prob2[x[,3], 3]*prob2[x[,4], 4]
  l <- q[,1]*log(l1) + q[,2]*log(l2)
  return(-sum(l))
}

while (diff > tol) {
  iter <- iter + 1
  q <- update.q(b)
  b.old <- b
  opt <- optim(b, logl)
  b <- opt$par
  pi <- apply(q, 2 , mean)
  diff <- max(abs(b-b.old))
  print(iter)
  print(pi)
  print(diff)
}

b1 <- b[1:(d*k)]
b2 <- b[(d*k+1):(2*d*k)]
prob1 <- exp(matrix(b1, d, k)) / matrix(rep(apply(exp(matrix(b1, d, k)), 2, sum), d), d, k, byrow=T)
prob2 <- exp(matrix(b2, d, k)) / matrix(rep(apply(exp(matrix(b2, d, k)), 2, sum), d), d, k, byrow=T)

