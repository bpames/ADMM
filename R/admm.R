#' Sample matrix
#'
#' Generates binary (M,N) - matrix sampled from dense (m,n) - submatrix
#' Let U* and V* be m and n index sets.
#' For each i in U*, j in V* we let a_ij = 1 with probability q and 0  otherwise.
#' For each remaining (ij) we set a_ij = 1 with probability p < q and take a_ij = 0 otherwise.
#'
#' @param M number of rows in sampled matrix
#' @param N number of rows in sampled matrix
#' @param c natural number used to calculate number of rows in dense submatrix
#' @param n natural number used to calculate number of rows in dense submatrix
#' @param p density outside planted submatrix
#' @param q density inside planted submatrix
#' @return Matrix G sampled from the planted dense (mn)-submatrix model, dense sumbatrix X0, matrix Y0 used to count the number of disagreements between G and X0
#' @export

matrix_plantsubm <- function(M,N,c,n,p,q){
  m <-c*n
  gamma<-3/((q-p)*m*n)

  #Initialize G
  G <- matrix(runif(M*N),M,N)
  G <- ceiling(G-(1-p))

  #Make dense submatrix
  start_position <- matrix(runif(m*n),m,n)
  G[1:m, 1:n]  <- ceiling(start_position-(1-q))

  #Get X, fill the dense submatrix
  X0 <- matrix(0L, nrow=M, ncol=N)
  X0[1:m, 1:n] <- matrix(1L, nrow=m, ncol=n)

  #Get Y
  Y0 <- matrix(0L, nrow=M, ncol=N)
  Y0[1:m,1:n] <- matrix(1L, nrow=m,ncol=n) - G[1:m,1:n]

  list(sampled_matrix = G, dense_submatrix = X0, disagreements = Y0)
}




#' Soft threshholding operator.
#'
#' Applies the shrinkage operator for singular value tresholding.
#' @param K matrix
#' @param tau regularization parameter
#' @return Matrix
#' @export

mat_shrink <- function (K, tau){

  r <- dim(K)[1]
  c <- dim(K)[2]

  s <- svd(K, nu=r, nv=c)
  L <- pmax(s$d-tau,0)

  if (r < c) {
    K <- s$u*diag(L)*t(s$v[,1:r])
  } else {
    K <- s$u[,1:c]*diag(L)*t(s$v)
  }
  return(K)
}






#' ADMM
#'
#' Iteratively solves the convex optimization problem using ADMM.
#'
#'       min     |X|_* + gamma* |Y|_1 + 1_Omega_W (W) + 1_Omega_Q (Q) + 1_Omega_Z (Z)
#'       st      X - Y = 0
#'               X = W
#'               X = Z
#'
#'
#' where Omega_W (W), Omega_Q (Q), Omega_Z (Z) are the sets:
#'       Omega_W = {W in R^MxN | e^TWe = mn}
#'       Omega_Q = {Q in R^MxN | Projection of Q on not N = 0}
#'       Omega_Z = {Z in R^MxN | Z_ij <= 1 for all (i,j) in M x N}
#' 1_S is the indicator function of the set S in R^MxN such that 1_S(X) = 0 if X in S and +infinity otherwise
#'
#' @param G: sampled binary matrix
#' @param c: natural number used to calculate number of rows in dense submatrix
#' @param n: number of columns in dense submatrix
#' @param gamma:  l1 regularization parameter
#' @param tau:    penalty parameter for equality constraint violation
#' @param tol:    stopping tolerance in algorithm
#' @param maxiter: maximum number of iterations of the algorithm to run
#' @param quiet: toggles between displaying intermediate statistics
#' @return Rank one matrix with mn nonzero entries, matrix Y that is used to count the number of disagreements between G and X
#' @export

ADMM <- function(G, c, n, tau = 0.35, gamma = 0.75, tol = 1.0e-4,maxiter, quiet = TRUE){

  m <- c*n
  mu <- 1/tau
  M <- dim(G)[1]
  N <- dim(G)[2]

  for(i in 1:M)
  {
    for(j in 1:M)
    {
      if(is.na(G[i,j])|(!(G[i,j] == 1 | G[i,j] == 0)))
      {
        stop('The argument "G" must be binary matrix')
      }
    }
  }

  if(m > M | n > N){
    warning ('Dimensions of dense submatrix exceed dimensions of  initial matrix G')
  }

  #Initialize
  W <- matrix(1L, nrow = M, ncol = N)*m*n/(M*N)
  X <- W
  Y <- X
  Z <- X
  Q <- X - Y
  LambdaQ <- matrix(0L, nrow=nrow(X), ncol=ncol(X))
  LambdaZ <- matrix(0L, nrow=nrow(X), ncol=ncol(X))
  LambdaW <- matrix(0L, nrow=nrow(X), ncol=ncol(X))


convergence <- 0
iter <- 0

    while (iter < maxiter && convergence==0){

      #Update Q
      Q_old <- Q
      Q <- X - Y + mu*LambdaQ
      Q <- Q*G

      #Update X
      mat<- 1/3*(Y + Q + Z + W - mu*(LambdaQ + LambdaW + LambdaZ))
      X = mat_shrink(K=mat, tau=1/(3*tau))

      #Update Y
      A <-X - Q - gamma*matrix(1L, nrow = M, ncol = N)*mu + LambdaQ*mu
      B <- matrix(0L, nrow = M, ncol = N)
      Y <- pmax(A, B)

      #Update W
      W_old <- W
      newW <- X + mu*LambdaW
      alfa <- (tail(m, n=1)*tail(n, n=1) - sum(as.vector(newW)))/(M*N)
      W <- newW + alfa*matrix(1L, nrow = M, ncol = N)

      #Update Z
      Z_old <- Z
      Z <- X + mu*LambdaZ
      D <- matrix(0L, M, N)
      Z1 <- pmax(Z, D)
      E <- matrix(1L, M, N)
      Z <- pmin(Z1, E)

      #Update dual variables
      LambdaQ <- LambdaQ + tau*(X - Y - Q)
      LambdaW <- LambdaW + tau*(X - W)
      LambdaZ <- LambdaZ + tau*(X - Z)

      #Check convergence
      #primal
      NZ <- norm(X - Z, type="F")
      NW <- norm(X - W, type="F")
      NQ <- norm(X - Y - Q, type="F")

      errP <- max(NZ, NW, NQ)/norm(X, type="F")

      #Dual feasibility
      NDz <- norm(Z - Z_old, type="F")
      NDw <- norm(W - W_old, type="F")
      NDp <- norm(Q - Q_old, type="F")

      errD <- max(NDz, NDw, NDp)/norm(X, type="F")

      if(errP < tol && errD < tol) {
        convergence=1
        break
      }
      else {
        convergence=0
      }
iter <- iter +1

      print(sprintf('iter: %d, primal_gap: %e, dual_gap: %e',iter, errP, errD))


    }
  results <- list(X=X, Y=Y)
  return(results)
}


