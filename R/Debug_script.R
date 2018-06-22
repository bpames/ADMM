c=2
n=2
p=0.25
q=1-p
m=c*n
tau=0.35
tol=1.0e-4
gamma=3/((q-p)*m*n)
G=matrix(c(1,1,1,0,1,1,1,1,0,1,0,1,0,0,0,1),nrow=4,ncol=4)


# TEST MATRIX SHRINK COMMAND.
S = mat_shrink(K = G, tau = tau)

# Reconstruct matrix.
u <- S$s$u
v <- S$s$v
L <- diag(S$L)

testK <- u %*% L %*% t(v)
K <- u[,1:4]*diag(L)*t(v)

res = ADMM(G = G, c =c, n = n,maxiter = 30, quiet = FALSE)
