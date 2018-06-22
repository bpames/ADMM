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
