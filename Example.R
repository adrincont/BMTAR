# Example mtaregim:
orders = list(p = 2,q = 1,d = 1)
Phi = list(phi2 = matrix(c(0.1,0.6,-0.4,0.5),2,2, byrow = T))
Beta = list(beta1 = matrix(c(0.3,-0.4),2, 1))
Delta = list(delta1 = matrix(c(0.6,1),2,1))
Sigma = matrix(c(1,0.6,0.6,1.5),2,2,byrow = T)
cs = matrix(c(1,-1),nrow = 2)
Ri = mtaregim(orders = orders,Phi = Phi,Beta = Beta,Delta = Delta,Sigma = Sigma,cs = cs)
