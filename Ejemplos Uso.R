library(forecast) #objeto ts y autoplot
library(gridExtra) #graficas ggplot juntas
library(mvtnorm) #Normal multivariada
library(LaplacesDemon) #Wishart inversa
library(ggplot2) #graficas
library(expm) #raiz y potencia de matrices
library(tsDyn) # simulacion VAR
library(ks) #funcion inversa de vec
library(Brobdingnag) #manejo de grandes numeros 
rm(list = ls()[!{ls() %in% c('mtaregim','mtarsim','mtarns','mtarstr')}])
# Simulacion proceso Ut ####
Tlen = 1000
nu = 1
sigmaru = matrix(c(1,0.4,0.4,2),nu + 1,nu + 1)
Au = matrix(c(0.5,0.1,0.4,0.5),nu + 1,nu + 1,byrow = T)
ut = tsDyn::VAR.sim(B = Au,n = Tlen,lag = 1,include = c("none"),varcov = sigmaru)
forecast::autoplot(ts(ut),facets = TRUE) + theme_bw()
# Primer ejemplo Model 2 pg.26 ####
R1 = mtaregim(orders = list(p = 2,q = 1,d = 1),
            Phi = list(phi1 = matrix(c(0.5,-0.2,-0.2,0.8),2,2,byrow = T),
                     phi2 = matrix(c(0.1,0.6,-0.4,0.5),2,2,byrow = T)),
            Beta = list(beta1 = matrix(c(0.3,-0.4),2,1)),
            Delta = list(delta1 = matrix(c(0.6,1),2,1)),
            Sigma = matrix(c(1,0.6,0.6,1.5),2,2,byrow = T),
            cs = matrix(c(1,-1),nrow = 2)
            )
R2 = mtaregim(orders = list(p = 1,q = 0,d = 0),
            Phi = list(phi1 = matrix(c(0.3,0.5,0.2,0.7),2,2,byrow = T)),
            Sigma = matrix(c(2.5,0.5,0.5,1),2,2,byrow = T),
            cs = matrix(c(5,2),nrow = 2)
            )
yt = mtarsim(N = Tlen,Rg = list(R1 = R1,R2 = R2),r = -0.22,Ut = Ut)
forecast::autoplot(ts(yt$Yt),facets = TRUE)
rm(list = ls()[!{ls() %in% c('yt','R1','R2','ut','mtaregim','mtarsim','mtarns','mtarstr')}])
# Estimacion con Sigma, r y estructurales conocidos ####
estim1 = mtarns(Yt = yt$Yt,Ut = ut,l = 2,r = qnorm(0.4),pj = yt$pj,dj = yt$dj,qj = yt$qj, 
             Sigma = list(R1 = R1$sigma,R2 = R2$sigma),niter = 5000,chain = TRUE)
# Estimacion con r y estructurales conocidos ####
estim2 = mtarns(Yt = yt$Yt,Ut = ut,l = 2,r = qnorm(0.4),
                pj = yt$pj,dj = yt$dj,qj = yt$dj,niter = 5000,chain = TRUE)
# Estimacion con sigma y estructurales conocidos ####
estim3 = mtarns(Yt = yt$Yt,Ut = ut,l = 2,pj = yt$pj,dj = yt$dj,qj = yt$qj, 
                Sigma = list(R1 = R1$sigma,R2 = R2$sigma),niter = 5000,chain = TRUE)
# Estimacion con estructurales conocidos ####
estim4 = mtarns(Yt = yt$Yt,Ut = ut,l = 2,pj = c(yt$pj),dj = yt$dj,qj = yt$dj,niter = 5000,chain = TRUE)

# Estimacion con estructurales conocidos dando iniciales ####
#estim = mtarns(Yt = yt$Yt,Ut = ut,l = 2,r = qnorm(0.4),pj = yt$pj,dj = yt$dj,qj = yt$qj, 
#              sigmaini = list(list(diag(2),2),list(diag(2),2)),
#              thetaini = list(list(matrix(0,14,1),diag(14)),list(matrix(0,10,1),diag(10))) ,niter = 1000)

# Estimacion con l conocido ####
estruc = mtarstr(Yt = yt$Yt,Ut = ut,l = 2,orders = list(pjmax = c(3,3),qjmax = c(2,2),djmax = c(2,2)),niter = 2000,method = 'SSVS',chain = T)

estimKUO = mtarstr(Yt = yt$Yt,Ut = ut,l = 2,pjmax = c(2,2),qjmax = c(2,2),djmax = c(2,2),niter = 2000,method = 'KUO',chain = TRUE) 
estimSSVS = mtarstr(Yt = yt$Yt,Ut = ut,l = 2,pjmax = c(2,2),qjmax = c(2,2),djmax = c(2,2),niter = 2000,method = 'SSVS',chain = TRUE)

save.image("~/Val/Unal/MTAR/nsstr.RData")
# Segundo ejemplo Model 3 pg.27 ####
R1=mtaregim(orders=list(p=1,q=0,d=0),
            Phi=list(phi1=matrix(c(-0.9,0,0.2,-0.5),2,2,byrow = T)),
            Sigma=matrix(c(1,0.3,0.3,4),2,2,byrow = T),
            cs=matrix(c(2,1),2,1))
R2=mtaregim(orders=list(p=2,q=1,d=0), 
            Phi=list(phi1=matrix(c(0.7,0,0,0.6),2,2,byrow = T),
                     phi2=matrix(c(0.8,0.2,0,-0.4),2,2,byrow = T)),
            Beta=list(beta1=matrix(c(1.2,-0.8),2,1)),
            Sigma=matrix(c(1,0,0,1),2,2,byrow = T),
            cs=matrix(c(0.4,-4),2,1))
R3=mtaregim(orders=list(p=3,q=2,d=1),
            Phi=list(phi3=matrix(c(-0.8,0,0.2,0.8),2,2,byrow = T)),
            Beta=list(beta2=matrix(c(-0.6,0.7),2,1)),
            Delta=list(delta1=matrix(c(0.6,2),2,1)),
            Sigma=matrix(c(2,-0.4,-0.4,1),2,2,byrow = T),
            cs=matrix(c(-3,2),2,1))
yt=mtarsim(N=Tlen,Rg=list(R1=R1,R2=R2,R3=R3),r=c(qnorm(0.25),qnorm(0.75)),Ut=Ut)
forecast::autoplot(ts(yt$Yt),facets = T)
