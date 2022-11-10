#### #### #### #### #### #### #### ####
#### Stima modelli Lineare - MCMC
#### #### #### #### #### #### #### ####

# definiamo la funzione che ci permette di simulare da una
# normale multivariata e calcolarne la densità
rmnorm=function(n = 1, mean = rep(0, d), varcov)
{
  d <- if (is.matrix(varcov))
    ncol(varcov)
  else 1
  z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
  y <- t(mean + t(z))
  return(y)
}
dmnorm=function (x, mean = rep(0, d), varcov, log = FALSE)
{
  d <- if (is.matrix(varcov))
    ncol(varcov)
  else 1
  if (d > 1 & is.vector(x))
    x <- matrix(x, 1, d)
  n <- if (d == 1)
    length(x)
  else nrow(x)
  X <- t(matrix(x, nrow = n, ncol = d)) - mean
  Q <- apply((solve(varcov) %*% X) * X, 2, sum)
  logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
  logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
  if (log)
    logPDF
  else exp(logPDF)
}

## ## ## ## ## ##
## Modello Lneare - MCMC
## ## ## ## ## ##

ModLin = function(formula,  beta.mean, beta.variance,sigma2.a, sigma2.b, start.beta, start.sigma2, iter, burnin, thin, Data, sd.prop = 1)
{

  ### Definiamo l'osservazioni e la matrice X
  Cov     = model.matrix(formula,Data)
  NameObs = as.character(terms(formula)[[2]])
  Obs     = Data[,NameObs,drop=F] # drop serve a fare in modo che il risultato sia ancora una matrice

  nbeta = ncol(Cov)
  nobs  = nrow(Cov)

  ### definisco quanti campioni a posteriori devo salvare
  ### e gi oggetti che ritorna la funzione
  NsampleSave     = floor((iter-burnin)/thin)
  BetaSave        = matrix(NA,ncol=nbeta, nrow=NsampleSave )
  sigma2Save      = matrix(NA,ncol=1, nrow=NsampleSave )

  ### definisco i valori correnti di beta e sigma2
  betaMCMC   = start.beta
  sigma2MCMC = start.sigma2



  ## oggetti utili per l'MCMC
  appSamp = burnin
  i = 1
  XtX = t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
  Xty = t(Cov[i,,drop=F])*Obs[i,1]
  for(i in 2:nobs)
  {
    XtX = XtX+t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
    Xty = Xty+t(Cov[i,,drop=F])*Obs[i,1]
  }
  for(iMCMC in 1:NsampleSave)
  {
    for(jMCMC in 1:appSamp)
    {
      ## ## ## ## ## ## ## ##
      ## campioniamo Beta
      ## ## ## ## ## ## ## ##

      # prior parameters
      Vpost = solve(beta.variance)
      Mpost = solve(beta.variance)%*%matrix(beta.mean,ncol=1)

      # contributi delle osservazioni
      Vpost = Vpost+XtX/sigma2MCMC
      Mpost = Mpost+Xty/sigma2MCMC

      Vpost = solve(Vpost)
      Mpost = Vpost%*%Mpost

      #simlazione Gibbs
      betaMCMC = rmnorm(1,Mpost,Vpost)

      ## ## ## ## ## ## ## ##
      ## campioniamo sigma2 _ Metropolis
      ## ## ## ## ## ## ## ##

      # propongo  un nuovo valore
      sigma2prop = exp(rnorm(1,log(sigma2MCMC),sd.prop))

      # calcolo l'acceptance rate
      # potete scrivelo in due modi
      # 1) la distribuzione proposal è un log normale (esponenziale di una normale) e la prior gamma
      # 2) invece di sigma2 lavoriamo con tau=log(sigma2), che ha proposal simmetrica, e la prior la dobbiamo trovare
      # con la regola di trasformazioni di variabili:
      # f(tau) = f(sigma2(tau)) |d sigma2(tau) / d tau|
      # il second metodo, che è quello che usiamo, ha il vantaggio
      # che essendo la proposa simmetrica, si semplifica nel
      # calcolo di alpha

      # prior contribtion di tau
      logalpha = 0
      logalpha = logalpha+dgamma(sigma2prop, shape=sigma2.a, rate= sigma2.b,log=T)+log(sigma2prop)
      logalpha = logalpha-(dgamma(sigma2MCMC, shape=sigma2.a, rate= sigma2.b,log=T)+log(sigma2MCMC))

      # likelihood contribution
      for(i in 1:nobs)
      {
        logalpha = logalpha+dnorm(Obs[i,],Cov[i,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2prop^0.5, log=T)
        logalpha = logalpha-dnorm(Obs[i,],Cov[i,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2MCMC^0.5, log=T)
      }

      alpha = min(1,exp(logalpha))

      u = runif(1,0,1)
      if(u<alpha)
      {
        sigma2MCMC = sigma2prop
      }
    }
    appSamp = thin

    BetaSave[iMCMC,]    = betaMCMC
    sigma2Save[iMCMC,]  = sigma2MCMC
  }
  return(list(Beta =BetaSave, sigma2 = sigma2Save ))
}


#### #### #### #### #### #### ####
#### TESTIAMO IL CODICE CON DATI SIMULATI
#### #### #### #### #### #### ####
set.seed(10)
n       = 100
beta    = c(1,2,1,-1)
x1      = runif(n, 0,1)
x2      = runif(n, -1,1)
x3      = rnorm(n, 0,1)
sigma2  = 1
y       = beta[1]+beta[2]*x1+beta[3]*x2+beta[4]*x3+rnorm(n, 0,sigma2^0.5)

Data = data.frame(y =y, x1=x1,x2=x2,x3=x3)


Results = ModLin(
  y ~ x1+x2+x3,
  beta.mean = rep(0,4),
  beta.variance=diag(100,4),
  sigma2.a=1,
  sigma2.b=1,
  start.beta=rep(0,4),
  start.sigma2=1,
  iter=3000,
  burnin=1,
  thin=1,
  Data=Data
)

par(mfrow=c(2,2))
plot(Results$Beta[,1], type="l")
acf(Results$Beta[,1])

plot(Results$Beta[,2], type="l")
acf(Results$Beta[,2])

plot(Results$Beta[,3], type="l")
acf(Results$Beta[,3])

plot(Results$Beta[,4], type="l")
acf(Results$Beta[,4])

plot(Results$sigma2[,1], type="l")
acf(Results$sigma2[,1])
par(mfrow=c(1,1))

#### #### #### #### #### #### ####
#### Utilizziamo dei dati reali
#### #### #### #### #### #### ####

library(datasets)
library(coda)


Dataset = trees
?trees

# la formulazione giusta è probabilmente
# Volume ~= c * Height * Girth^2
# noi quindi la modelliamo come
# E(log(Volume)) = b_1+b_2log(Height)+b_3log(Girth)

## facciamo dei plot
plot(log(Dataset))

# creiamo il dataset

DataLog = log(Dataset)

# il modello lineare
Results = ModLin(
  Volume ~ Girth+Height,
  beta.mean = rep(0,3),
  beta.variance=diag(100,3),
  sigma2.a=1,
  sigma2.b=1,
  start.beta=rep(0,3),
  start.sigma2=1,
  iter=3000,
  burnin=1,
  thin=1,
  Data=DataLog
)

par(mfrow=c(4,2))
plot(Results$Beta[,1], type="l")
acf(Results$Beta[,1])

plot(Results$Beta[,2], type="l")
acf(Results$Beta[,2])

plot(Results$Beta[,3], type="l")
acf(Results$Beta[,3])

plot(Results$sigma2[,1], type="l")
acf(Results$sigma2[,1])
par(mfrow=c(1,1))

summary(DataLog)
## trasformiamo l'output in mcmc (coda)
MCMC_out = cbind(Results$Beta,Results$sigma2)
colnames(MCMC_out) = c(paste("Beta",1:3,sep=""), "sigma")
MCMC_out = as.mcmc(MCMC_out)
summary(MCMC_out)

## calcoliamo i residui
residui_mat = matrix(NA, ncol=nrow(DataLog), nrow=nrow(MCMC_out))
for(i in 1:nrow(MCMC_out))
{
  residui_mat[i,] = DataLog$Volume - as.matrix(cbind(1,(DataLog[,c("Girth","Height")])))%*%matrix(Results$Beta[i,], ncol=1)
}
par(mfrow=c(4,3))
for(i in 1:12)
{
  plot(density(residui_mat[,i]), main=paste("res", i))
}

# calcoliamo la media dei residui e facciamo dei plot
res_mean = colMeans(residui_mat)

par(mfrow=c(1,3))
plot(density(res_mean)) # sembrano bimodali
plot(DataLog[,"Girth"],res_mean)
plot(DataLog[,"Height"],res_mean)

## calcoliamo la matrice della "previsione"
prev_mat = matrix(NA, ncol=nrow(DataLog), nrow=nrow(MCMC_out))
for(i in 1:nrow(MCMC_out))
{
  prev_mat[i,] = rnorm(nrow(DataLog),as.matrix(cbind(1,(DataLog[,c("Girth","Height")])))%*%matrix(Results$Beta[i,], ncol=1),Results$sigma2[i,]^0.5 )
}

# vediamo se i dati cadono nelle distribuzione predittive
par(mfrow=c(4,3))
for(i in 1:12)
{
  plot(density(prev_mat[,i]), main=paste("res", i))
  abline(v=DataLog$Volume[i], col=2)
}

# calcoliamo la forma della distribuzione predittiva dell'osservazione is
is = 1
yseq = seq(2,3,by=0.01)
pred_dens =matrix(0,ncol=length(yseq), nrow=nrow(MCMC_out))
for(i in 1:nrow(MCMC_out))
{
  pred_dens[i,] = dnorm(yseq,sum(cbind(1,DataLog[is,c("Girth","Height")])*Results$Beta[i,]), Results$sigma2[i,]^0.5 )
}
# plottiamo i risultati e confrontiamoli con i capioni dalla predittiva
plot(yseq,colMeans(pred_dens))
lines(density(prev_mat[,is]), col=2)
