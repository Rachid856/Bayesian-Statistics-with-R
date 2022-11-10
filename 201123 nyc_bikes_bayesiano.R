#### #### #### #### #### #### ####
#### Poisson bayesiano e Model Selection
#### #### #### #### #### #### ####

rm(list=ls())

#### librerie
library(R2jags)
library(sjmisc)

#### #### #### #### #### ####
#### The data
#### #### #### #### #### ####
#### #### DESCRIZIONE
#### Numero (in migliagia) di bici osservate passare sul ponte Williamsburg
#### di New York, dal prima aprile al 31 ottobre. Ordinate temporalmente
####
#### weekend: giorno della settimana
#### hightemp: temperatura massima giornaliera
#### lowtemp: temperatura massima giornaliera
#### precip_rain: quantità di pioggia giornaliera
#### precip_snow: quantità di neve giornaliera
#### time: variabile che indica la differenza tra il giorno osservato e
#### il 31 marzo
#### count: numero di biciclette (in migliaia)
#### #### #### #### #### ####

# # # # # # # # # # # # # # #
# l'.Rdata nyc_bike_bayesiano contiene i risultati dei modelli.
# # # # # # # # # # # # # # #

setwd("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/didattica/Modelli Statistici/Codici R/PoissonBayesiano/")
nyc_bikes = read.csv( "nyc_bike_counts.csv")
summary(nyc_bikes)

# facciamo dei plot
plot(nyc_bikes)
plot(nyc_bikes$hightemp,nyc_bikes$lowtemp)

### ### ### ### ### ###
### Modello frequentista
### ### ### ### ### ###

# abbiam problemi di collinearità nelle temperature
mod_freq1 = glm(count ~ weekday+hightemp+lowtemp+precip_rain+precip_snow+time, data = nyc_bikes, family = "poisson")
summary(mod_freq1)

mod_freq2 = glm(count ~ weekday+lowtemp+precip_rain+precip_snow+time, data = nyc_bikes, family = "poisson")
summary(mod_freq2)

mod_freq3 = glm(count ~ weekday+hightemp+precip_rain+precip_snow+time, data = nyc_bikes, family = "poisson")
summary(mod_freq3)

# vediamo quale modello è migliore, secondo l'AIC
AIC(mod_freq1,mod_freq2,mod_freq3)
# e proviamo a fare model selection con la funzione step
mod_step = step(mod_freq1)
summary(mod_step)


#### #### #### #### #### #### ####
#### Modello Bayes - stime e model selection
#### #### #### #### #### #### ####

# stimiamo il modello con tutte le covariate
# Matrice X
X = model.matrix(~ weekday+hightemp+lowtemp+precip_rain+precip_snow+time ,data=  nyc_bikes)
# Vettore y
Y = nyc_bikes$count[,drop=F]
# indici vari
n = nrow(X)
p = ncol(X)
nchain = 2
# lista con tutti gli elementi che servono al modello
dataList = list(
	Y   = Y,
	X   = X,
	n   = n
)

#### Modello Bayesiano
mod1_string <- 'model {
  for(i in 1:n) {
    Y[i]  ~  dpois(lambda[i])
		lambda[i] = exp(mu[i])
    mu[i] <- 	b[1]+ #(Intercept)
							b[2]*X[i,2]+ #weekdayMonday
							b[3]*X[i,3]+ #weekdaySaturday
							b[4]*X[i,4]+ #weekdaySunday
							b[5]*X[i,5]+ #weekdayThursday
							b[6]*X[i,6]+ #weekdayTuesday
							b[7]*X[i,7]+ #weekdayWednesday
							b[8]*X[i,8]+ #hightemp
							b[9]*X[i,9]+ #lowtemp
							b[10]*X[i,10]+ #precip_rain
							b[11]*X[i,11]+ #precip_snow
							b[12]*X[i,12]  # time
  }
	b[1] ~ dnorm(0,1/1000)
	b[2] ~ dnorm(0,1/1000)
	b[3] ~ dnorm(0,1/1000)
	b[4] ~ dnorm(0,1/1000)
	b[5] ~ dnorm(0,1/1000)
	b[6] ~ dnorm(0,1/1000)
	b[7] ~ dnorm(0,1/1000)
	b[8] ~ dnorm(0,1/1000)
	b[9] ~ dnorm(0,1/1000)
	b[10] ~ dnorm(0,1/1000)
	b[11] ~ dnorm(0,1/1000)
	b[12] ~ dnorm(0,1/1000)
}'

# inizializzazioni
set.seed(1)
inits =  list(
  list("b" = runif(p,-1,1)),
  list("b" = runif(p,-1,1))
)

# parametri da salvare
SavePar = c("b", "lambda")

# fittiamo il modello
set.seed(1)
fit_glm1 = jags(
			data 								= dataList,
			inits 							= inits,
			parameters.to.save 	= SavePar,
			model.file 					= textConnection(mod1_string), # il model.file dovrebbe essere un file,
			 																									# "textConnection" crea il file e lo carica
			n.chains 						= nchain,
			n.iter 							= 10000,
			n.burnin 						= 5000,
			n.thin 							= 2,
			DIC 								= T
)
#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####
beta_post_chain1 = fit_glm1$BUGSoutput$sims.array[-c(1:100),1,1:p]
beta_post_chain2 = fit_glm1$BUGSoutput$sims.array[-c(1:100),2,1:p]
colnames(beta_post_chain1) = colnames(beta_post_chain2) = colnames(X)

# confrontiamo le catene per vedere se il modello è a convergenza
par(mfrow=c(3,2))
for(i in 1:p)
{
	plot(beta_post_chain1[,i], type="l", main= colnames(beta_post_chain1)[i])
	lines(beta_post_chain2[,i],col=2)
	# if(p%%6)
	# {
	# 	plot.new()
	# }
}

# stimiamo le medie e gli intervalli di credibilità
beta_mean = round(colMeans(beta_post_chain1),5)
beta_CI1  = round(apply(beta_post_chain1,2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1,2,quantile,prob=0.975),5)

Beta_Res = cbind(beta_mean,beta_CI1,beta_CI2)
Beta_Res
# abbiamo lo stesso problema interpretativo dovuto alla collinearità

### vediamo le distribuzioni di alcuni parametri
plot(density(beta_post_chain1[,"time"]))
plot(density(beta_post_chain1[,"precip_snow"]))

# c'è molta autocorrelazione tra le temperature
plot(beta_post_chain1[,"hightemp"],beta_post_chain1[,"lowtemp"])
acf(beta_post_chain1[,-c(1:7)])

#### #### #### #### #### ####
#### Prior Spike-and-Slab -  Modello gerarchico
#### #### #### #### #### ####
# se theta è un parametro e
# come prior usiamo
# theta| z ~ N(0, zs1+(1-z)s2)
# con s1 piccola e s2 grande e
# z ~ Bin(1,pi).
# allora la prior su theta è  "Spike-and-Slab"

mod2_string <- 'model {
  for(i in 1:n) {
    Y[i]  ~  dpois(lambda[i])
		lambda[i] = exp(mu[i])
    mu[i] <- 	b[1]+ #(Intercept)
							b[2]*X[i,2]+ #weekdayMonday
							b[3]*X[i,3]+ #weekdaySaturday
							b[4]*X[i,4]+ #weekdaySunday
							b[5]*X[i,5]+ #weekdayThursday
							b[6]*X[i,6]+ #weekdayTuesday
							b[7]*X[i,7]+ #weekdayWednesday
							b[8]*X[i,8]+ #hightemp
							b[9]*X[i,9]+ #lowtemp
							b[10]*X[i,10]+ #precip_rain
							b[11]*X[i,11]+ #precip_snow
							b[12]*X[i,12]  # time
  }
	var1 = 10000
	var2 = 0.0000000001

	b_app[1] ~ dnorm(0,1/var1)
	b_app[2] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[3] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[4] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[5] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[6] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[7] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[8] ~ dnorm(0,1/var1)
	b_app[9] ~ dnorm(0,1/var2)
	b_app[10] ~ dnorm(0,z[10]*1/var1+(1-z[10])*1/var2)
	b_app[11] ~ dnorm(0,z[11]*1/var1+(1-z[11])*1/var2)
	b_app[12] ~ dnorm(0,z[12]*1/var1+(1-z[11])*1/var2)

	b[1]  	= b_app[1]
	b[2]   	= b_app[2]
	b[3]   	= b_app[3]
	b[4]   	= b_app[4]
	b[5]   	= b_app[5]
	b[6]   	= b_app[6]
	b[7]   	= b_app[7]
	b[8]   	= (z[8]*			b_app[8]+(1-z[8])*b_app[9])
	b[9]   	= (1-z[8])*	b_app[8]+z[8]*b_app[9]
	b[10]   = b_app[10]
	b[11]   = b_app[11]
	b[12]   = b_app[12]


	# Alcune z non servono e non vengono utlizzate,
	# ma le tengo per comodità
	z[1] ~ dbinom(prob[1],1)
	z[2] ~ dbinom(prob[1],1)
	z[3] ~ dbinom(prob[1],1)
	z[4] ~ dbinom(prob[1],1)
	z[5] ~ dbinom(prob[1],1)
	z[6] ~ dbinom(prob[1],1)
	z[7] ~ dbinom(prob[1],1)
	z[8] ~ dbinom(prob[2],1)
	z[9] ~ dbinom(prob[1],1)
	z[10] ~ dbinom(prob[1],1)
	z[11] ~ dbinom(prob[1],1)
	z[12] ~ dbinom(prob[1],1)

	prob[1] ~ dunif(0,1)
	prob[2] ~ dunif(0,1)
	prob[3] ~ dunif(0,1)



}'

# inizializzazioni
set.seed(123)
inits =  list(
  list("b_app" = runif(0,-1,1), "z" = rep(1,p), "prob"=rep(1,3)),
  list("b_app" = runif(0,-1,1), "z" = rep(1,p), "prob"=rep(1,3))
)

# parametri da salvare
SavePar = c("b", "lambda", "z","prob")

set.seed(123)
# fittiamo il modello
fit_glm2 = jags(
			data 								= dataList,
			inits 							= inits,
			parameters.to.save 	= SavePar,
			model.file 					= textConnection(mod2_string), # il model.file dovrebbe essere un file,
			 																									# "textConnection" crea il file e lo carica
			n.chains 						= nchain,
			n.iter 							= 10000,
			n.burnin 						= 5000,
			n.thin 							= 2,
			DIC 								= T # simile all'AIC, serve per il model selection
)

#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####
beta_post_chain1 = fit_glm2$BUGSoutput$sims.array[,1,1:p]
beta_post_chain2 = fit_glm2$BUGSoutput$sims.array[,2,1:p]
colnames(beta_post_chain1) = colnames(beta_post_chain2) = colnames(X)

# confrontiamo le catene per vedere se il modello è a convergenza
par(mfrow=c(3,2))
for(i in 1:p)
{
	plot(beta_post_chain1[,i], type="l", main= colnames(beta_post_chain1)[i])
	lines(beta_post_chain2[,i],col=2)
}

# stimiamo le medie e gli intervalli di credibilità
beta_mean = round(colMeans(beta_post_chain1),5)
beta_CI1  = round(apply(beta_post_chain1,2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1,2,quantile,prob=0.975),5)

Beta_Res = cbind(beta_mean,beta_CI1,beta_CI2)
Beta_Res
#### #### #### ####
#### vediamo la variabile z
z_post = fit_glm2$BUGSoutput$sims.array[,1,c("z[1]","z[8]","z[10]","z[11]","z[12]")]
colnames(z_post) = c("Season","Temperature","Rain","Snow","Time")
colMeans(z_post)
# 1 in temperatura ci dice che sceglie hightemp


### ### ### ### ### ### ### ### ###
### possiamo usare qualcosa di simile per
### decidere che trasformazione utilizzare per time

# Matrice X  con time in scala log e radice quadrata
X = model.matrix(~ weekday+hightemp+lowtemp+precip_rain+precip_snow+I(log(time))+I(time^0.5) ,data=  nyc_bikes)
# Vettore y
Y = nyc_bikes$count[,drop=F]
# indici vari
n = nrow(X)
p = ncol(X)
nchain = 2

# lista con tutti gli elementi che servono al modello
dataList = list(
	Y   = Y,
	X   = X,
	n   = n
)

mod3_string <- 'model {
  for(i in 1:n) {
    Y[i]  ~  dpois(lambda[i])
		lambda[i] = exp(mu[i])
    mu[i] <- 	b[1]+ #(Intercept)
							b[2]*X[i,2]+ #weekdayMonday
							b[3]*X[i,3]+ #weekdaySaturday
							b[4]*X[i,4]+ #weekdaySunday
							b[5]*X[i,5]+ #weekdayThursday
							b[6]*X[i,6]+ #weekdayTuesday
							b[7]*X[i,7]+ #weekdayWednesday
							b[8]*X[i,8]+ #hightemp
							b[9]*X[i,9]+ #lowtemp
							b[10]*X[i,10]+ #precip_rain
							b[11]*X[i,11]+ #precip_snow
							b[12]*X[i,12]+  # log(time)
							b[13]*X[i,13]  # time^0.5
  }
	var1 = 10000
	var2 = 0.0000000001

	b_app[1] ~ dnorm(0,1/var1)
	b_app[2] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[3] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[4] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[5] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[6] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[7] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[8] ~ dnorm(0,1/var1)
	b_app[9] ~ dnorm(0,1/var2)
	b_app[10] ~ dnorm(0,z[10]*1/var1+(1-z[10])*1/var2)
	b_app[11] ~ dnorm(0,z[11]*1/var1+(1-z[11])*1/var2)
	b_app[12] ~ dnorm(0,1/var1)
	b_app[13] ~ dnorm(0,1/var2)

	b[1]  	= b_app[1]
	b[2]   	= b_app[2]
	b[3]   	= b_app[3]
	b[4]   	= b_app[4]
	b[5]   	= b_app[5]
	b[6]   	= b_app[6]
	b[7]   	= b_app[7]
	b[8]   	= (z[8]*			b_app[8]+(1-z[8])*b_app[9])
	b[9]   	= (1-z[8])*	b_app[8]+z[8]*b_app[9]
	b[10]   = b_app[10]
	b[11]   = b_app[11]
	b[12]   = (z[12]*			b_app[12]+(1-z[12])*b_app[13])
	b[13]   = (1-z[12])*	b_app[12]+z[12]*b_app[13]


	z[1] ~ dbinom(prob[1],1)
	z[2] ~ dbinom(prob[1],1)
	z[3] ~ dbinom(prob[1],1)
	z[4] ~ dbinom(prob[1],1)
	z[5] ~ dbinom(prob[1],1)
	z[6] ~ dbinom(prob[1],1)
	z[7] ~ dbinom(prob[1],1)
	z[8] ~ dbinom(prob[2],1)
	z[9] ~ dbinom(prob[1],1)
	z[10] ~ dbinom(prob[1],1)
	z[11] ~ dbinom(prob[1],1)
	z[12] ~ dbinom(prob[3],1)
	z[13] ~ dbinom(prob[1],1)

	prob[1] ~ dunif(0,1)
	prob[2] ~ dunif(0,1)
	prob[3] ~ dunif(0,1)



}'

# inizializzazioni
set.seed(1234)
inits =  list(
  list("b_app" = runif(0,-1,1), "z" = rep(1,p), "prob"=rep(1,3)),
  list("b_app" = runif(0,-1,1), "z" = rep(1,p), "prob"=rep(1,3))
)

# parametri da salvare
SavePar = c("b", "lambda", "z","prob")


# fittiamo il modello
set.seed(123)
fit_glm3 = jags(
			data 								= dataList,
			inits 							= inits,
			parameters.to.save 	= SavePar,
			model.file 					= textConnection(mod3_string), # il model.file dovrebbe essere un file,
			 																									# "textConnection" crea il file e lo carica
			n.chains 						= nchain,
			n.iter 							= 10000,
			n.burnin 						= 5000,
			n.thin 							= 2,
			DIC 								= T # simile all'AIC, serve per il model selection
)


#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####
beta_post_chain1 = fit_glm3$BUGSoutput$sims.array[,1,1:p]
beta_post_chain2 = fit_glm3$BUGSoutput$sims.array[,2,1:p]
colnames(beta_post_chain1) = colnames(beta_post_chain2) = colnames(X)

# confrontiamo le catene per vedere se il modello è a convergenza
par(mfrow=c(3,2))
for(i in 1:p)
{
	plot(beta_post_chain1[,i], type="l", main= colnames(beta_post_chain1)[i])
	lines(beta_post_chain2[,i],col=2)
}

# stimiamo le medie e gli intervalli di credibilità
beta_mean = round(colMeans(beta_post_chain1),5)
beta_CI1  = round(apply(beta_post_chain1,2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1,2,quantile,prob=0.975),5)

Beta_Res = cbind(beta_mean,beta_CI1,beta_CI2)
Beta_Res

#### #### #### ####
#### vediamo la variabile z
z_post = fit_glm3$BUGSoutput$sims.array[,1,c("z[1]","z[8]","z[10]","z[11]","z[12]")]
colnames(z_post) = c("Season","Temperature","Rain","Snow","Trans-Time")
colMeans(z_post)

## vediamo se c'è dipendenza tra i valori di z di temperatura e time
TT = table(z_post[,"Temperature"], z_post[,"Trans-Time"])
TT
# la tabele ci dice che quanto la temperatura è "high" (seconda riga), il modello tende a preferire
# la trasformata log del tempo, invece quanto la temperature è "low", la trasformazione radice ha prob approx 0.


# per le stime di beta potremmo vederle anche condizionate al valore di z

# z trans_time
W         = z_post[,"Trans-Time"]==1
beta_mean = round(colMeans(beta_post_chain1[W,]),5)
beta_CI1  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.975),5)
Beta_Res_z1 = cbind(beta_mean,beta_CI1,beta_CI2)

W         = z_post[,"Trans-Time"]==0
beta_mean = round(colMeans(beta_post_chain1[W,]),5)
beta_CI1  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.975),5)
Beta_Res_z0 = cbind(beta_mean,beta_CI1,beta_CI2)

cbind(Beta_Res_z1,NA,Beta_Res_z0)

# z temp
W         = z_post[,"Temperature"]==1
beta_mean = round(colMeans(beta_post_chain1[W,]),5)
beta_CI1  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.975),5)
Beta_Res_z1 = cbind(beta_mean,beta_CI1,beta_CI2)

W         = z_post[,"Temperature"]==0
beta_mean = round(colMeans(beta_post_chain1[W,]),5)
beta_CI1  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1[W,],2,quantile,prob=0.975),5)
Beta_Res_z0 = cbind(beta_mean,beta_CI1,beta_CI2)

cbind(Beta_Res_z1,NA,Beta_Res_z0)

# Il valore di log temp non è significativo nel modelo generale, non lo è
# se condiziono solo quando z lo sceglie come trasformazione, ma lo è quando
# z della temparatura sceglie il low (ma questa stima si basa co solo 51 iterazioni, quindi poco affidabile)
# i.e.
sum(z_post[,"Temperature"]==0)
### ### ### ### ### ### ### ### ###
### facciamo la stessa analisi ma con una prior su pi informativa e vicino a 0
mod4_string <- 'model {
  for(i in 1:n) {
    Y[i]  ~  dpois(lambda[i])
		lambda[i] = exp(mu[i])
    mu[i] <- 	b[1]+ #(Intercept)
							b[2]*X[i,2]+ #weekdayMonday
							b[3]*X[i,3]+ #weekdaySaturday
							b[4]*X[i,4]+ #weekdaySunday
							b[5]*X[i,5]+ #weekdayThursday
							b[6]*X[i,6]+ #weekdayTuesday
							b[7]*X[i,7]+ #weekdayWednesday
							b[8]*X[i,8]+ #hightemp
							b[9]*X[i,9]+ #lowtemp
							b[10]*X[i,10]+ #precip_rain
							b[11]*X[i,11]+ #precip_snow
							b[12]*X[i,12]+  # log(time)
							b[13]*X[i,13]  # time^0.5
  }
	var1 = 10000
	var2 = 0.0000000001

	b_app[1] ~ dnorm(0,1/var1)
	b_app[2] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[3] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[4] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[5] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[6] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[7] ~ dnorm(0,z[1]*1/var1+(1-z[1])*1/var2)
	b_app[8] ~ dnorm(0,1/var1)
	b_app[9] ~ dnorm(0,1/var2)
	b_app[10] ~ dnorm(0,z[10]*1/var1+(1-z[10])*1/var2)
	b_app[11] ~ dnorm(0,z[11]*1/var1+(1-z[11])*1/var2)
	b_app[12] ~ dnorm(0,1/var1)
	b_app[13] ~ dnorm(0,1/var2)

	b[1]  	= b_app[1]
	b[2]   	= b_app[2]
	b[3]   	= b_app[3]
	b[4]   	= b_app[4]
	b[5]   	= b_app[5]
	b[6]   	= b_app[6]
	b[7]   	= b_app[7]
	b[8]   	= (z[8]*			b_app[8]+(1-z[8])*b_app[9])
	b[9]   	= (1-z[8])*	b_app[8]+z[8]*b_app[9]
	b[10]   = b_app[10]
	b[11]   = b_app[11]
	b[12]   = (z[12]*			b_app[12]+(1-z[12])*b_app[13])
	b[13]   = (1-z[12])*	b_app[12]+z[12]*b_app[13]


	z[1] ~ dbinom(prob[1],1)
	z[2] ~ dbinom(prob[1],1)
	z[3] ~ dbinom(prob[1],1)
	z[4] ~ dbinom(prob[1],1)
	z[5] ~ dbinom(prob[1],1)
	z[6] ~ dbinom(prob[1],1)
	z[7] ~ dbinom(prob[1],1)
	z[8] ~ dbinom(prob[2],1)
	z[9] ~ dbinom(prob[1],1)
	z[10] ~ dbinom(prob[1],1)
	z[11] ~ dbinom(prob[1],1)
	z[12] ~ dbinom(prob[3],1)
	z[13] ~ dbinom(prob[1],1)

	prob[1] ~ dbeta(0.0001,1)
	prob[2] ~ dunif(0.0001,1) ## preferisco low-temp
	prob[3] ~ dunif(0.0001,1) ## preferisco time^0.5



}'

# inizializzazioni
set.seed(1234)
inits =  list(
  list("b_app" = runif(0,-1,1), "z" = rep(1,p), "prob"=rep(1,3)),
  list("b_app" = runif(0,-1,1), "z" = rep(1,p), "prob"=rep(1,3))
)

# parametri da salvare
SavePar = c("b", "lambda", "z","prob" )


# fittiamo il modello
set.seed(123)
fit_glm4 = jags(
			data 								= dataList,
			inits 							= inits,
			parameters.to.save 	= SavePar,
			model.file 					= textConnection(mod4_string), # il model.file dovrebbe essere un file,
			 																									# "textConnection" crea il file e lo carica
			n.chains 						= nchain,
			n.iter 							= 10000,
			n.burnin 						= 5000,
			n.thin 							= 2,
			DIC 								= T # simile all'AIC, serve per il model selection
)


#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####
beta_post_chain1 = fit_glm4$BUGSoutput$sims.array[,1,1:p]
beta_post_chain2 = fit_glm4$BUGSoutput$sims.array[,2,1:p]
colnames(beta_post_chain1) = colnames(beta_post_chain2) = colnames(X)

# confrontiamo le catene per vedere se il modello è a convergenza
par(mfrow=c(3,2))
for(i in 1:p)
{
	plot(beta_post_chain1[,i], type="l", main= colnames(beta_post_chain1)[i])
	lines(beta_post_chain2[,i],col=2)
}

# stimiamo le medie e gli intervalli di credibilità
beta_mean = round(colMeans(beta_post_chain1),5)
beta_CI1  = round(apply(beta_post_chain1,2,quantile,prob=0.025),5)
beta_CI2  = round(apply(beta_post_chain1,2,quantile,prob=0.975),5)

Beta_Res = cbind(beta_mean,beta_CI1,beta_CI2)
Beta_Res

#### #### #### ####
#### vediamo la variabile z
z_post = fit_glm4$BUGSoutput$sims.array[,1,c("z[1]","z[8]","z[10]","z[11]","z[12]")]
colnames(z_post) = c("Season","Temperature","Rain","Snow","Trans-Time")
colMeans(z_post)
