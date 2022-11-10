### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Si vuole studiare come gli studenti decidano 
### il tipo di programma da seguire
### tra generale, accademico e tecnico (vocation)
### 
### ses: stato economico
### schtyp: tipo di scuola
### read, write, math e science: valutazione nelle rispettive materie
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

## 
rm(list =ls())
## 

## librerie
require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
##

# DIrectory
DIR = "C:/Users/marti/OneDrive/Desktop/modelli_statistici/LEZIONE 9"
setwd(DIR)
# carichiamo i dati
load('DataBinomial.RData')
summary(Data)
head(Data)


## ## ## ## ## ## ## ## 
## Analisi Descrittive
## ## ## ## ## ## ## ## 
?with
with(Data, table(ses, prog))
table(Data$ses, Data$prog)

?do.call
?with
?tapply
with(Data, do.call(rbind, tapply(write, prog, function(x) c(M = mean(x), SD = sd(x)))))
mean(Data$write[Data$prog=='general'])
sd(Data$write[Data$prog=='general'])


## Qualche boxplot
boxplot(Data$write~ Data$prog)
quantile(Data$write[Data$prog=='general'],prob=c(0.25,0.5,0.75))
boxplot(Data$science~ Data$prog)
boxplot(Data$math~ Data$prog)


## volendo possiamo testare delle ipotesi con il t-test o equivalentemente
## con un nodello lineare o anova

t.test(Data$write[Data$prog=="general"],Data$write[Data$prog=="academic"],var.equal=T)
summary(lm(write~ prog, data =Data[Data$prog!="vocation",] ))
summary(aov(write ~ prog, data=Data[Data$prog!="vocation",]))


### ### ### ### ### ### ### 
### testiamo i modelli nominali
### ### ### ### ### ### ### 

# attenzione che le x devono essere indipendenti
plot(Data$math,Data$science)
plot(Data$math,Data$write)
plot(Data$science,Data$write)
plot(Data$read,Data$write)

Data$V1V3 = ifelse(Data$prog=="vocation",1,0)
Mod13 = glm(V1V3 ~ ses + schtyp +read+ write+ math, data =Data[Data$prog!="academic",], family="binomial" )
summary(Mod13) # vocation=1, general = 0



Data$V1V2 = ifelse(Data$prog=="academic",1,0)
Mod12 = glm(V1V2 ~ ses + schtyp +read+ write+ math, data =Data[Data$prog!="vocation",], family="binomial" )
summary(Mod12) # academic=1, general = 0



### VEdiamo i parametri
cbind(summary(Mod13)$coefficients[,c(1,4)],summary(Mod12)$coefficients[,c(1,4)])

### Scegliamo i subset dei parametri con step
?step
StepMod13 = step(Mod13)
summary(StepMod13)

StepMod12 = step(Mod12)
summary(StepMod12)

## Possiamo calcolare le previsioni
Exp13 = exp(predict(StepMod13, type="link",newdata=Data))
Exp12 = exp(predict(StepMod12, type="link",newdata=Data))


pigeneral = 1/(1+Exp13+Exp12) # GENERAL
pivocation = Exp13/(1+Exp13+Exp12)# VOCATION
piacademic = Exp12/(1+Exp13+Exp12)# ACADEMIC

Prev = cbind(pigeneral,pivocation,piacademic)
rowSums(Prev) 

## Qualche interpretazione grafica

## write
plot(0,0, xlim=c(min(Data$write), max(Data$write)), ylim=c(0,1))
points(Data$write,Prev[,1], col=2, pch=20)
points(Data$write,Prev[,2], col=3, pch=20)
points(Data$write,Prev[,3], col=4, pch=20)

## math
plot(0,0, xlim=c(min(Data$math), max(Data$math)), ylim=c(0,1))
points(Data$math,Prev[,1], col=2, pch=20)
points(Data$math,Prev[,2], col=3, pch=20)
points(Data$math,Prev[,3], col=4, pch=20)


par(mfrow=c(1,3))
boxplot(Prev[,1] ~ Data$ses, col=1:3)
boxplot(Prev[,2] ~ Data$ses, col=1:3)
boxplot(Prev[,3] ~ Data$ses, col=1:3)

par(mfrow=c(1,1))


### ### ### ### ### ### ### 
### Possiamo testare altri tipi di modelli
### ### ### ### ### ### ### 

Data$V1 = ifelse(Data$prog=="general",1,0)

ModV1 = glm(V1 ~ ses + schtyp +read+ write+ math, data =Data, family="binomial" )
summary(ModV1) # vocation=1, general = 0

Prev=ModV1$fitted.values
predict(ModV1,newdata=Data)
predict(ModV1,newdata=Data,type='response')


Data$V2 = ifelse(Data$prog=="general",1,0)
ModV2 = glm(V2 ~ ses +read+ write+ math, data =Data, family="binomial" )
summary(ModV1) # vocation=1, general = 0
exp(coef(ModV2))

testdata=data.frame(ses=c('low','middle','high'),write=mean(Data$write),math=mean(Data$math),read=mean(Data$read))
testdata$prob=predict(ModV2,newdata=testdata,type='response')
testdata


testdata=data.frame(ses=c('low'),write=c(30,40,50),math=mean(Data$math),read=mean(Data$read))
testdata$prob=predict(ModV2,newdata=testdata,type='response')
testdata
