##### SUPPLY CHAIN ######
if (!require(PhaseTypeR)) install.packages("PhaseTypeR")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(igraph)) install.packages("igraph")

library(PhaseTypeR)
library(ggplot2)
library(tidyr)
library(igraph)
########################
# Example
no.gdos<-3
a<-b<-d<-c(0)
d=c(0.15,0.08)
a=c(0.6,0.8)
b=0.05
g=0.8

# Parameters system
A=matrix(c(1-a[1]-d[1],0,0,a[1],1-a[2]-d[2],b,0,a[2],1-b-g),nrow=no.gdos,byrow=T);A
B=matrix(c(1,0,0),nrow=no.gdos,byrow=T);B
C=matrix(c(0,0,g),ncol=no.gdos);C

# Parameters DPH
x<-solve(diag(no.gdos)-A)%*%B;x
M<-diag(c(x));M
# DPH
alpha0=C%*%M;alpha0
S=solve(M)%*%A%*%M;S
s=solve(M)%*%B;s
#(diag(no.gdos)-S)%*%rep(1,no.gdos)
  
beta<-alpha0[3];beta
vproba=1-beta;vproba
alpha=alpha0/beta;alpha ## we re-define: probability vector

################
# Differential equation
anos=15
y=rep(0,anos-1)
uu=100
u=rep(uu,anos)
x=matrix(0,ncol=3,nrow=anos)
for(k in 1:(anos-1))
{
  x[k+1,1]<-(1-a[1]-d[1])*x[k,1]+u[k]
  x[k+1,2]<-a[1]*x[k,1]+(1-a[2]-d[2])*x[k,2]+b*x[k,3]
  x[k+1,3]<-a[2]*x[k,2]+(1-b-g)*x[k,3]
  y[k]<-g*x[k,3]
  print(y[k])
}
y

################
# We generate the DPH object: # librería PhaseTypeR
dph <- DPH(S,alpha)#dph <- DPH(round(S,2),alpha)  
dph

net <- phase_type_to_network(dph)
plot(net, edge.curved=.1, edge.label = E(net)$weight,
     edge.color = ifelse(as_data_frame(net, what="edges")$from == 'V0',
                         'purple', 'grey'),
     layout = layout_with_fr(net,  weights = rep(1, length(E(net)))))


x <- seq(0,13)
pdf <- dDPH(x, dph) # densidad  pmf
cdf <- pDPH(x, dph) # cumulative

# Plot PDF and CDF
df <- data.frame(x = x, PMF = pdf, CDF = cdf)
ggplot(df, aes(x)) +
  geom_step(aes(y = CDF, color = "CDF"), size = 1) +
  geom_point(aes(y = PMF, color = "PMF"), size = 2) +
  labs(title = "", # Eliminamos el título
       x = "x", y = "Probability") +
  scale_x_continuous(breaks = x) +
  scale_color_manual(name = "Function", 
                     values = c("PMF" = "blue", "CDF" = "red"),
                     labels = c("CDF", "PMF")) +
  theme(legend.position = "right")  # Posición de la leyenda

################
# random sample
set.seed(0)
tam=1000
muestra<-rDPH(tam, dph)
# Mean
mean(muestra)

Tabla<-data.frame(x,round(pdf*100,1)) # as a percentage
colnames(Tabla)<-c("Año","Porcentaje")
Tabla

# Probabilities to return to the case (beta, 0, 0)
probBeta<-pdf*beta 
probBeta[1]<-vproba
sum(probBeta) #has to be 1
TablaBeta<-data.frame(x,round(probBeta*100,1))
colnames(TablaBeta)<-c("Año","Porcentaje")
TablaBeta
#sum(TablaBeta$Porcentaje)                      

# Data frame
dfproba <- data.frame(x = x,
                      pmf_IniProb = pdf,  # pdf
                      pmf = probBeta)  # probBeta
# For ggplot
dfproba_long <- pivot_longer(dfproba, cols = c(pmf_IniProb, pmf), 
                             names_to = "Function", values_to = "Probability")
# Plot
ggplot(dfproba_long, aes(x = factor(x), y = Probability, fill = Function)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "",
       x = "Month",
       y = "Probability",
       fill = "Function") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#########################
# Comparison
FFn=cdf*beta*uu

# Plot
dff <- data.frame(x = x, y = y, CDFN = FFn)
ggplot(dff, aes(x)) +
  geom_step(aes(y = CDFN, color = "CDFN"), size = 1)  +
  geom_point(aes(y =y, color = "y") , size = 2)+
  labs(title = "", x = "Month", y = "Products sold to customers") +
  scale_x_continuous(breaks = x) +
  scale_color_manual(name = "", 
                     values = c("y" = "black", "CDFN" = "orange"),
                     labels = c("y(k) DPH", "y(k)")) +
  theme(legend.position = "right")  


