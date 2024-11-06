###### STUDENT DYNAMICS ##########

if (!require(PhaseTypeR)) install.packages("PhaseTypeR")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(igraph)) install.packages("igraph")

library(PhaseTypeR)
library(ggplot2)
library(tidyr)
library(igraph)

########################
# EXAMPLE
no.gdos<-3
a<-b<-matrix(0,nrow=no.gdos)
a=c(0.6,0.8,0.9)
b=c(0.2,0.15,0.08)

beta<-(a[1]*a[2]*a[3])/((1-b[1])*(1-b[2])*(1-b[3]))
vproba=1-beta
alpha=c(0,0,1)
S=matrix(c(b[1],0,0,1-b[2],b[2],0,0,1-b[3],b[3]), nrow=3,byrow=T)
s=c(1-b[1],0,0)
  
################
# Differential equation
anos=12
y=rep(0,anos-1)
uu=50
u=rep(uu,anos)
x=matrix(0,ncol=3,nrow=anos)
for(k in 1:(anos-1))
{
  x[k+1,1]<-b[1]*x[k,1]+u[k]
  x[k+1,2]<-a[1]*x[k,1]+b[2]*x[k,2]
  x[k+1,3]<-a[2]*x[k,2]+b[3]*x[k,3]
  y[k]<-a[3]*x[k,3]
  print(y[k])
}
y

################
# We generate the DPH object: # librería PhaseTypeR
dph <- DPH(S,alpha) 
dph

net <- phase_type_to_network(dph)
plot(net, edge.curved=.1, edge.label = E(net)$weight,
     edge.color = ifelse(as_data_frame(net, what="edges")$from == 'V0',
                         'purple', 'grey'),
     layout = layout_with_fr(net,  weights = rep(1, length(E(net)))))
#
x <- seq(0,10)
pdf <- dDPH(x, dph) # densidad  pmf
cdf <- pDPH(x, dph) # cumulative


# Plot PDF and CDF
df <- data.frame(x = x, PMF = pdf, CDF = cdf)
ggplot(df, aes(x)) +
  geom_step(aes(y = CDF, color = "CDF"), size = 1) +
  geom_point(aes(y = PMF, color = "PMF"), size = 2) +
  labs(title = "", 
       x = "x", y = "Probability") +
  scale_x_continuous(breaks = x) +
  scale_color_manual(name = "Function", 
                     values = c("PMF" = "blue", "CDF" = "red"),
                     labels = c("CDF", "PMF")) +
  theme(legend.position = "right")  

################
# Random sample
set.seed(2018)
tam=1000
muestra<-rDPH(tam, dph)

# Plots random sample with ggplot and shows density
df_sample <- data.frame(x = muestra)
ggplot(df_sample, aes(x)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(color = "blue") +
  labs(title = "Random Sample from DPH",
       x = "Values", y = "Density") +
  scale_x_continuous(breaks = x)

# Mean
mean(muestra)

#########################
# Data frame
Tabla<-data.frame(0:10,round(pdf*100,1)) # as a percentage
colnames(Tabla)<-c("Año","Porcentaje")
Tabla

# Probabilities to return to the case (beta, 0, 0)
probBeta<-pdf*beta 
probBeta[1]<-vproba
sum(probBeta) # approx 1
TablaBeta<-data.frame(0:10,round(probBeta*100,1))
colnames(TablaBeta)<-c("Año","Porcentaje")
TablaBeta
#sum(TablaBeta$Porcentaje)                      

#########################
# data frame
dfproba <- data.frame(x = 0:10,
                      pmf_IniProb = pdf,  # Datos de ejemplo para pdf
                      pmf = probBeta)  # Datos de ejemplo para probBeta
# ggplot
dfproba_long <- pivot_longer(dfproba, cols = c(pmf_IniProb, pmf), 
                             names_to = "Function", values_to = "Probability")
ggplot(dfproba_long, aes(x = factor(x), y = Probability, fill = Function)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Comparison of probabilities",
       x = "Years",
       y = "Probability",
       fill = "Function") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#########################
# Comparison
FFn=cdf*beta*uu

dff <- data.frame(x = x, y = y, CDFN = FFn)
ggplot(dff, aes(x)) +
  geom_step(aes(y = CDFN, color = "CDFN"), size = 1)  +
  geom_point(aes(y =y, color = "y") , size = 2)+
  labs(title = "", x = "Year", y = "Number of graduates") +
  scale_x_continuous(breaks = x) +
  scale_color_manual(name = "", 
                     values = c("y" = "black", "CDFN" = "orange"),
                     labels = c("y(k) DPH", "y(k)")) +
  theme(legend.position = "right")  
  



