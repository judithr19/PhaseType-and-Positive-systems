##### CONTINUOUS SYSTEMS-CPH ######

# Install packages
if (!require(PhaseTypeR)) install.packages("PhaseTypeR")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(igraph)) install.packages("igraph")
if (!requireNamespace("deSolve", quietly = TRUE)) {install.packages("deSolve")}

library(PhaseTypeR)
library(ggplot2)
library(tidyr)
library(igraph)
library(deSolve)

########################
# Dimension
n<-3
# Parameters of the system
A=matrix(c(-2,1,0,0,-1,1,0,0,-1),byrow=T,nrow=n);
B=matrix(c(1,1,1),byrow=T,nrow=n);
C=matrix(c(1,0,0),ncol=n);

# Define
Delta1<-rbind(cbind(A,B),rep(0,n+1));Delta1
uu<-50
u=rep(uu,n)

#################################################
# Realization of the system: EDO
# Define the function
sistema <- function(t, state, parms) {
  dxdt <- A %*% state + B*u
  return(list(dxdt))
}
# Define the times
tiempos <- seq(0, 10, by = 0.1)
# Define the initial contition
x0 <- rep(0,n)#c(0, 0, 0)
# Solve the system
solucion <- ode(y = x0, times = tiempos, func = sistema, parms = NULL)
solucion
# Output y(t)
y <- as.vector(t(C %*% t(solucion[, -1])))
y
# Plot the output y(t)
plot(tiempos, y, type = "l", xlab = "Time", ylab = "y(t)", main = "Output of the system")

data.frame(tiempos,y)
#################################################
# Relation POSITIVE SYSTEM-CPH
# Parameters of the system
A;B;C

# Eigenvalues de A
a<-eigen(A);a
eta=max(abs(a$values));eta
Delta<-Delta1+eta*diag(n+1);Delta # \Tilde{A}
# Parameters of the CPH
v <- eigen(Delta);v
# Find the index of the largest eigenvalue
vmax=max(v$values)
indice_max <- which.max(v$values)
# Obtain the eigenvector associated with the largest eigenvalue
eigenvector_max <- v$vectors[, indice_max]
U=diag(c(eigenvector_max[1:n]))/eigenvector_max[n+1];U
  
# CPH
alpha0=C%*%U;alpha0
# Constructing the probability vector by normalizing the original vector
vector_probabilidad <- alpha0 / sum(alpha0)
print(vector_probabilidad)
S=solve(U)%*%A%*%U;S
s=solve(U)%*%B;s
#-S%*%rep(1,n)

beta<-alpha0[1];beta
alpha=vector_probabilidad # we re-define: probability vector

############################################
# Generate the DPH object: # PhaseTypeR library
ph <- PH(S,alpha) 
ph

net <- phase_type_to_network(ph)
plot(net, edge.curved=.1, edge.label = E(net)$weight,#vertex.size = 20,
     edge.color = ifelse(as_data_frame(net, what="edges")$from == 'V0',
                         'purple', 'grey'),
     layout = layout_with_fr(net,  weights = rep(1, length(E(net)))))

par(mfrow=c(3,3))
for (i in 0:8) {
  net <- phase_type_to_network(ph, i)
  plot(net, edge.arrow.size=0.5, edge.curved=.1, 
       edge.label = E(net)$weight, 
       main=paste0('Time = ', i),
       edge.color = ifelse(as_data_frame(net, what="edges")$from == 'V0',
                           'purple', 'grey'),
       layout = layout_with_fr(net,  weights = rep(1, length(E(net)))))
}

t <- seq(0, 10, by = 0.1)
pdf <- dPH(t, ph) # densidad  pmf
cdf <- pPH(t, ph) # cumulative

# Graphics la PDF y la CDF
df <- data.frame(x = t, PDF = pdf, CDF = cdf)
dev.new(4,4)
ggplot(df, aes(t)) +
  geom_step(aes(y = CDF, color = "CDF"), linewidth = 1) +
  geom_line(aes(y = PDF, color = "PDF"), size = 2) +
  labs(title = "CPH", # Eliminamos el título
       x = "Time", y = "Probability") +
  scale_x_continuous(breaks = t*10) +
  scale_color_manual(name = "Function", 
                     values = c("CDF" = "red", "PDF" = "blue"),
                     labels = c("CDF", "PDF")) +
  theme(legend.position = "right")  

#####################################
# COMPARISON OF DPH AND DIFFERENTIAL EQ. DIFFERENTIAL
FFn=cdf*beta*uu

# The case where 0<beta<1:
if(beta<1)
{
  FFn=(cdf-(1-beta))*uu/pdf[1]
}

# Gráfica 
dff <- data.frame(x = t, y = y, CDFN = FFn)
#dev.new(4,4)
ggplot(dff, aes(x)) +
  geom_step(aes(y = CDFN, color = "CDFN"), size = 1)  +
  geom_point(aes(y =y, color = "y") , size = 2)+
  labs(title = "", x = "Time", y = "Output") +
  scale_x_continuous(breaks = t*10) +
  scale_color_manual(name = "", 
                     values = c("y" = "black", "CDFN" = "orange"),
                     labels = c("y(t) CPH", "y(t)")) +
  theme(legend.position = "right")  

###################################
# random sample
set.seed(2024)
tam=1000
muestra<-rPH(tam, ph)
# Plots random sample with ggplot and shows density
df_sample <- data.frame(x = muestra)
#dev.new(4,4)
ggplot(df_sample, aes(x)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_density(color = "blue") +
  labs(title = "Random Sample from CPH",
       x = "t", y = "Density") +
  scale_x_continuous(breaks = t*5)

###################################








