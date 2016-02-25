library(lattice)
library(mvtnorm)
library(numDeriv)
library(ellipse)
library(grDevices)
library(wesanderson)


plotTrajectory <- function (start) {
  current_gen = start
  polygon(ellipse(0.8, centre = current_gen, level = 0.3), col = "darkgray")
  segments(current_gen[1] - v1[1], current_gen[2] - v1[2], 
           current_gen[1] + v1[1], current_gen[2] + v1[2], lwd = 2) 
  segments(current_gen[1] - v2[1], current_gen[2] - v2[2], 
           current_gen[1] + v2[1], current_gen[2] + v2[2], lwd = 2) 
  
  points(current_gen[1], current_gen[2], pch = 19)
  net_beta = c(0, 0)
  for(i in 1:gen){
    beta = grad(W_bar, t(current_gen))
    net_beta = net_beta + beta
    next_gen = current_gen + G%*%beta
    arrows(current_gen[1], current_gen[2], 
           next_gen[1], next_gen[2], pch = 18, length = 0.14, lwd = 2.5)
    arrows(current_gen[1], current_gen[2], 
           current_gen[1] + beta[1]/5, current_gen[2] + beta[2]/5, 
           pch = 18, length = 0.14, lwd = 2.5, col = 'purple')
    current_gen = next_gen
  }
  arrows(start[1], start[2], 
         start[1] + net_beta[1]/5, start[2] + net_beta[2]/5, 
         pch = 18, length = 0.14, lwd = 2.5, col = 'red')
#  net_delta = G %*% net_beta
#  arrows(start[1], start[2], 
#         start[1] + net_delta[1], start[2] + net_delta[2], 
#         pch = 18, length = 0.14, lwd = 2.5, col = 'black')
}

w_cov = matrix(c(1, 0.7, 0.7, 1), ncol = 2)
G = matrix(c(1, 0.8, 0.8, 1)/2, ncol = 2)
cov2cor(G)
W_bar = function(x) {
    log(
    dmvnorm(x, mean = c(4, 5), sigma = w_cov) + 
    dmvnorm(x, mean = c(1, 5), sigma = w_cov) + 
    dmvnorm(x, mean = c(2, 2), sigma = w_cov) + 
    dmvnorm(x, mean = c(7, 5), sigma = w_cov) + 
    dmvnorm(x, mean = c(5, 2), sigma = w_cov))
}
step = 0.1
x <- seq(-1.5, 8.5, step) ## valores para mu
y <- seq(-1.5, 8.5, step)
X <- as.matrix( expand.grid(x, y))
colnames(X) <- c("mu","var")
Z <- vector()
for(i in 1:nrow(X)){
  Z[i] <- W_bar(c(X[i,1], X[i,2]))
}
Z = exp(Z - log(sum(exp(Z))))
b <- matrix(Z, length(x))

grad = grad(W_bar, c(5, 5)) 
gen = 10
es = eigen(cov2cor(G))$values
v1 = sqrt(es[1])/1.2 * eigen(cov2cor(G))$vectors[,1]
v2 = sqrt(es[2])/1.2 * eigen(cov2cor(G))$vectors[,2]

mypalette = colorRampPalette(wes_palette(10, name = "Zissou", type = "continuous"))
png("multipeaklandscape.png", width = 1000, height = 900)
filled.contour(x, y, z = b, color.palette = mypalette,
               plot.axes = { 
                 axis(1); 
                 axis(2);
                 plotTrajectory(c(1.8,4))
                 plotTrajectory(c(1.5,0))
                 plotTrajectory(c(5,4))
               }
)
dev.off()


#################################33
# Single Peak
#################################


W_bar = function(x) {
  log(
    dmvnorm(x, mean = c(3, 3), sigma = w_cov))
}

x <- seq(0, 6, 0.1) ## valores para mu
y <- seq(-1, 5, 0.1)
X <- as.matrix( expand.grid(x, y))
colnames(X) <- c("mu","var")
Z <- vector()
for(i in 1:nrow(X)){
  Z[i] <- W_bar(c(X[i,1], X[i,2]))
}
Z = exp(Z - log(sum(exp(Z))))
b <- matrix(Z, length(x))

es = eigen(cov2cor(G))$values
v1 = sqrt(es[1])/1.2 * eigen(cov2cor(G))$vectors[,1]
v2 = sqrt(es[2])/1.2 * eigen(cov2cor(G))$vectors[,2]

png("singlepeaklandscape.png", width = 1000, height = 900)
filled.contour(x, y, z = b, color.palette = mypalette,
               plot.axes = { 
                 axis(1); 
                 axis(2);
                 plotTrajectory(c(1.8,4))
                 plotTrajectory(c(1.5,0))
                 plotTrajectory(c(5,4))
               }
)
dev.off()
