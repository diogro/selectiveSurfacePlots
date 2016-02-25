library(lattice)
library(mvtnorm)
library(numDeriv)
library(ellipse)

w_cov = matrix(c(1, 0.7, 0.7, 1), ncol = 2)
G = matrix(c(1, 0.8, 0.8, 1)/2, ncol = 2)

W_bar = function(x) {
    log(
    dmvnorm(x, mean = c(4, 5), sigma = w_cov) + 
    dmvnorm(x, mean = c(1, 5), sigma = w_cov) + 
    dmvnorm(x, mean = c(2, 2), sigma = w_cov) + 
    dmvnorm(x, mean = c(7, 5), sigma = w_cov) + 
    dmvnorm(x, mean = c(5, 2), sigma = w_cov))
}

x <- seq(-1, 9, 0.1) ## valores para mu
y <- seq(-1, 8, 0.1)
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

png("multipeaklandscape.png", width = 948, height = 567)
filled.contour(x, y, z = b,color = terrain.colors,
               plot.axes = { 
                 axis(1); 
                 axis(2);
                 current_gen = c(1.8,4)
                 polygon(ellipse(0.8, centre = current_gen, level = 0.3), type = 'l', col = "darkgray")
                 points(current_gen[1], current_gen[2], pch = 19)
                 for(i in 1:gen){
                   next_gen = current_gen + G%*%grad(W_bar, t(current_gen))
                   arrows(current_gen[1], current_gen[2], next_gen[1], next_gen[2], pch = 18, length = 0.14)
                   current_gen = next_gen
                 }
                 current_gen = c(1.5,0)
                 polygon(ellipse(0.8, centre = current_gen, level = 0.3), type = 'l', col = "darkgray")
                 points(current_gen[1], current_gen[2], pch = 19)
                 for(i in 1:gen){
                   next_gen = current_gen + G%*%grad(W_bar, t(current_gen))
                   arrows(current_gen[1], current_gen[2], next_gen[1], next_gen[2], pch = 18, length = 0.14)
                   current_gen = next_gen
                 }
                 current_gen = c(5,4)
                 polygon(ellipse(0.8, centre = current_gen, level = 0.3), type = 'l', col = "darkgray")
                 points(current_gen[1], current_gen[2], pch = 19)
                 for(i in 1:gen){
                   next_gen = current_gen + G%*%grad(W_bar, t(current_gen))
                   arrows(current_gen[1], current_gen[2], next_gen[1], next_gen[2], pch = 18, length = 0.14)
                   current_gen = next_gen
                 }
               }
)
dev.off()
