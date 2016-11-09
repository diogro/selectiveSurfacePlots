if(!require(lattice)) { install.packages("lattice"); library(lattice) }
if(!require(mvtnorm)) { install.packages("mvtnorm"); library(mvtnorm) }
if(!require(ellipse)) { install.packages("ellipse"); library(ellipse) }
if(!require(grDevices)) { install.packages("grDevices"); library(grDevices) }
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }

G = matrix(c(1.0, 0.8,
             0.8, 1.0), ncol = 2)
plotEllipse <- function(center, s, color, G) {
  es = eigen(cov2cor(G))$values
  v1 = s*sqrt(es[1])/1.185 * eigen(cov2cor(G))$vectors[,1]
  v2 = s*sqrt(es[2])/1.185 * eigen(cov2cor(G))$vectors[,2]
  polygon(ellipse(G[1,2], scale = c(s, s), centre = center, level = 0.3), 
          col = color)
  segments(center[1] - v1[1], center[2] - v1[2],
           center[1] + v1[1], center[2] + v1[2], lwd = 2)
  segments(center[1] - v2[1], center[2] - v2[2],
           center[1] + v2[1], center[2] + v2[2], lwd = 2)
  points(center[1], center[2], pch = 19)
}

png("drift.png", width = 1080, height = 900)
plot(1, type="n", axes=F, xlab="", ylab="", 
     xlim=c(-3, 10), 
     ylim=c(-3, 10))
center = c(5, 5)
s = 5
color = adjustcolor(wes_palette("Royal1")[2], alpha.f = 0.5)
plotEllipse(center, 5, color, G)
for(i in 1:7){
  cov_shift = rnorm(1, 0, 0.05)
  drift = matrix(c(0.0, cov_shift,
                   cov_shift, 0.0), ncol = 2)
  plotEllipse(rmvnorm(1, center, G*4), 1, wes_palette("Royal1")[1], G + drift)
}
new_center = center
for(j in 1:8){
  shift = rmvnorm(1, c(0, 0), G/2)
  segments(new_center[1], new_center[2],
           new_center[1] + shift[1], new_center[2] + shift[2], lwd = 2, 
           col = "tomato3")
  new_center = new_center + shift
  
}
cov_shift = rnorm(1, 0, 0.05)
drift = matrix(c(0.0, cov_shift,
                 cov_shift, 0.0), ncol = 2)
plotEllipse(new_center, 1, wes_palette("Royal1")[1], G + drift)
plotEllipse(center, 1, wes_palette("Royal1")[3], G)
dev.off(dev.cur())
