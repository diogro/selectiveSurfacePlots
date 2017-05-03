if(!require(lattice)) { install.packages("lattice"); library(lattice) }
if(!require(mvtnorm)) { install.packages("mvtnorm"); library(mvtnorm) }
if(!require(numDeriv)) { install.packages("numDeriv"); library(numDeriv) }
if(!require(ellipse)) { install.packages("ellipse"); library(ellipse) }
if(!require(grDevices)) { install.packages("grDevices"); library(grDevices) }
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }

w_cov = matrix(c(1.0, 0.0,
                 0., 1.0), ncol = 2)

W_bar = function(x) {
  log((x[1]^2 - x[2]^2) + 101)
}

x <- seq(-10, 10, step) ## valores para mu
y <- seq(-10, 10, step)
X <- as.matrix(expand.grid(x, y))
Z <- vector()
for(i in 1:nrow(X)){
  Z[i] <- W_bar(c(X[i,1], X[i,2]))
}
Z = exp(Z - log(sum(exp(Z))))
b <- matrix(Z, length(x))

mypalette = colorRampPalette(c(wes_palette(10, name = "Zissou", type = "continuous"), "darkred"))
#mypalette = colorRampPalette(c(wes_palette(10, name = "Zissou", type = "continuous"), "darkred")[10:1])


png("~/Dropbox/Cursos/ModularidadeRibeirao2017/selecaoMultivariada/figures/saddle.png", width = 1080, height = 900)
filled.contour(x, y, z = b, color.palette = mypalette,
               plot.axes = {
                 axis(1);
                 axis(2);
               }
)
dev.off(dev.cur())


require(lattice)
g <- expand.grid(x = -10:10, y = -10:10, gr = 1:2)
g$z <- g$x^2 - g$y^2
png("~/Dropbox/Cursos/ModularidadeRibeirao2017/selecaoMultivariada/figures/saddle3D.png", width = 1080, height = 900)
wireframe(z ~ x * y, data = g, groups = gr,
          scales = list(arrows = FALSE),
          drape = TRUE, col.regions = mypalette(100),
          screen = list(z = 30, x = -60))
dev.off(dev.cur())

