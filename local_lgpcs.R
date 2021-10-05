rm(list = ls())

#load packages

#spatstat contains the main functions for dealing with spatial point processes
library(spatstat.core) #envelope(), berman.test(), rhohat(), ssf(), ppm(), Smooth(), simulate.ppm(), simulate.profilepl()
library(spatstat.geom) #owin(), distfun(), as.ppp(), as.im(), eval.im()

#spatstat.local contains the main functions to perform local fitting and diagnostics
library(spatstat.local) #bw.locppm(), locppm(), ttestmap(), homtestmap(), bw.loccit(), loccit()

#packages for plotting
library(mapproj) #map()
library(colorspace) #heat_hcl(), diverge_hcl(), terrain_hcl()
library(latex2exp) #TeX()

#load data
load("dati.RData")

ls()

#
# Descriptions of loaded objects:
#
# fz:                 distance from the faults_italy_image
# min_dist:           distance from the nearest seismic station_italy_image
# ppp_gr:             earthquakes_greece_planar point pattern
# ppp_v:              volcanoes_greece_planar point pattern
# psp_cs:             faults_greece_planar line segment pattern
# psp_pl:             plate boundary_greece_planar point pattern
# ppp_it:             earthquakes_italy_planar point pattern
#

#define colors and curve levels limits for graphs
g0 <- rev(heat_hcl(100, c = c(80, 30), l = c(30, 90), power = c(1 / 5, 1.3)))
g025 <- colourmap(rev(heat_hcl(33, c = c(80, 30), l = c(30, 90), power = c(1 / 5, 1.3))), range = c(0, 0.25)) 
g05 <- colourmap(rev(heat_hcl(33, c = c(80, 30), l = c(30, 90),  power = c(1 / 5, 1.3))), range = c(0, 0.5)) 
g06 <- colourmap(rev(heat_hcl(33, c = c(80, 30), l = c(30, 90), power = c(1 / 5, 1.3))), range = c(0, 0.06))
col_div <- colourmap(diverge_hcl(100, palette = "Purple-Brown"), range = c( - 4, 4)) 
plevels <- c(0.0001, 0.001, 0.01, 0.05)


#define windows
#greece
w_gr <- owin(poly = list(x = c(23.5, 28, 28, 20, 20), y = c(33.5, 33.5, 40.5, 40.5, 36)))
#greece: Ithaki, Kefalonia and Zakynthos
w_is <- owin(c(20, 21), c(37.5, 38.5))
#italy: Abruzzo region
w_ab <- owin(xrange = c(13, 14.8), yrange = c(41.6, 42.9))




###############################################
# Section 6. Analysis of the Greek seismicity #
###############################################



#compute the envelopes of the homogeneous K-funcion 
E <- envelope(ppp_gr, Kest, nsim = 39, rank = 1) 

#Figure1
par(mfrow = c(1, 2))
#a
map("world", xlim = c(20.02, 27.98), ylim = c(33.75, 40.45), border = NA, fill = TRUE, col = "lightgrey")
mtext(c("Latitude", "Longitude"), side = c(1, 2), line = 2.5, cex = 1.5)
plot(ppp_gr, add = TRUE, pch = "+", xlab = "Latitude", ylab = "Longitude")
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)
plot(ppp_v, add = TRUE, col = "black", bg = "green", pch = 21)
plot(psp_cs, add = TRUE, col = "blue")
plot(psp_pl, add = TRUE, col = "red")
#b
plot(E, legend = F, main = "", cex.axis = 1.5, xlab = "", ylab = "") 
mtext(c("r", "K(r)"), side = c(1, 2), line = 2.5, cex = 2)


#compute the distance variables 
d_pl <- distfun(psp_pl)
d_cs <- distfun(psp_cs)
d_v <- distfun(ppp_v)


#Figure2
par(mfrow = c(1, 3))
#a
plot(d_pl, main = "", ribargs = list(cex.axis = 1.5), col = terrain_hcl(100))
mtext(c("Latitude", "Longitude"), side = c(1, 2), line = 2.5, cex = 2)
axis(1, cex.axis = 2)
axis(2, cex.axis = 2)
#b
plot(d_cs, main = "", ribargs = list(cex.axis = 1.5), col = terrain_hcl(100))
mtext(c("Latitude", "Longitude"), side = c(1, 2), line = 2.5, cex = 2)
axis(1, cex.axis = 2)
axis(2, cex.axis = 2)
#c
plot(d_v, main = "", ribargs = list(cex.axis = 1.5), col = terrain_hcl(100))
mtext(c("Latitude", "Longitude"), side = c(1, 2), line = 2.5, cex = 2)
axis(1, cex.axis = 2)
axis(2, cex.axis = 2)

#compute the berman tests
ber_v <- berman.test(ppp_gr, d_v)
ber_cs <- berman.test(ppp_gr, d_cs) 
ber_pl <- berman.test(ppp_gr, d_pl)

par(mfrow = c(1, 3))
plot(ber_v)
plot(ber_cs)
plot(ber_pl)

#compute the smooth functions
#Figure3
par(mfrow = c(1, 3))
#a
plot(rhohat(ppp_gr, d_pl), main = "", legend = FALSE, xlab = "", ylab = "", cex.axis = 1.5)
mtext(c(TeX("$D_{pb}(\\textbf{u})$"), TeX("$f(D_{pb}(\\textbf{u}))$")), side = c(1, 2), line = 2.5, cex = 1.5)
#b
plot(rhohat(ppp_gr, d_v), main = "", legend = FALSE, xlab = "", ylab = "", cex.axis = 1.5)
mtext(c(TeX("$D_{v}(\\textbf{u})$"), TeX("$f(D_{v}(\\textbf{u}))$")), side = c(1, 2), line = 2.5, cex = 1.5)
#c
plot(rhohat(ppp_gr, d_cs), main = "", legend = FALSE, xlab = "", ylab = "", cex.axis = 1.5)
mtext(c(TeX("$D_{f}(\\textbf{u})$"), TeX("$f(D_{f}(\\textbf{u}))$")), side = c(1, 2), line = 2.5, cex = 1.5)

#select bandwidth for the local composite likelihood
sigma_gr <- bw.locppm(ppp_gr)

#Fit local model
loc_gr <- locppm(ppp_gr ~ d_cs + d_pl + d_v, sigma = sigma_gr, vcalc = "full") 

#compute local tests 
t_cs <- ttestmap(loc_gr, "d_cs")
t_pl <- ttestmap(loc_gr, "d_pl")
t_v <- ttestmap(loc_gr, "d_v")

#Figure4
par(mfrow = c(2, 3))
#a
plot(loc_gr, which = 2, main = "", ribside = "bottom", cex.axis = 1.5, col = g0)
contour(loc_gr, which = 2, add = TRUE)
plot(ppp_gr, add = TRUE, pch = "+", cex = 0.5)
#b
plot(loc_gr, which = 3, main = "", ribside = "bottom", cex.axis = 1.5, col = g0)
contour(loc_gr, which = 3, add = TRUE)
plot(ppp_gr, add = TRUE, pch = "+", cex = 0.5)
#c
plot(loc_gr, which = 4, main = "", ribside = "bottom", cex.axis = 1.5, col = g0)
contour(loc_gr, which = 4, add = TRUE)
plot(ppp_gr, add = TRUE, pch = "+", cex = 0.5)
#d
plot(t_cs, main = "", ribside = "bottom", cex.axis = 1.5, col = g0)
contour(t_cs, add = TRUE, levels = c( - 2, 2))
plot(ppp_gr, add = TRUE, pch = "+", cex = 0.5)
#e
plot(t_pl, main = "", ribside = "bottom", cex.axis = 1.5, col = g0)
contour(t_pl, add = TRUE, levels = c( - 2, 2))
plot(ppp_gr, add = TRUE, pch = "+", cex = 0.5)
#f
plot(t_v, main = "", ribside = "bottom", cex.axis = 1.5, col = g0)
contour(t_v, add = TRUE, levels = c( - 2, 2))
plot(ppp_gr, add = TRUE, pch = "+", cex = 0.5)


#compute the 'homogeneity' test
hom_loc_gr <- homtestmap(loc_gr, what = "pvalue")

#create the grid
hom_grid <- as.ppp(cbind(hom_loc_gr$x, hom_loc_gr$y), W = owin(range(hom_loc_gr$x), range(hom_loc_gr$y)))
npoints(hom_grid)

#Perform bootstrap procedure

local_boot <- function(mod, n_boot){
  
  datfr <- cbind(mod$x, mod$y)
  out <- matrix(NA, ncol = n_boot, nrow = 92)
  
  for (i in 1:n_boot){
    
    sam <- datfr[sample.int(nrow(datfr), replace = T), ]
    pps <- as.ppp(cbind(sam[, 1], sam[, 2]), W = owin(range(sam[, 1]), range(sam[, 2])))[w_gr]
    fit <- locppm(pps ~ d_cs + d_pl + d_v, sigma = bw.locppm(pps), vcalc = "full")
    hom <- homtestmap(fit, what = "pvalue")
    out[, i] <- hom$marks
    
  }
  
  return(out)
  
}

res <- local_boot(mod = ppp_gr, n_boot = 100)

hom_mean <- apply(res, 1, mean)
hom_sd <- apply(res, 1, sd)


#Figure 5
par(mfrow = c(1, 2))
#a
plot(ssf(hom_grid, hom_mean), col = g06, main = "", ribside = "bottom", cex.axis = 1.5)
contour(ssf(hom_grid, hom_mean), levels = plevels, add = TRUE)
#b
plot(ssf(hom_grid, hom_sd), col = g06, main = "", ribside = "bottom", cex.axis = 1.5)
contour(ssf(hom_grid, hom_sd), levels = plevels, add = TRUE)


#Fit nested models
loc_nest_1 <- locppm(ppp_gr ~ 1, sigma = bw.locppm(ppp_gr), vcalc = "full") 
loc_nest_2 <- locppm(ppp_gr ~ d_cs, sigma = bw.locppm(ppp_gr), vcalc = "full") 
loc_nest_3 <- locppm(ppp_gr ~ d_cs + d_pl, sigma = bw.locppm(ppp_gr), vcalc = "full") 

#perform 'homogeneity' tests
hom_nest_1 <- homtestmap(loc_nest_1, what = "pvalue") 
hom_nest_2 <- homtestmap(loc_nest_2, what = "pvalue") 
hom_nest_3 <- homtestmap(loc_nest_3, what = "pvalue") 

#Figure6 
par(mfrow = c(2, 2))
#a
plot(hom_nest_1, col = g025, main = "", ribside = "bottom", cex.axis = 1.5)  
contour(hom_nest_1, add = TRUE)
plot(ppp_gr, add = TRUE, pch = "+")
#b
plot(hom_nest_2, col = g025, main = "", ribside = "bottom", cex.axis = 1.5)
contour(hom_nest_2, add = TRUE)
plot(ppp_gr, add = TRUE, pch = "+")
#c
plot(hom_nest_3, col = g025, main = "", ribside = "bottom", cex.axis = 1.5)
contour(hom_nest_3, add = TRUE)
plot(ppp_gr, add = TRUE, pch = "+")
#d
plot(loc_gr, col = g025, main = "", ribside = "bottom", cex.axis = 1.5)
contour(loc_gr, add = TRUE)
plot(ppp_gr, add = TRUE, pch = "+")




#Fit global model
glo_gr <- ppm(ppp_gr ~ d_v + d_cs + d_pl + s(x, y, bs = "tp", k = 30), use.gam = TRUE, sigma = bw.diggle(ppp_gr))      

#compute smoothed raw residuals of both the global and local model
lambda <- density(ppp_gr, sigma = bw.diggle(ppp_gr)) 
#global
pr.glo_gr <- predict(glo_gr)            
lambda.glo_gr <- as.im(pr.glo_gr, W = as.owin(lambda)) 
rr.glo_gr <- eval.im(lambda - lambda.glo_gr)
#local
pr.loc_gr <- predict(loc_gr)            
lambda.loc_gr <- as.im(pr.loc_gr, W = as.owin(lambda)) 
rr.loc_gr <- eval.im(lambda - lambda.loc_gr)

#Figure7
par(mfrow = c(1, 2))
#a
plot(Smooth(rr.loc_gr), col = col_div, main = "", ribside = "bottom", cex.axis = 1.5) 
contour(Smooth(rr.loc_gr), add = TRUE)
#b
plot(Smooth(rr.glo_gr), col = col_div, main = "", ribside = "bottom", cex.axis = 1.5) 
contour(Smooth(rr.glo_gr), add = TRUE)


#compute envelopes of the inhomogeneous K-functions
#global
sim.glo_gr <- simulate.ppm(glo_gr, nsim = 10000, seed = 1)
#local
sim.loc_gr <- simulate.profilepl(loc_gr, nsim = 10000, seed = 1)


#Figure8
par(mfrow = c(1, 2))
#a
plot(envelope(ppp_gr, Kinhom, simulate = sim.loc_gr, nsim = 100, correction = "best"), main = "", ylim = c(0, 10), legend = FALSE, cex.axis = 1.5, xlab = "", ylab = "")
mtext(c("r", "Kinhom(r)"), side = c(1, 2), line = 2.5, cex = 2)
#b
plot(envelope(ppp_gr, Kinhom, simulate = sim.glo_gr, nsim = 100, correction = "best"), main = "", ylim = c(0, 10), legend = FALSE, cex.axis = 1.5, xlab = "", ylab = "")
mtext(c("r", "Kinhom(r)"), side = c(1, 2), line = 2.5, cex = 2)





#################################################################################
# Section 6.2 Local Log-Gaussian Cox Process model for a Greek seismic sequence #
#################################################################################

#define point process for greek analysis of Ithaki, Kefalonia and Zakynthos 
ppp_is <- ppp_gr[w_is]

#compute the distance variables 
psp_cs_is <- psp_cs[w_is]
psp_pl_is <- psp_pl[w_is]
d_cs_is <- distfun(psp_cs_is)
d_pl_is <- distfun(psp_pl_is)

#Figure9
par(mfrow = c(1, 1))
map("world", xlim = c(20, 21), ylim = c(37.5, 38.5), fill = TRUE, col = "grey", border = NA)
mtext(c("Latitude", "Longitude"), side = c(1, 2), line = 2.5, cex = 1.5)
plot(ppp_is, add = TRUE, pch = "+", cex = 1.25)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)
plot(psp_cs_is, add = TRUE, col = "blue")
plot(psp_pl_is, add = TRUE, col = "red")

#select bandwidth for the local Palm likelihood
sigma_is_palm <- bw.loccit(ppp_is)

#Fit local model
loc_is <- loccit(ppp_is, ~ d_cs_is * d_pl_is, clusters = "LGCP", sigma = sigma_is_palm) 


#Figure10
par(mfrow = c(1, 3))
#a
plot(with(loc_is[["modelpar"]], sqrt(sigma2)), main = "", log = TRUE, ribside = "bottom", cex.axis = 1.5, col = g0)
plot(ppp_is, add = TRUE, pch = "+")
#b
plot(with(loc_is[["modelpar"]], alpha),  main = "", log = TRUE, ribside = "bottom", cex.axis = 1.5, col = g0)
plot(ppp_is, add = TRUE, pch = "+")
#c
plot(with(loc_is[["modelpar"]], sqrt(mu)), main="", log = TRUE, ribside = "bottom", cex.axis = 1.5, col = g0)
plot(ppp_is, add = TRUE, pch = "+")





##################################################
# Appendix A. Analysis of the Italian seismicity #
##################################################

#define point process for italian analysis of Abruzzo 
pun_ab <- ppp_it[w_ab]


#compute the envelopes of the homogeneous K-funcion 
E_ab <- envelope(pun_ab, Kest, nsim = 39, rank = 1) 


#Figure11
par(mfrow = c(1, 2))
#a
map("italy", xlim = c(13, 14.8), ylim = c(41.6, 42.9),  border = NA, fill = TRUE, col = "lightgrey")
map("italy", xlim = c(13, 14.8), ylim = c(41.6, 42.9), add = TRUE)
mtext(c("Latitude", "Longitude"), side = c(1, 2), line = 2.5, cex = 1.5)
plot(pun_ab, add = TRUE, pch = "+", col = "black", cex = 2)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)
#b
plot(E_ab, legend = F, main = "", cex.axis = 1.5, xlab = "", ylab = "") 
mtext(c("r", "K(r)"), side = c(1, 2), line = 2.5, cex = 2)




#select bandwidth for the local composite likelihood
sigma_ab <- bw.locppm(pun_ab)

#Fit local model
loc_ab <- locppm(pun_ab, ~ min_dist + fz, sigma = sigma_ab, vcalc = "full")

#compute local tests 
t_ns <- ttestmap(loc_ab, "min_dist") 
t_fz <- ttestmap(loc_ab, "fz")

#Figure12
par(mfrow = c(2, 2))
#a
plot(loc_ab, which = 3, main = "", ribside = "bottom", cex.axis = 2, col = g0) 
contour(loc_ab, which = 3, add = TRUE, labcex = 2)
plot(pun_ab, add = TRUE, pch = "+")
#b
plot(t_fz, main = "", ribside = "bottom", cex.axis = 2, col = g0)      
contour(t_fz, add = TRUE, levels = c( - 2, 2), labcex = 2)
plot(pun_ab, add = TRUE, pch = "+")
#c
plot(loc_ab, which = 2, main = "", ribside = "bottom", cex.axis = 2, col = g0) 
contour(loc_ab, which = 2, add = TRUE, labcex = 2)
plot(pun_ab, add = TRUE, pch = "+")
#d
plot(t_ns, main = "", ribside = "bottom", cex.axis = 2, col = g0)      
contour(t_ns, add = TRUE, levels = c( - 2, 2), labcex = 2)
plot(pun_ab, add = TRUE, pch = "+")




#Fit nested models
loc_1 <- locppm(pun_ab, ~ 1, sigma = sigma_ab, vcalc = "full") 
loc_fz <- locppm(pun_ab, ~ min_dist, sigma = sigma_ab, vcalc = "full")

#perform 'homogeneity' tests
hom_1 <- homtestmap(loc_1, what = "pvalue") 
hom_fz <- homtestmap(loc_fz, what = "pvalue") 
hom_ab <- homtestmap(loc_ab, what = "pvalue") 


#Figure13
par(mfrow = c(1, 3))
#a
plot(hom_1, col = g05, main = "", ribside = "bottom", cex.axis = 1.5)  
contour(hom_1, add = TRUE, labcex = 1.5)
plot(pun_ab, add = TRUE, pch = "+")
#b
plot(hom_fz, col = g05, main = "", ribside = "bottom", cex.axis = 1.5)
contour(hom_fz, add = TRUE, labcex = 1.5)
plot(pun_ab, add = TRUE, pch = "+")
#c
plot(hom_ab, col = g05, main = "", ribside = "bottom", cex.axis = 1.5)
contour(hom_ab, add = TRUE, labcex = 1.5)
plot(pun_ab, add = TRUE, pch = "+")


#select bandwidth for the local Palm likelihood
sigma_ab_palm <- bw.loccit(pun_ab)

#Fit local model
loc_ab <- loccit(pun_ab, ~ fz + min_dist, clusters = "LGCP", sigma = sigma_ab_palm) 

#Figure14
par(mfrow = c(1, 2))
#a
plot(with(loc_ab[["modelpar"]], sqrt(sigma2)), main = "", log = TRUE, ribside = "bottom", cex.axis = 2, col = g0)
plot(pun_ab, add = TRUE, pch = "+")
#b
plot(with(loc_ab[["modelpar"]], alpha), main = "", log = TRUE, ribside = "bottom", cex.axis = 2, col = g0)
plot(pun_ab, add = TRUE, pch = "+")


