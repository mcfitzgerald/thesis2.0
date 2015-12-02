# R script works on two column csv files containing free ligand concentration
# (col header "L") and bound fraction (col header "B" w/ values 0..1). The 
# variables "name.trunk" and "data.file" determine the output filenames 
# and data source respectively. This script can be run as part of a shell
# script that loops over all *.csv in a directory. It will generate 5 files:
# A pdf of a semi-log plot of the binding data, a pdf of a direct plot with 
# two fitted lines for nls solutions to the Adair equation for 1 and 2 classes
# of binding sides, one pdf each for the residuals of the fits, and one .Rout
# file that contains the fit parameters, and and ANOVA comparison of the fits

# EDIT EACH RUN - Data source and file naming

name.trunk <- "yard_no-linker"

data.file <-"yard_no-linker.csv"

# Output messages

sink(paste(name.trunk,".Rout",sep=""))

# Load libraries

library(minpack.lm)  # For better nls fitting

#Plot options
graph.xax <- "Free EGF [nM]"

graph.yax <- "Bound Fraction"

#Open data file

data.loaded <- read.csv(data.file)

#Generate SemiLog Plot

pdf(file=paste(name.trunk, "SemiLog.pdf"), height=5, width=5) 

par(cex.main=0.9)  #Scale title font

plot(data.loaded$L, data.loaded$B, log="x", xlab=graph.xax, ylab=graph.yax, 
main=name.trunk)

abline(h=c(0.1,0.9),col="red")

dev.off()

# Build Models (remember this assumes L and B headers)

adair1 <- nlsLM(B ~ (k1*L)/(1 + (k1*L)), data=data.loaded, start=list(k1=1), 
trace=T)

summary(adair1)

adair2 <- nlsLM(B ~ (((k1*L)+(2*(L^2)*k2*k1))/(2*(1 + (k1*L) + (k1*k2*(L^2))))),
data=data.loaded, start=list(k1=1,k2=0.5), trace=T)

summary(adair2)

# Compare models 
anova(adair1,adair2)

AIC(adair1,adair2)

# Build predictor L values

pred.lvals <- seq((min(data.loaded$L)), (max(data.loaded$L)), length.out=100)

# Model prediction

adair1_fit <- predict(adair1, newdata=list(L=pred.lvals))

adair2_fit <- predict(adair2,newdata=list(L=pred.lvals))

# Plot fits
pdf(file=paste(name.trunk, "DirectFit.pdf"), height=5, width=5) 

par(cex.main=0.9)  #Scale title font

plot(data.loaded$L, data.loaded$B, xlab=graph.xax, ylab=graph.yax, 
main=name.trunk)

lines(pred.lvals, adair1_fit, col="blue")

lines(pred.lvals, adair2_fit, col="red", lty=2)

legend("bottomright", c("raw data","adair1","adair2"), cex=0.8, col=c("black",
"blue", "red"), pch=c("o","-","-"))

dev.off()

# Residual Plots

adair1_resid <- resid(adair1)
adair2_resid <- resid(adair2)

pdf(file=paste(name.trunk, "Residual1.pdf"), height=5, width=5) 

par(cex.main=0.9)  #Scale title font

plot(data.loaded$L, adair1_resid, log="x", xlab=graph.xax, ylab="residuals",
main=paste(name.trunk, "Adair1 Residual"))

dev.off()

pdf(file=paste(name.trunk, "Residual2.pdf"), height=5, width=5) 

par(cex.main=0.9)  #Scale title font

plot(data.loaded$L, adair2_resid, log="x",  xlab=graph.xax, ylab="residuals", 
main=paste(name.trunk, "Adair2 Residual"))

dev.off()
