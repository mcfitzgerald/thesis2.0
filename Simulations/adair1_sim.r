# Simulates binding isotherm with Adair equation

# Define binding constants

k1 <- 5.3e+9

# Generate free ligand concentrations

# Outer product of (1,2..10) and (0.001,0.01..1000)
L <- c((seq(1,10,by=0.1)) %o% 10^(-6:6))

# Generate bound fraction values

B <- (k1*L)/(1 + (k1*L))

# Plot results

plot(L,B,log="x")
