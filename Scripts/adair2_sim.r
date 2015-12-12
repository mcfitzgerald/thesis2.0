# Simulates binding isotherm with Adair equation

# Define binding constants

k1 <- 1 

k2 <- 0

# Generate free ligand concentrations

L <- c(1:10 %o% 10^(-3:3))  # Outer product of (1,2..10) and (0.001,0.01..1000)

# Generate bound fraction values

B <- ((k1*L)+(2*(L^2)*k1*k2))/((1 + (k1*L) + (k1*k2*(L^2))))

# Plot results

plot(L,B)
