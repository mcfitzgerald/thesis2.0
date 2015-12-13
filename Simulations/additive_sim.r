# Simulates binding isotherm with Adair equation

# Define binding constants

k1 <- 0.1

k2 <- 10 

# Generate free ligand concentrations

L <- c(1:10 %o% 10^(-3:3))  # Outer product of (1,2..10) and (0.001,0.01..1000)

# Generate bound fraction values

B <- (1/2)*(((k1*L)/(1 + (k1*L))) + ((k2*L)/(1 + (k2*L))))

# Plot results

plot(L,B,log="x")

# Compute binding intervale

sat90_ind <- min(which(abs(B-0.9)==min(abs(B-0.9))))

sat10_ind <- min(which(abs(B-0.1)==min(abs(B-0.1))))

bind_int <- (log10(L[sat90_ind]) - log10(L[sat10_ind]))

print(sprintf("The binding interval is %f log units", bind_int))

print("Thanks!")
