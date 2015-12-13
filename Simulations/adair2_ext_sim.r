# Simulates binding isotherm with Adair equation

# Define binding constants

k1 <- as.numeric(readline(prompt="Enter k1: ")) 

k2 <- as.numeric(readline(prompt="Enter k2: "))

# Generate free ligand concentrations

L <- c((seq(1,10,by=0.1)) %o% 10^(-5:5))  # Outer product of (1,2..10) 
                                          # and (0.001,0.01..1000)

# Generate bound fraction values

B <- (1/2)*(((k1*L)+(2*(L^2)*k1*k2))/(((1 + (k1*L) + (k1*k2*(L^2))))))

# Plot results

plot(L,B, log="x")

sat90_ind <- min(which(abs(B-0.9)==min(abs(B-0.9))))

print(sat90_ind)

sat10_ind <- min(which(abs(B-0.1)==min(abs(B-0.1))))

print(sat10_ind)

bind_int <- (log10(L[sat90_ind]) - log10(L[sat10_ind]))

print(sprintf("The binding interval is %f log units", bind_int))
