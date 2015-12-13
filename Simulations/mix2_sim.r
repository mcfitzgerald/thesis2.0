# Simulates binding of system containing a distribution of monomers
# and dimers with a possible total of three distinguishable binding sites
# the monomer is described by the adair equation for one site, and the dimer
# dimer is described by a klotz equation for two sites allowing if existence
# of cooperativity

# Define binding constants

k1 <- as.numeric(readline(prompt="Enter k1 (k for monomer): "))
print(k1)

k2 <- as.numeric(readline(prompt="Enter k2 (first k of dimer): ")) 
print(k2)

k3 <- as.numeric(readline(prompt="Enter k3 (second k of dimer): "))
print(k3)

m1 <- as.numeric(readline(prompt="Enter mix factor 1 (0..1): "))
print(m1)


m2 <- (1.0 - m1)

print(m2)
# Generate free ligand concentrations

# Outer product of (1,2..10) and (0.001,0.01..1000)
L <- c((seq(1,10,by=0.1)) %o% 10^(-5:5)) 

# Generate bound fraction values

B <- ((m1*((k1*L)/(1 + (k1*L)))) + 
      (m2*(((k2*L)+(2*(L^2)*k2*k3))/(2*(((1 + (k2*L) + (k2*k3*(L^2)))))))))

# Plot results

plot(L,B,log="x")

# Compute binding interval

sat90_ind <- min(which(abs(B-0.9)==min(abs(B-0.9))))
print(sat90_ind)
sat10_ind <- min(which(abs(B-0.1)==min(abs(B-0.1))))
print(sat10_ind)
bind_int <- (log10(L[sat90_ind]) - log10(L[sat10_ind]))

print(sprintf("The binding interval is %f log units", bind_int))
