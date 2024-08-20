# Given parameters
E_S <- 750  # Target species richness
J <- 250000     # Total number of individuals (example value)
tau <- 10    # Interaction parameter (example value)

# Define the function to find the root
f <- function(mu) {
  if (mu <= 0 || mu >= 1) {
    return(Inf)  # Avoid division by zero or negative values
  }
  return((J * mu / (1 - mu)) * log((1 + tau * mu) / (mu + tau * mu)) - E_S)
}

# Initial guess for the interval where the root lies
lower_bound <- 0.0000000000001
upper_bound <- 0.999

# Solve for mu using uniroot
solution <- uniroot(f, lower = lower_bound, upper = upper_bound)

# Estimated speciation rate (mu)
mu_solution <- solution$root

print(paste("Estimated speciation rate (mu):", mu_solution))