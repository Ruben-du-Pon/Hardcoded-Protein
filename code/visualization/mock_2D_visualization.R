# Define the protein string
protein.string <- "HHPHHHPHPHHHPH"

# Determine the sequence length
sequence_length <- nchar(protein.string)

# Set the starting point
start_position <- sequence_length %/% 2

# Create a vector of colors based on the protein string
colors <- ifelse(strsplit(protein.string, '')[[1]] == 'H', 'red', 'blue')

# Create a function to simulate the folding
fold_protein <- function(start, end, level) {
  lines(c(start, end), c(level, level), col = 'black', lwd = 2)
}

# Plot the protein folding
plot(1:sequence_length, rep(1, sequence_length), type = "n", xlab = "Position", ylab = "",
     main = "Protein Folding", xlim = c(1, sequence_length*2), ylim = c(0, 2))

# Add points with different colors
points(start_position:(start_position + sequence_length - 1), rep(1, sequence_length),
       pch = 16, col = colors)

# Simulate the folding by connecting points with lines
for (i in 1:(sequence_length - 1)) {
  fold_protein(start_position + i - 1, start_position + i, 1)
}

legend("topright", inset = c(0.8, 0.1),                    
       legend = c("Hydrofoob","Polair"), 
       pch = 16, col = c('red', 'blue'))
