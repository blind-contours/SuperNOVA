library(data.table)
# Test Script

# Create a data frame
test_df <- data.frame(
  noshift = c(0.1, 0.5, 0.9),
  downshift = c(0.2, 0.4, 0.8),
  upshift = c(0.3, 0.6, 0.7),
  upupshift = c(0.4, 0.7, 0.6)
)

# Run the function
test_est_hn <- est_hn(test_df)

# Compare the output
expected_output <- data.table(
  noshift = c(0.6, 0.8, 0.8888889),
  shift = c(1.1111111, 0.8333333, 1.2857143)
)

# Check if the outputs are equivalent
stopifnot(all.equal(test_est_hn, expected_output, tolerance = 0.1))

# Print a success message
print("Test successful!")
