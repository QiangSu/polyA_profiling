# Input series obtained from R2 reads against each transcript
polyT_series <- seq(10, 143)

# Calculate the output series using the linear equation as scalar product
polyA_series <- (polyT_series + 6.13) / 0.65

# Print the output series
print(polyA_series)

# Writing the polyA length to a file
write.table(norm_wf, '~~~', sep = '\t', col.names = NA, quote = FALSE)