# Input series that is obtained from R2 reads against each transcript
polyT_series <- seq(10, 143)

# Calculate the output series using the linear equation
weighting_factor <- (input_series + 11.4) / 0.49

# Print the output series
print(weighting_factor)

#normalizing the weighting factor
sum_wf <- sum(weighting_factor)

norm_wf <- weighting_factor/sum_wf

# Writing the normalized weighting factor to a file
write.table(norm_wf, '~~~', sep = '\t', col.names = NA, quote = FALSE)