
#Step 1: Reading the CSV files

a <- read.csv("~~~/transcritp_poly(T)-length_raw_data.csv", header = T, row.name = NULL, check.names=F, fileEncoding = "UTF-8")

#This line reads the CSV file transcritp_poly(T)-length_raw_data.csv using the read.csv() function. Adjust the file path to the actual location of your file. The header = T parameter specifies that the file has a header row.

#Step 2: Creating a matrix with diagonal frequency elements

matrix_nxn <- diag(a[2:134,4], nrow = 133, ncol = 133)
#This line creates a square matrix matrix_nxn with dimensions 133x133. It takes the column 4 data from rows 2 to 134 of a and uses it as the diagonal elements of the matrix. All other elements outside the diagonal are set to 0.

#Step 3: Printing the resulting matrix

print(matrix_nxn)

#This line prints the resulting matrix matrix_nxn.

#Step 4: Reading another CSV file for the weighting factor

wf <- read.csv("~~~/revision_data/t12-wf.csv", header = T, row.name = NULL, check.names=F, fileEncoding = "UTF-8")

#This line reads the CSV file t12-wf.csv using the read.csv() function. Adjust the file path to match the actual location of your file. The header = T parameter indicates that the file has a header row.

#Step 5: Creating a matrix for the weighting factor

matrix_nx1 <- matrix(wf[4:136,4], nrow = 133, ncol = 1)
matrix_nx1 <- as.numeric(matrix_nx1)

#These lines create a matrix matrix_nx1 with dimensions 133x1. It takes the column 4 data from rows 4 to 136 of wf and assigns it to the matrix. The as.numeric() function is used to convert the elements of the matrix to numeric type.

#Step 6: Performing matrix multiplication

product <- matrix_nxn %*% matrix_nx1

#This line performs matrix multiplication between matrix_nxn and matrix_nx1 using the %*% operator. The resulting matrix is stored in the variable product.

#Step 7: Writing the result of weighted relative frequency corresponding to each length size to a file

write.table(product, '~~~', sep = '\t', col.names = NA, quote = FALSE)

#This line writes the matrix product to a file. Replace '~~~' with the desired file path. The resulting output will be tab-separated with no column names and without quotes.



