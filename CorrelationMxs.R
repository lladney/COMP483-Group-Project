# Load in and Install Library
install.packages("corrplot")
library(corrplot)

# Load the csv file using the read.csv function.
data <- read.csv("protein.csv")

# Create a correlation matrix using the cor function.
corr_matrix <- cor(data[,2:21])

# Transpose the matrix 
tranposedmatrix = cor(t(data[,2:21]))
colnames(tranposedmatrix)=data[,1]
rownames(tranposedmatrix) = data[,1]

# Create the heatmap using the corrplot function.
corrplot(tranposedmatrix, type = "upper", method = "color", 
                  tl.col = "black", tl.srt = 45)

# This will create a boxplot for each amino acid, showing the distribution of frequencies for each one. 
# You can compare the boxplots to see if there are any significant differences in the frequency distributions between the amino acids.
boxplot(data[,2], data[,3], data[,4], data[,5], data[,6], data[,7], data[,8], data[,9], data[,10], data[,11], 
        data[,12], data[,13], data[,14], data[,15], data[,16], data[,17], data[,18], data[,19], data[,20], data[,21],
        col = c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"))

'''
# Transpose the matrix 
tranposedmatrix = t(corr_matrix)

# Create the heatmap using the corrplot function.
corrplot(corr_matrix, type = "upper", method = "color", 
         tl.col = "black", tl.srt = 45)
'''
