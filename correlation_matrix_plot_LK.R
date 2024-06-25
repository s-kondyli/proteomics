library(corrplot)
library(ggplot2)


# Select the columns with ratios for the correlation plots(assuming columns ending with '_ratio' contain the ratios)->adjust each time
ratio_columns <- data[ , grep("_p.adj$", names(data))]

# Step 3: Ensure all selected columns are numeric
ratio_columns <- data.frame(lapply(ratio_columns, as.numeric))

# Step 4: Compute the correlation matrix
cor_matrix <- cor(ratio_columns, use="complete.obs")


# New row and column names
new_labels <- c("1k stimulated vs unstimulated", "2k stimulated vs unstimulated", "5k stimulated vs unstimulated")

# Assign new row and column names to the correlation matrix
rownames(cor_matrix) <- new_labels
colnames(cor_matrix) <- new_labels

# Step 5: Create the correlation plot
corrplot(cor_matrix, method = "square", addCoef.col = "black", number.cex = 0.7, tl.col = "black", #tl.col for the color of the labels
         main = "Correlation matrix of the padj values") 


