#Scatter plot

library(ggplot2)

# Convert columns of interest to numeric

data$X5k.stimulated_vs_X5k.untreated_ratio <- as.numeric(data$X5k.stimulated_vs_X5k.untreated_ratio)

data$X2k.stimulated_vs_X2k.untreated_ratio <- as.numeric(data$X2k.stimulated_vs_X2k.untreated_ratio)


#create a df to plot the scatter plot with the columns you need -> adjust every time

df_scatter = data.frame(gene = data$name, x_axis = data$X2k.stimulated_vs_X2k.untreated_ratio, y_axis = data$X5k.stimulated_vs_X5k.untreated_ratio)

#This line calculates the correlation coefficients (r) between x and y values
r = round(cor(data$X2k.stimulated_vs_X2k.untreated_ratio, data$X5k.stimulated_vs_X5k.untreated_ratio), digits = 2)

#Add new columns to the dataframe based on a logical vector in the initial data
df_scatter <- df_scatter %>%
  mutate(Significant_in   = case_when(
    data$X5k.stimulated_vs_X5k.untreated_significant == T & data$X2k.stimulated_vs_X2k.untreated_significant == T  ~ "Both",
    data$X5k.stimulated_vs_X5k.untreated_significant == T ~ "significant in the 5k",
    data$X2k.stimulated_vs_X2k.untreated_significant == T ~ "significant in the 2k",
    data$X5k.stimulated_vs_X5k.untreated_significant == F & data$X2k.stimulated_vs_X2k.untreated_significant == F ~ "None"
  ))


#plot it 

ggplot(df_scatter, aes(x = x_axis, y = y_axis, color=Significant_in)) +
  geom_point() +
  labs(x = "2k T cells", y = "5k T cells", title = "Scatter plot of -log2 ratio of identified proteins") +
  scale_color_brewer(palette = "Spectral")  # Change to any ColorBrewer palette name you prefer
