


ggplot(df_scatter, aes(x = DIANN, y = PD, color= Sig) ) +
  geom_point( size=1.5) +
  theme_minimal() +
  geom_smooth(method=lm, color="black", se =F) +
  annotate(geom="text", x= min(df$DIANN) * 0.9, y=max(df$PD) * 0.9, label=paste0('r = ', r),
           color= "black", size = 4) +
  labs(title = "-log2 fold-change correlation plot of mutual proteins in DIANN 1.8.1 and PD", x = "DIANN", y = "PD") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )

df_scatter <- df %>%
  mutate(Sig = case_when(
    PD_significant == T & DIANN_significant == T  ~ "Both",
    PD_significant== T ~ "PD",
    DIANN_significant == T ~ "DIANN",
    PD_significant == F & DIANN_significant == F ~ "None"
  ))

### Ina's code

ggplot(df_scatter, aes(x = DIANN, y = PD, color= Sig) ) +
  geom_point( size=1.5) +
  theme_minimal() +
  geom_smooth(method=lm, color="black", se =F) +
  annotate(geom="text", x= min(df$DIANN) * 0.9, y=max(df$PD) * 0.9, label=paste0('r = ', r),
           color= "black", size = 4) +
  labs(title = "-log2 fold-change correlation plot of mutual proteins in DIANN 1.8.1 and PD", x = "DIANN", y = "PD") +
  theme_minimal() +
  theme (
    plot.title = element_text(size = 10),  # Adjust title size here
    plot.subtitle = element_text(size = 10),  # Adjust subtitle size here
    plot.caption = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0))  # Adjust margin here
  )


df_scatter <- df %>%
  mutate(Sig = case_when(
    PD_significant == T & DIANN_significant == T  ~ "Both",
    PD_significant== T ~ "PD",
    DIANN_significant == T ~ "DIANN",
    PD_significant == F & DIANN_significant == F ~ "None"
  )) %>% 
  filter(Sig != "None")
