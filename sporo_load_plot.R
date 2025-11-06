
library(foreach)
library(doParallel)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)

sporo_load <- matrix(0,1000*1000,3)
idx <- c(18,20,40)*24
for (j in 1:1000) {
  result <- foreach(i=1:3,.packages=c("rootSolve")) %dopar% {
    sporozoite_record <- sporozoite_count(1000,round(game_male_record[j,idx[i]],0), round(game_female_record[j,idx[i]],0),3.44,0.8, 0.029,13.6, 2.7, 17.8, 2.2,4427)
    sporozoite_record
  }
  sporo_load[((j-1)*1000+1):(j*1000),1] <-ifelse(result[[1]][,4]>15,result[[1]][,7],0)
  sporo_load[((j-1)*1000+1):(j*1000),2] <-ifelse(result[[2]][,4]>15,result[[2]][,7],0)
  sporo_load[((j-1)*1000+1):(j*1000),3] <-ifelse(result[[3]][,4]>15,result[[3]][,7],0)
  print(j)
}


write.csv(sporo_load,"Sporo_load_10_7.csv",row.names = FALSE)



Sporo_count_10_5 <- read.csv("Sporo_load_10_5.csv")
Sporo_count_10_6 <- read.csv("Sporo_load_10_6.csv")
Sporo_count_10_7 <- read.csv("Sporo_load_10_7.csv")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
df <- data.frame(
  value = c(Sporo_count_10_5[,1],
            Sporo_count_10_6[,1],
            Sporo_count_10_7[,1]),
  group = factor(rep(c("10_5","10_6","10_7"),
                     c(nrow(Sporo_count_10_5),
                       nrow(Sporo_count_10_6),
                       nrow(Sporo_count_10_7))))
)

xb <- seq(-10000, 80000, by = 10000)
lab_intervals <- paste0("(", head(xb, -1), ",", tail(xb, -1), "]")
lab_intervals_math <- paste0("(", head(xb, -1)/10000, ",", tail(xb, -1)/10000, "]")
lab_intervals_math[1] <- "Uninfected\nor dead" 
lab_intervals_math[9] <- "(7,+∞)" 

df_binned <- df %>%
  mutate(bin = cut(value,
                   breaks = xb,
                   include.lowest = TRUE,
                   right = TRUE,
                   labels = lab_intervals_math)) %>%
  count(group, bin) %>%
  complete(group, bin, fill = list(n = 0))

# --- Plot ---
ggplot(df_binned, aes(x = bin, y = n, fill = group)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  scale_x_discrete(drop = FALSE, limits = lab_intervals_math) +
  scale_y_continuous(
    breaks = c(0, 500000, 1000000),   # set the ticks you want
    labels = c(
      expression(0),
      expression(5   %*% 10^5),
      expression(1 %*% 10^6)
    )           
  ) +
  scale_fill_manual(
    values = c("10_5" = "#E5E5E5",
               "10_6" = "#A9A9A9",
               "10_7" = "#404040"),
    labels = c(expression(10^5), expression(10^6), expression(10^7))
  ) +
  labs(title = "40 Days post infection",  
       x = expression(Number~of~Sporozoites~"(×10"^4*")"),
       y = "Frequency",
       fill = "Parasitemia") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(angle = 90, size = 20, vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.ticks=element_line(),
    panel.grid = element_blank(),
    panel.border = element_blank(),                    
    axis.line.x  = element_line(color = "black", size = 0.5), 
    axis.line.y  = element_line(color = "black", size = 0.5), 
    legend.position = c(0.95, 0.95),           
    legend.justification = c("right", "top"), 
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 10, l = 10),
    #plot.title = element_text(hjust = 0.5, size = 15)
  )


library(dplyr)
library(tidyr)
library(ggplot2)

lab_intervals_math_gap <- append(lab_intervals_math, "gap", after = 1)

df_binned_gap <- df_binned %>%
  mutate(bin = factor(bin, levels = lab_intervals_math_gap)) %>%
  complete(group, bin, fill = list(n = 0))

ggplot(df_binned_gap, aes(x = bin, y = n, fill = group)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  scale_x_discrete(
    drop = FALSE,
    limits = lab_intervals_math_gap,
    labels = ifelse(lab_intervals_math_gap == "gap", "", lab_intervals_math_gap)
  ) +
  scale_y_continuous(
    breaks = c(0, 500000, 1000000),
    labels = c(
      expression(0),
      expression(5 %*% 10^5),
      expression(1 %*% 10^6)
    )
  ) +
  scale_fill_manual(
    values = c("10_5" = "#E5E5E5",
               "10_6" = "#A9A9A9",
               "10_7" = "#404040"),
    labels = c(expression(10^5), expression(10^6), expression(10^7))
  ) +
  labs(
    title = "18 Days post infection",
    x = expression(Number~of~Sporozoites~"(×10"^4*")"),
    y = "Frequency",
    fill = "Parasitemia"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(angle = 90, size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.ticks=element_line(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x  = element_line(color = "black", size = 0.5),
    axis.line.y  = element_line(color = "black", size = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 10, l = 10)
  )

