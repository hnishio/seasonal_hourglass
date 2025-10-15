
# Load packages
library(tidyverse)
library(patchwork)
library(data.table)
library(viridis)

# Create output directory
out <- "03_HAI_seasonal/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}



##### H3K27me3 #####

### Plotting (daily temp)
df_hai <- as.data.frame(fread(paste0("03_K27AI_250427/", "K27ai_rawdata_maxRNAover2_repmean.csv")))
attr(df_hai$date, "tzone") <- "Asia/Tokyo"
obs_var <- var(df_hai$hai)

df_permuted_HAI <- as.data.frame(fread(file=paste0("03_K27AI_250427/", "permutedK27ai_rawdata_maxRNAover2_repmean.csv"), header = T))
permuted_var <- apply(df_permuted_HAI, 1, var)
p_val <- sum(permuted_var >= obs_var) / length(permuted_var)
if(p_val == 0){
  p_val <- 1 / length(permuted_var)
}

# Load Nishiwaki dailytemp
df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
df_temp2 <- df_temp %>%
  filter(date >= df_hai$date[1] & date <= df_hai$date[nrow(df_hai)])

# Setting var_nameiables
labelx<-c("","11","12",
          "1","2","3","4","5","6","7","8","9","10","")
length(labelx)

local({
  xmax <- max(df_hai[,6], na.rm = T)
  xmin <- min(df_hai[,6], na.rm = T)
  ymax <- 10.9
  ymin <- min(df_hai[,2:3], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_hai <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_hai_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
  df_temp2$temp_hai <- scale_to_hai(df_temp2$dailytemp)
  
  g_ts_K27 <<- ggplot() +
    geom_line(data = df_temp2, aes(x = date, y = temp_hai), col = "gray50", linewidth = 0.2) + 
    geom_ribbon(data = df_hai, aes(x = date, ymin = q2.5, ymax = q97.5),
                alpha = 0.3, fill = "gray30") +
    geom_line(data = df_hai, aes(x = date, y = hai), col = "black", linewidth = 0.6) +
    # geom_point(data = df_hai, aes(x = date, y = hai), col = "black", alpha = 0.6, size = 0.2) +
    # scale_x_continuous(breaks = 1:12) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    annotate("text", 
             x=as.POSIXct("2012-12-31"), y=yceiling-yrange*0.13,
             label = bquote(italic(p) < .(format(p_val, scientific = FALSE))), 
             size = 6/ggplot2::.pt) +
    theme_test(base_size=6) +
    theme(plot.title = element_text(size = 6, face = "bold"),
          axis.title=element_text(size=6, colour = "black"), 
          axis.text=element_text(size=6, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(0,0.2,0.5,0.2), "mm"),
          plot.tag = element_text(size = 8, face = "bold")) +
    labs(title = "H3K27me3", x = "", y = "Age index") +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       breaks = seq(0, 40, by = 0.5),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_hai_rev(.), name = "Temperature (°C)")) +
    annotate(geom = "text", label = "Young", 
             x=as.POSIXct("2012-09-20"), y=yceiling+yrange*0.05,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "Old", 
             x=as.POSIXct("2012-09-20"), y=yfloor-yrange*0.05,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2012",
             x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.41,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013",
             x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.41,
             size = 6/ggplot2::.pt, vjust=0) +
    geom_segment(data = df_hai, mapping = aes(x = date, y = hai),
                 x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.24,
                 xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.24, linewidth = 0.2) +
    geom_segment(data = df_hai, mapping = aes(x = date, y = hai),
                 x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.24,
                 xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.24, linewidth = 0.2)
})





##### Each PS mean accumulation #####

### Comparison between PSs
df_merge2 <- as.data.frame(fread(file=paste0("03_K27AI_250427/", "PS_log2rpkm_maxRNAover2_K27repmean.csv")))
df_merge2[1:5,1:5]

# Calculate mean expression for each PS at each date
df_eachPS <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), mean, .names = "{.col}"))

dates <- c("2012-11-6 12:00 JST", "2012-12-4 12:00 JST",
           "2013-1-8 12:00 JST", "2013-2-5 12:00 JST",
           "2013-3-5 12:00 JST", "2013-4-2 12:00 JST",
           "2013-4-30 12:00 JST", "2013-5-28 12:00 JST",
           "2013-7-2 12:00 JST", "2013-7-30 12:00 JST",
           "2013-8-27 12:00 JST", "2013-9-24 12:00 JST")

names(df_eachPS) <- c("PSnum", dates)

# Calculate SE expression for each PS at each date
se <- function(x){sd(x) / sqrt(length(x))}
df_eachPS_se <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), se, .names = "{.col}"))
df_eachPS_minusse <- cbind(df_eachPS_se[,1], (df_eachPS[,2:ncol(df_eachPS)] - df_eachPS_se[,2:ncol(df_eachPS_se)]))
df_eachPS_plusse <- cbind(df_eachPS_se[,1], (df_eachPS[,2:ncol(df_eachPS)] + df_eachPS_se[,2:ncol(df_eachPS_se)]))
names(df_eachPS_minusse) <- c("PSnum", dates)
names(df_eachPS_plusse) <- c("PSnum", dates)

df_plot <- df_eachPS %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean"
  ) %>%
  mutate(date = as.POSIXct(date))  

df_plot_minusse <- df_eachPS_minusse %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean_minus_se"
  ) %>%
  mutate(date = as.POSIXct(date))

df_plot_plusse <- df_eachPS_plusse %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean_plus_se"
  ) %>%
  mutate(date = as.POSIXct(date))

df_plot <- inner_join(df_plot, df_plot_minusse, by = c("PSnum", "date"))
df_plot <- inner_join(df_plot, df_plot_plusse, by = c("PSnum", "date"))


### PS1 and PS18 with temperature
# Load Nishiwaki dailytemp
df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
df_temp2 <- df_temp %>%
  filter(date >= df_plot$date[1] & date <= df_plot$date[nrow(df_plot)])

dates <- as.POSIXct(names(df_eachPS)[-1])
temp_sample <- df_temp2$dailytemp[df_temp2$date %in% dates]

# Plotting
glist <- list()

# PS1
local({
  i <- 1
  df_plot_rev <- df_plot %>% filter(PSnum == i)
  xmax <- max(df_temp2$date, na.rm = T)
  xmin <- min(df_temp2$date, na.rm = T)
  ymax <- max(df_plot_rev[,3:5], na.rm = T)
  ymin <- min(df_plot_rev[,3:5], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  yceiling <-  yceiling + yrange * 0.1
  
  # Scale adjustment
  scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
  df_temp2$temp_rna <- scale_to_rna(df_temp2$dailytemp)
  
  # Setting labels
  labelx<-c("","11","12",
            "1","2","3","4","5","6","7","8","9","10","")
  length(labelx)
  PSname <- df_merge2$PSname[min(which(df_merge2$PSnum==i))]
  
  # correlation
  cor_res <- cor.test(df_plot_rev$mean, temp_sample, method = "pearson")
  mark <- ifelse(cor_res$p.value < 0.0001, "****", 
                 ifelse(cor_res$p.value < 0.001, "***", 
                        ifelse(cor_res$p.value < 0.01, "**",
                               ifelse(cor_res$p.value < 0.05, "*", "ns"))))
  r <- format(cor_res$estimate, digits=2, nsmall = 2)
  rlab <- bquote(paste(italic(r), " = ", .(r), .(mark), sep=""))
  
  glist[[i]] <<- ggplot() +
    geom_line(data = df_temp2, aes(x = date, y = temp_rna), col = "gray50", linewidth = 0.2) + 
    geom_ribbon(data = df_plot_rev, aes(x = date, ymin = mean_minus_se, ymax = mean_plus_se),
                alpha = 0.3, fill = "gray30") +
    geom_line(data = df_plot_rev, aes(x = date, y = mean), color = "black", linewidth = 0.5) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_rna_rev(.), name = "")) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    theme_test(base_size = 6) +
    theme(plot.title = element_text(size = 6, face = "bold"),
          axis.title=element_text(size = 6, colour = "black"), 
          axis.title.y.right = element_blank(), 
          axis.text=element_text(size = 6, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "",
         y = "Mean H3K27me3") +
    annotate("text", label=rlab,
             x=as.POSIXct("2013-08-15"), y=yceiling-yrange*0.08, 
             size = 6/ggplot2::.pt, hjust = 0) +
    annotate(geom = "text", label = "2012", 
             x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013", 
             x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    geom_segment(x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.2, linewidth = 0.2)
})

# PS18
local({
  i <- 18
  df_plot_rev <- df_plot %>% filter(PSnum == i)
  xmax <- max(df_temp2$date, na.rm = T)
  xmin <- min(df_temp2$date, na.rm = T)
  ymax <- max(df_plot_rev[,3:5], na.rm = T)
  ymin <- min(df_plot_rev[,3:5], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  yceiling <-  yceiling + yrange * 0.1
  
  # Scale adjustment
  scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
  df_temp2$temp_rna <- scale_to_rna(df_temp2$dailytemp)
  
  # Setting labels
  labelx<-c("","11","12",
            "1","2","3","4","5","6","7","8","9","10","")
  length(labelx)
  PSname <- df_merge2$PSname[min(which(df_merge2$PSnum==i))]
  
  # correlation
  cor_res <- cor.test(df_plot_rev$mean, temp_sample, method = "pearson")
  mark <- ifelse(cor_res$p.value < 0.0001, "****", 
                 ifelse(cor_res$p.value < 0.001, "***", 
                        ifelse(cor_res$p.value < 0.01, "**",
                               ifelse(cor_res$p.value < 0.05, "*", "ns"))))
  r <- format(cor_res$estimate, digits=2, nsmall = 2)
  rlab <- bquote(paste(italic(r), " = ", .(r), .(mark), sep=""))
  
  glist[[i]] <<- ggplot() +
    geom_line(data = df_temp2, aes(x = date, y = temp_rna), col = "gray50", linewidth = 0.2) + 
    geom_ribbon(data = df_plot_rev, aes(x = date, ymin = mean_minus_se, ymax = mean_plus_se),
                alpha = 0.3, fill = "gray30") +
    geom_line(data = df_plot_rev, aes(x = date, y = mean), color = "black", linewidth = 0.5) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_rna_rev(.), name = "")) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    theme_test(base_size = 6) +
    theme(plot.title = element_text(size = 6, face = "bold"),
          axis.title=element_text(size = 6, colour = "black"), 
          axis.title.y.right = element_blank(), 
          axis.text=element_text(size = 6, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "",
         y = "Mean H3K27me3") +
    annotate("text", label=rlab,
             x=as.POSIXct("2013-08-15"), y=yceiling-yrange*0.08, 
             size = 6/ggplot2::.pt, hjust = 0) +
    annotate(geom = "text", label = "2012", 
             x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013", 
             x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    geom_segment(x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.2, linewidth = 0.2)
})



### All PS
xmax <- max(df_temp2$date, na.rm = T)
xmin <- min(df_temp2$date, na.rm = T)
ymax <- max(df_plot[,3], na.rm = T)
ymin <- min(df_plot[,3], na.rm = T)
yrange <- (ymax - ymin)
yceiling <-  ymax + yrange * 0.05
yfloor <- ymin - yrange * 0.05

# Scale adjustment
scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
df_temp2$temp_rna <- scale_to_rna(df_temp2$dailytemp)

# Setting labels
labelx<-c("","11","12",
          "1","2","3","4","5","6","7","8","9","10","")

g_K27 <- ggplot(df_plot, aes(x = date, y = mean, group = PSnum, color = PSnum)) +
  geom_line() +
  scale_color_viridis(name="PS", option="turbo") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                     expand = c(0, 0)) +
  scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
  coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
  theme_test(base_size = 6) +
  theme(plot.title = element_text(size = 6, face = "bold"),
        axis.title=element_text(size = 6, colour = "black"), 
        axis.text=element_text(size = 6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin=unit(c(0, 6.5, 0.5, 0.5),"mm"),
        plot.tag = element_text(size = 8, face = "bold"),
        legend.title=element_text(size=6, color = "black"),
        legend.text=element_text(size=6, color = "black"),
        legend.position=c(1.13, 0.5),
        legend.key.height=unit(3, "mm"),
        legend.key.width = unit(1.5, "mm"),
        legend.background = element_rect(fill="transparent"),
        legend.spacing.x = unit(0.7, 'mm')) +
  labs(title = "H3K27me3",
       x = "",
       y = "Mean level") +
  annotate(geom = "text", label = "2012",
           x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.41,
           size = 6/ggplot2::.pt, vjust=0) +
  annotate(geom = "text", label = "2013",
           x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.41,
           size = 6/ggplot2::.pt, vjust=0) +
  geom_segment(x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.24,
               xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.24, linewidth = 0.2,
               colour = "black") +
  geom_segment(x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.24,
               xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.24, linewidth = 0.2,
               colour = "black")






##### H3K4me3 #####

### Plotting (daily temp)
df_hai <- as.data.frame(fread(paste0("02_K4AI_250427/", "K4ai_rawdata_maxRNAover2_repmean.csv")))
attr(df_hai$date, "tzone") <- "Asia/Tokyo"
obs_var <- var(df_hai$hai)

df_permuted_HAI <- as.data.frame(fread(file=paste0("02_K4AI_250427/", "permutedK4ai_rawdata_maxRNAover2_repmean.csv"), header = T))
permuted_var <- apply(df_permuted_HAI, 1, var)
p_val <- sum(permuted_var >= obs_var) / length(permuted_var)
if(p_val == 0){
  p_val <- 1 / length(permuted_var)
}

# Load Nishiwaki dailytemp
df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
df_temp2 <- df_temp %>%
  filter(date >= df_hai$date[1] & date <= df_hai$date[nrow(df_hai)])

# Setting var_nameiables
labelx<-c("","11","12",
          "1","2","3","4","5","6","7","8","9","10","")
length(labelx)

local({
  xmax <- max(df_hai[,6], na.rm = T)
  xmin <- min(df_hai[,6], na.rm = T)
  ymax <- 9.26
  ymin <- 8.96
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_hai <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_hai_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
  df_temp2$temp_hai <- scale_to_hai(df_temp2$dailytemp)
  
  g_ts_K4 <<- ggplot() +
    geom_line(data = df_temp2, aes(x = date, y = temp_hai), col = "gray50", linewidth = 0.2) + 
    geom_ribbon(data = df_hai, aes(x = date, ymin = q2.5, ymax = q97.5),
                alpha = 0.3, fill = "gray30") +
    geom_line(data = df_hai, aes(x = date, y = hai), col = "black", linewidth = 0.6) +
    # geom_point(data = df_hai, aes(x = date, y = hai), col = "black", alpha = 0.6, size = 0.2) +
    # scale_x_continuous(breaks = 1:12) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    annotate("text", 
             x=as.POSIXct("2012-12-31"), y=yceiling-yrange*0.13,
             label = bquote(italic(p) < .(format(p_val, scientific = FALSE))), 
             size = 6/ggplot2::.pt) +
    theme_test(base_size=6) +
    theme(plot.title = element_text(size = 6, face = "bold"),
          axis.title=element_text(size=6, colour = "black"), 
          axis.text=element_text(size=6, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(0,0.2,0.5,0.2), "mm"),
          plot.tag = element_text(size = 8, face = "bold")) +
    labs(title = "H3K4me3", x = "", y = "Age index") +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       breaks = seq(0, 40, by = 0.1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_hai_rev(.), name = "Temperature (°C)")) +
    annotate(geom = "text", label = "Young", 
             x=as.POSIXct("2012-09-20"), y=yceiling+yrange*0.05,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "Old", 
             x=as.POSIXct("2012-09-20"), y=yfloor-yrange*0.05,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2012",
             x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.41,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013",
             x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.41,
             size = 6/ggplot2::.pt, vjust=0) +
    geom_segment(data = df_hai, mapping = aes(x = date, y = hai),
                 x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.24,
                 xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.24, linewidth = 0.2) +
    geom_segment(data = df_hai, mapping = aes(x = date, y = hai),
                 x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.24,
                 xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.24, linewidth = 0.2)
})





##### Each PS mean accumulation #####

### Comparison between PSs
df_merge2 <- as.data.frame(fread(file=paste0("02_K4AI_250427/", "PS_log2rpkm_maxRNAover2_K4repmean.csv")))
df_merge2[1:5,1:5]

# Calculate mean expression for each PS at each date
df_eachPS <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), mean, .names = "{.col}"))

dates <- c("2012-11-6 12:00 JST", "2012-12-4 12:00 JST",
           "2013-1-8 12:00 JST", "2013-2-5 12:00 JST",
           "2013-3-5 12:00 JST", "2013-4-2 12:00 JST",
           "2013-4-30 12:00 JST", "2013-5-28 12:00 JST",
           "2013-7-2 12:00 JST", "2013-7-30 12:00 JST",
           "2013-8-27 12:00 JST", "2013-9-24 12:00 JST")

names(df_eachPS) <- c("PSnum", dates)

# Calculate SE expression for each PS at each date
se <- function(x){sd(x) / sqrt(length(x))}
df_eachPS_se <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), se, .names = "{.col}"))
df_eachPS_minusse <- cbind(df_eachPS_se[,1], (df_eachPS[,2:ncol(df_eachPS)] - df_eachPS_se[,2:ncol(df_eachPS_se)]))
df_eachPS_plusse <- cbind(df_eachPS_se[,1], (df_eachPS[,2:ncol(df_eachPS)] + df_eachPS_se[,2:ncol(df_eachPS_se)]))
names(df_eachPS_minusse) <- c("PSnum", dates)
names(df_eachPS_plusse) <- c("PSnum", dates)

df_plot <- df_eachPS %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean"
  ) %>%
  mutate(date = as.POSIXct(date))  

df_plot_minusse <- df_eachPS_minusse %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean_minus_se"
  ) %>%
  mutate(date = as.POSIXct(date))

df_plot_plusse <- df_eachPS_plusse %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean_plus_se"
  ) %>%
  mutate(date = as.POSIXct(date))

df_plot <- inner_join(df_plot, df_plot_minusse, by = c("PSnum", "date"))
df_plot <- inner_join(df_plot, df_plot_plusse, by = c("PSnum", "date"))


### PS1 and PS18 with temperature
# Load Nishiwaki dailytemp
df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
df_temp2 <- df_temp %>%
  filter(date >= df_plot$date[1] & date <= df_plot$date[nrow(df_plot)])

dates <- as.POSIXct(names(df_eachPS)[-1])
temp_sample <- df_temp2$dailytemp[df_temp2$date %in% dates]

# Plotting
glist <- list()

# PS1
local({
  i <- 1
  df_plot_rev <- df_plot %>% filter(PSnum == i)
  xmax <- max(df_temp2$date, na.rm = T)
  xmin <- min(df_temp2$date, na.rm = T)
  ymax <- max(df_plot_rev[,3:5], na.rm = T)
  ymin <- min(df_plot_rev[,3:5], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  yceiling <-  yceiling + yrange * 0.1
  
  # Scale adjustment
  scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
  df_temp2$temp_rna <- scale_to_rna(df_temp2$dailytemp)
  
  # Setting labels
  labelx<-c("","11","12",
            "1","2","3","4","5","6","7","8","9","10","")
  length(labelx)
  PSname <- df_merge2$PSname[min(which(df_merge2$PSnum==i))]
  
  # correlation
  cor_res <- cor.test(df_plot_rev$mean, temp_sample, method = "pearson")
  mark <- ifelse(cor_res$p.value < 0.0001, "****", 
                 ifelse(cor_res$p.value < 0.001, "***", 
                        ifelse(cor_res$p.value < 0.01, "**",
                               ifelse(cor_res$p.value < 0.05, "*", "ns"))))
  r <- format(cor_res$estimate, digits=2, nsmall = 2)
  rlab <- bquote(paste(italic(r), " = ", .(r), .(mark), sep=""))
  
  glist[[i]] <<- ggplot() +
    geom_line(data = df_temp2, aes(x = date, y = temp_rna), col = "gray50", linewidth = 0.2) + 
    geom_ribbon(data = df_plot_rev, aes(x = date, ymin = mean_minus_se, ymax = mean_plus_se),
                alpha = 0.3, fill = "gray30") +
    geom_line(data = df_plot_rev, aes(x = date, y = mean), color = "black", linewidth = 0.5) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_rna_rev(.), name = "")) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    theme_test(base_size = 6) +
    theme(plot.title = element_text(size = 6, face = "bold"),
          axis.title=element_text(size = 6, colour = "black"), 
          axis.title.y.right = element_blank(), 
          axis.text=element_text(size = 6, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "",
         y = "Mean H3K4me3") +
    annotate("text", label=rlab,
             x=as.POSIXct("2013-08-15"), y=yceiling-yrange*0.08, 
             size = 6/ggplot2::.pt, hjust = 0) +
    annotate(geom = "text", label = "2012", 
             x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013", 
             x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    geom_segment(x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.2, linewidth = 0.2)
})

# PS18
local({
  i <- 18
  df_plot_rev <- df_plot %>% filter(PSnum == i)
  xmax <- max(df_temp2$date, na.rm = T)
  xmin <- min(df_temp2$date, na.rm = T)
  ymax <- max(df_plot_rev[,3:5], na.rm = T)
  ymin <- min(df_plot_rev[,3:5], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  yceiling <-  yceiling + yrange * 0.1
  
  # Scale adjustment
  scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
  df_temp2$temp_rna <- scale_to_rna(df_temp2$dailytemp)
  
  # Setting labels
  labelx<-c("","11","12",
            "1","2","3","4","5","6","7","8","9","10","")
  length(labelx)
  PSname <- df_merge2$PSname[min(which(df_merge2$PSnum==i))]
  
  # correlation
  cor_res <- cor.test(df_plot_rev$mean, temp_sample, method = "pearson")
  mark <- ifelse(cor_res$p.value < 0.0001, "****", 
                 ifelse(cor_res$p.value < 0.001, "***", 
                        ifelse(cor_res$p.value < 0.01, "**",
                               ifelse(cor_res$p.value < 0.05, "*", "ns"))))
  r <- format(cor_res$estimate, digits=2, nsmall = 2)
  rlab <- bquote(paste(italic(r), " = ", .(r), .(mark), sep=""))
  
  glist[[i]] <<- ggplot() +
    geom_line(data = df_temp2, aes(x = date, y = temp_rna), col = "gray50", linewidth = 0.2) + 
    geom_ribbon(data = df_plot_rev, aes(x = date, ymin = mean_minus_se, ymax = mean_plus_se),
                alpha = 0.3, fill = "gray30") +
    geom_line(data = df_plot_rev, aes(x = date, y = mean), color = "black", linewidth = 0.5) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_rna_rev(.), name = "")) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    theme_test(base_size = 6) +
    theme(plot.title = element_text(size = 6, face = "bold"),
          axis.title=element_text(size = 6, colour = "black"), 
          axis.title.y.right = element_blank(), 
          axis.text=element_text(size = 6, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "",
         y = "Mean H3K4me3") +
    annotate("text", label=rlab,
             x=as.POSIXct("2013-08-15"), y=yceiling-yrange*0.08, 
             size = 6/ggplot2::.pt, hjust = 0) +
    annotate(geom = "text", label = "2012", 
             x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013", 
             x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.33,
             size = 6/ggplot2::.pt, vjust=0) +
    geom_segment(x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.2, linewidth = 0.2)
})



### All PS
xmax <- max(df_temp2$date, na.rm = T)
xmin <- min(df_temp2$date, na.rm = T)
ymax <- max(df_plot[,3], na.rm = T)
ymin <- min(df_plot[,3], na.rm = T)
yrange <- (ymax - ymin)
yceiling <-  ymax + yrange * 0.05
yfloor <- ymin - yrange * 0.05

# Scale adjustment
scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
df_temp2$temp_rna <- scale_to_rna(df_temp2$dailytemp)

# Setting labels
labelx<-c("","11","12",
          "1","2","3","4","5","6","7","8","9","10","")

g_K4 <- ggplot(df_plot, aes(x = date, y = mean, group = PSnum, color = PSnum)) +
  geom_line() +
  scale_color_viridis(name="PS", option="turbo") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                     expand = c(0, 0)) +
  scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
  coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
  theme_test(base_size = 6) +
  theme(plot.title = element_text(size = 6, face = "bold"),
        axis.title=element_text(size = 6, colour = "black"), 
        axis.text=element_text(size = 6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin=unit(c(0, 6.5, 0.5, 0.5),"mm"),
        plot.tag = element_text(size = 8, face = "bold"),
        legend.title=element_text(size=6, color = "black"),
        legend.text=element_text(size=6, color = "black"),
        legend.position=c(1.13, 0.5),
        legend.key.height=unit(3, "mm"),
        legend.key.width = unit(1.5, "mm"),
        legend.background = element_rect(fill="transparent"),
        legend.spacing.x = unit(0.7, 'mm')) +
  labs(title = "H3K4me3",
    x = "",
    y = "Mean level") +
  annotate(geom = "text", label = "2012",
           x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.41,
           size = 6/ggplot2::.pt, vjust=0) +
  annotate(geom = "text", label = "2013",
           x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.41,
           size = 6/ggplot2::.pt, vjust=0) +
  geom_segment(x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.24,
               xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.24, linewidth = 0.2,
               colour = "black") +
  geom_segment(x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.24,
               xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.24, linewidth = 0.2,
               colour = "black")



### Integration of plots
void <- ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm"))

gall <- {(g_ts_K27 + labs(tag = "A")) + (g_ts_K4 + labs(tag = "C"))} /
  {(g_K27 + labs(tag = "B")) + (g_K4 + labs(tag = "D"))} +
  plot_layout(heights = c(1, 1))

ggsave(paste0(out, "Fig.3_HAI_seasonal.pdf"),
       gall, width = 100, height = 55, units = "mm")

