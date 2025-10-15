
# Load packages
library(tidyverse)
library(patchwork)
library(data.table)
library(viridis)
library(ggbiplot)

# Create output directory
out <- "04_MADS_seasonal/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}



### MADS genes
df_tf <- as.data.frame(fread(paste0("data/Ath_TF_list.txt")))
df_tf$Family <- str_replace(df_tf$Family, "/", ";")
tf_names <- sort(unique(df_tf$Family))

glist_all <- list()
min_dates <- NULL
name_tf <- NULL
mean_levels <- NULL
for(t in 1:length(tf_names)){
  
  ### Plotting
  if(!file.exists(paste0("23_MADS_RNA_251013/", tf_names[t], "_RNA_log2rpkm_maxRNAover2.csv"))){
    next
  }
  df_stat <- as.data.frame(fread(paste0("23_MADS_RNA_251013/", tf_names[t], "_RNA_log2rpkm_maxRNAover2.csv")))
  df_stat$date <- force_tz(df_stat$date, tzone = "Asia/Tokyo")
  
  # Load Nishiwaki dailytemp
  df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
  df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
  df_temp2 <- df_temp %>%
    filter(date >= df_stat$date[1] & date <= df_stat$date[nrow(df_stat)])
  
  # Setting var_nameiables
  labelx<-c("","7","","9","","11","",
            "1","","3","","5","","7","","9","","11","",
            "1","","3","","5","","7","")
  length(labelx)
  
  local({
    xmax <- max(df_stat$date, na.rm = T)
    xmin <- min(df_stat$date, na.rm = T)
    ymax <- max(df_stat$mean+df_stat$se, na.rm = T)
    ymin <- min(df_stat$mean-df_stat$se, na.rm = T)
    yrange <- (ymax - ymin)
    yceiling <-  ymax + yrange * 0.05
    yfloor <- ymin - yrange * 0.05
    
    # Scale adjustment
    scale_to_stat <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
    scale_to_stat_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
    df_temp2$temp_stat <- scale_to_stat(df_temp2$dailytemp)
    
    glist_all[[t]] <<- ggplot() +
      geom_line(data = df_temp2, aes(x = date, y = temp_stat), col = "gray50", linewidth = 0.2) + 
      geom_ribbon(data = df_stat, aes(x = date, ymin = mean-se, ymax = mean+se),
                  alpha = 0.3, fill = "gray30") +
      geom_line(data = df_stat, aes(x = date, y = mean), col = "black", linewidth = 0.6) +
      # geom_point(data = df_stat2, aes(x = date, y = stat), col = "black", alpha = 0.6, size = 0.2) +
      # scale_x_continuous(breaks = 1:12) +
      scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
      coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
      # annotate("text", 
      #          x=as.POSIXct("2011-12-01"), y=yceiling-yrange*0.08,
      #          label = bquote(italic(p) < .(format(p_val, scientific = FALSE))), 
      #          size=6/ggplot2::.pt) +
      theme_test(base_size=6) +
      theme(plot.title = element_text(size=6, colour = "black", face = "bold"),
            axis.title=element_text(size=6, colour = "black"), 
            axis.text=element_text(size=6, colour = "black"),
            axis.ticks = element_line(color = "black"),
            legend.position = "none",
            plot.margin = unit(c(1.5,0.2,1,0.2), "mm"),
            plot.tag = element_text(size = 8, face = "bold")) +
      labs(title = paste0("mRNA (", tf_names[t], ")"), 
           x = "", y = "Mean level") +
      scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                         expand = c(0, 0),
                         sec.axis = sec_axis(~ scale_to_stat_rev(.), name = "Temperature (°C)")) +
      annotate(geom = "text", label = "2011", 
               x=as.POSIXct("2011-09-25"), y=yfloor-yrange*0.27,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2012", 
               x=as.POSIXct("2012-07-01"), y=yfloor-yrange*0.27,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2013", 
               x=as.POSIXct("2013-04-10"), y=yfloor-yrange*0.27,
               size=6/ggplot2::.pt, vjust=0) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2011-6-26"), y=yfloor-yrange*0.16, 
                   xend=as.POSIXct("2011-12-25"), yend=yfloor-yrange*0.16, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2012-01-07"), y=yfloor-yrange*0.16, 
                   xend=as.POSIXct("2012-12-25"), yend=yfloor-yrange*0.16, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2013-01-07"), y=yfloor-yrange*0.16, 
                   xend=as.POSIXct("2013-07-06"), yend=yfloor-yrange*0.16, linewidth = 0.2)
  })
  
  min_dates <- c(min_dates, str_sub(df_stat$date[which.min(df_stat$mean)], 1, 10))
  name_tf <- c(name_tf, tf_names[t])
  mean_levels <- rbind(mean_levels, df_stat$mean)
}

# Remove TF families with little seasonal changes
range_c <- function(x){max(x)-min(x)}
tf_range <- apply(mean_levels, 1, range_c)
data.frame(tf_names = name_tf, range = tf_range)
thre <- quantile(tf_range, 0.25)
mean_levels2 <- mean_levels[tf_range > thre,]
name_tf2 <- name_tf[tf_range > thre]

# Standardaization
scale_c <- function(x){(x-mean(x))/sd(x)}
test <- apply(mean_levels2, 1, scale_c)
df_mean_levels <- as.data.frame(test)
names(df_mean_levels) <- name_tf2
row.names(df_mean_levels) <- paste0("t", 1:nrow(df_mean_levels))

# PCA analysis RNA
df_RNA_pca <- t(df_mean_levels)
RNA_pca_model <- prcomp(df_RNA_pca, center = T, scale=T)

summary(RNA_pca_model)
xmin <- -9
xmax <- 15
ymin <- -7
ymax <- 15
gpca_RNA <- ggbiplot(RNA_pca_model, obs.scale=1, var.scale=1,
                 labels=name_tf2, ellipse=F, circle=T, alpha=0.5, var.axes=F, 
                 labels.size=5/ggplot2::.pt) +
  coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax)) +
  # scale_colour_manual(name="Group", values= c(rgb(1,0.65,0), rgb(0.4,0.8,1), rgb(0,0,0.5), rgb(0.55,0,0)))+
  theme_bw()+
  theme(legend.title=element_text(size=6, colour = "black"),
        legend.text=element_text(size=6, colour = "black"),
        plot.title=element_text(size=6, colour = "black", face="bold", vjust = -0.7),
        axis.title=element_text(size=6, colour = "black"), 
        axis.text=element_text(size=6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 8, face = "bold")) + 
  labs(title="PCA of TF families based on mRNA")
ggsave(paste0(out, "PCA_RNA_251014.pdf"), 
       gpca_RNA, height = 60, width = 60, units = "mm")





##### H3K27me3 #####
glist_K27_all <- list()
min_dates <- NULL
name_tf <- NULL
mean_levels <- NULL
for(t in 1:length(tf_names)){
  
  ### Plotting
  if(!file.exists(paste0("22_MADS_K4K27_251012/K27/", tf_names[t], "_K27_log2rpkm_maxRNAover2.csv"))){
    next
  }
  df_stat <- as.data.frame(fread(paste0("22_MADS_K4K27_251012/K27/", tf_names[t], "_K27_log2rpkm_maxRNAover2.csv")))
  df_stat$date <- with_tz(df_stat$date, tzone = "Asia/Tokyo")
  
  # Load Nishiwaki dailytemp
  df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
  df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
  df_temp2 <- df_temp %>%
    filter(date >= df_stat$date[1] & date <= df_stat$date[nrow(df_stat)])
  
  # Setting var_nameiables
  labelx<-c("","11","12",
            "1","2","3","4","5","6","7","8","9","10","")
  length(labelx)
  
  local({
    xmax <- max(df_stat$date, na.rm = T)
    xmin <- min(df_stat$date, na.rm = T)
    ymax <- max(df_stat$mean+df_stat$se, na.rm = T)
    ymin <- min(df_stat$mean-df_stat$se, na.rm = T)
    yrange <- (ymax - ymin)
    yceiling <-  ymax + yrange * 0.05
    yfloor <- ymin - yrange * 0.05
    
    # Scale adjustment
    scale_to_stat <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
    scale_to_stat_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
    df_temp2$temp_stat <- scale_to_stat(df_temp2$dailytemp)
    
    glist_K27_all[[t]] <<- ggplot() +
      geom_line(data = df_temp2, aes(x = date, y = temp_stat), col = "gray50", linewidth = 0.2) + 
      geom_ribbon(data = df_stat, aes(x = date, ymin = mean-se, ymax = mean+se),
                  alpha = 0.3, fill = "gray30") +
      geom_line(data = df_stat, aes(x = date, y = mean), col = "black", linewidth = 0.6) +
      # geom_point(data = df_stat2, aes(x = date, y = stat), col = "black", alpha = 0.6, size = 0.2) +
      # scale_x_continuous(breaks = 1:12) +
      scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
      coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
      # annotate("text", 
      #          x=as.POSIXct("2011-12-01"), y=yceiling-yrange*0.08,
      #          label = bquote(italic(p) < .(format(p_val, scientific = FALSE))), 
      #          size=6/ggplot2::.pt) +
      theme_test(base_size=6) +
      theme(plot.title = element_text(size = 6, colour = "black", face = "bold"),
            axis.title=element_text(size=6, colour = "black"), 
            axis.text=element_text(size=6, colour = "black"),
            axis.ticks = element_line(color = "black"),
            legend.position = "none",
            plot.margin = unit(c(1.5,0.2,1,0.2), "mm"),
            plot.tag = element_text(size = 8, face = "bold")) +
      labs(title = paste0("H3K27me3 (", tf_names[t], ")"), 
           x = "", y = "Mean level") +
      scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                         expand = c(0, 0),
                         sec.axis = sec_axis(~ scale_to_stat_rev(.), name = "Temperature (°C)")) +
      annotate(geom = "text", label = "2012", 
               x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.27,
               size = 6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2013", 
               x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.27,
               size = 6/ggplot2::.pt, vjust=0) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.16, 
                   xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.16, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.16, 
                   xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.16, linewidth = 0.2)
    
  })
  
  min_dates <- c(min_dates, str_sub(df_stat$date[which.min(df_stat$mean)], 1, 10))
  name_tf <- c(name_tf, tf_names[t])
  mean_levels <- rbind(mean_levels, df_stat$mean)
}

data.frame(tf_names = name_tf, date = min_dates)
tf_sd <- apply(mean_levels, 1, sd)
data.frame(tf_names = name_tf, sd = tf_sd)

# Remove TF families with little seasonal changes
range_c <- function(x){max(x)-min(x)}
tf_range <- apply(mean_levels, 1, range_c)
data.frame(tf_names = name_tf, range = tf_range)
thre <- quantile(tf_range, 0.25)
mean_levels2 <- mean_levels[tf_range > thre,]
name_tf2 <- name_tf[tf_range > thre]

# library(season)
# mat_cos <- matrix(nrow=nrow(mean_levels), ncol=9)
# for(i in 1:nrow(mean_levels)){
#   df_cos <- data.frame(date=as.Date(df_stat$date), mean=mean_levels[i,])
#   res_cos <- cosinor(mean~1, date="date", data=df_cos)
#   phase <- gsub("Month = ", "2012-", summary(res_cos)$phase)
#   phase <- gsub("月 , day = ", "-", phase)
#   peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
#   lphase <- gsub("Month = ", "2012-", summary(res_cos)$lphase)
#   lphase <- gsub("月 , day = ", "-", lphase)
#   trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
#   mat_cos[i,] <- c(summary(res_cos)$amp*2,peak,trough,
#                summary(res_cos$glm)$coefficients[,"Estimate"],
#                summary(res_cos$glm)$coefficients[,"Pr(>|t|)"])
# }
# colnames(mat_cos) <- c("amp_cos","peak_cos","trough_cos",
#                       "coef_inter_cos","coef_cosw_cos","coef_sinw_cos",
#                       "pval_inter_cos","pval_cosw_cos","pval_sinw_cos")
# df_cos <- cbind(data.frame(tf = name_tf), as.data.frame(mat_cos))
# df_cos_0.05 <- df_cos[(df_cos$pval_cosw_cos<0.025 | df_cos$pval_sinw_cos<0.025),]

# Standardaization
scale_c <- function(x){(x-mean(x))/sd(x)}
test <- apply(mean_levels2, 1, scale_c)
df_mean_levels <- as.data.frame(test)
names(df_mean_levels) <- name_tf2
row.names(df_mean_levels) <- paste0("t", 1:12)


# PCA analysis H3K27me3
# Preparation of labels of plot
# mylabel <- c("C0_rep1","C0_rep2","C0_rep3",
#              "C2_rep1","C2_rep2","C2_rep3",
#              "C4_rep1","C4_rep2","C4_rep3",
#              "C4W1_rep1","C4W1_rep2","C4W1_rep3")
# mygroup <- c(rep("C0",3), rep("C2",3),
#              rep("C4",3),rep("C4W1",3))
# mygroup <- factor(mygroup, levels=c("C0", "C2", "C4", "C4W1"))

df_H3K27me3_pca <- t(df_mean_levels)
H3K27me3_pca_model <- prcomp(df_H3K27me3_pca, center = T, scale=T)

summary(H3K27me3_pca_model)
xmin <- -15
xmax <- 3
ymin <- -4
ymax <- 7
gpca_K27 <- ggbiplot(H3K27me3_pca_model, obs.scale=1, var.scale=1,
                 labels=name_tf2, ellipse=F, circle=T, alpha=0.5, var.axes=F, 
                 labels.size=5/ggplot2::.pt) +
  coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax)) +
  # scale_colour_manual(name="Group", values= c(rgb(1,0.65,0), rgb(0.4,0.8,1), rgb(0,0,0.5), rgb(0.55,0,0)))+
  theme_bw()+
  theme(legend.title=element_text(size=6, colour = "black"),
        legend.text=element_text(size=6, colour = "black"),
        plot.title=element_text(size=6, colour = "black", face="bold", vjust = -0.7),
        axis.title=element_text(size=6, colour = "black"), 
        axis.text=element_text(size=6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 8, face = "bold")) + 
  labs(title="PCA of TF families based on H3K27me3")
ggsave(paste0(out, "PCA_H3K27me3_251014.pdf"), 
       gpca_K27, height = 60, width = 60, units = "mm")







##### H3K4me3 #####
glist_K4_all <- list()
min_dates <- NULL
name_tf <- NULL
mean_levels <- NULL
for(t in 1:length(tf_names)){
  
  ### Plotting
  if(!file.exists(paste0("22_MADS_K4K27_251012/K4/", tf_names[t], "_K4_log2rpkm_maxRNAover2.csv"))){
    next
  }
  df_stat <- as.data.frame(fread(paste0("22_MADS_K4K27_251012/K4/", tf_names[t], "_K4_log2rpkm_maxRNAover2.csv")))
  df_stat$date <- with_tz(df_stat$date, tzone = "Asia/Tokyo")
  
  # Load Nishiwaki dailytemp
  df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
  df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
  df_temp2 <- df_temp %>%
    filter(date >= df_stat$date[1] & date <= df_stat$date[nrow(df_stat)])
  
  # Setting var_nameiables
  labelx<-c("","11","12",
            "1","2","3","4","5","6","7","8","9","10","")
  length(labelx)
  
  local({
    xmax <- max(df_stat$date, na.rm = T)
    xmin <- min(df_stat$date, na.rm = T)
    ymax <- max(df_stat$mean+df_stat$se, na.rm = T)
    ymin <- min(df_stat$mean-df_stat$se, na.rm = T)
    yrange <- (ymax - ymin)
    yceiling <-  ymax + yrange * 0.05
    yfloor <- ymin - yrange * 0.05
    
    # Scale adjustment
    scale_to_stat <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
    scale_to_stat_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
    df_temp2$temp_stat <- scale_to_stat(df_temp2$dailytemp)
    
    glist_K4_all[[t]] <<- ggplot() +
      geom_line(data = df_temp2, aes(x = date, y = temp_stat), col = "gray50", linewidth = 0.2) + 
      geom_ribbon(data = df_stat, aes(x = date, ymin = mean-se, ymax = mean+se),
                  alpha = 0.3, fill = "gray30") +
      geom_line(data = df_stat, aes(x = date, y = mean), col = "black", linewidth = 0.6) +
      # geom_point(data = df_stat2, aes(x = date, y = stat), col = "black", alpha = 0.6, size = 0.2) +
      # scale_x_continuous(breaks = 1:12) +
      scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
      coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
      # annotate("text", 
      #          x=as.POSIXct("2011-12-01"), y=yceiling-yrange*0.08,
      #          label = bquote(italic(p) < .(format(p_val, scientific = FALSE))), 
      #          size=6/ggplot2::.pt) +
      theme_test(base_size=6) +
      theme(plot.title = element_text(size = 6, colour = "black", face = "bold"),
            axis.title=element_text(size=6, colour = "black"), 
            axis.text=element_text(size=6, colour = "black"),
            axis.ticks = element_line(color = "black"),
            legend.position = "none",
            plot.margin = unit(c(1.5,0.2,1,0.2), "mm"),
            plot.tag = element_text(size = 8, face = "bold")) +
      labs(title = paste0("H3K4me3 (", tf_names[t], ")"), 
           x = "", y = "Mean level") +
      scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                         expand = c(0, 0),
                         sec.axis = sec_axis(~ scale_to_stat_rev(.), name = "Temperature (°C)")) +
      annotate(geom = "text", label = "2012", 
               x=as.POSIXct("2012-12-01"), y=yfloor-yrange*0.27,
               size = 6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2013", 
               x=as.POSIXct("2013-05-10"), y=yfloor-yrange*0.27,
               size = 6/ggplot2::.pt, vjust=0) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2012-11-01"), y=yfloor-yrange*0.16, 
                   xend=as.POSIXct("2012-12-28"), yend=yfloor-yrange*0.16, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2013-01-04"), y=yfloor-yrange*0.16, 
                   xend=as.POSIXct("2013-09-30"), yend=yfloor-yrange*0.16, linewidth = 0.2)
    
  })
  
  min_dates <- c(min_dates, str_sub(df_stat$date[which.min(df_stat$mean)], 1, 10))
  name_tf <- c(name_tf, tf_names[t])
  mean_levels <- rbind(mean_levels, df_stat$mean)
}

# Remove TF families with little seasonal changes
range_c <- function(x){max(x)-min(x)}
tf_range <- apply(mean_levels, 1, range_c)
data.frame(tf_names = name_tf, range = tf_range)
thre <- quantile(tf_range, 0.25)
mean_levels2 <- mean_levels[tf_range > thre,]
name_tf2 <- name_tf[tf_range > thre]

# Standardaization
scale_c <- function(x){(x-mean(x))/sd(x)}
test <- apply(mean_levels2, 1, scale_c)
df_mean_levels <- as.data.frame(test)
names(df_mean_levels) <- name_tf2
row.names(df_mean_levels) <- paste0("t", 1:12)

# PCA analysis H3K4me3
df_H3K4me3_pca <- t(df_mean_levels)
H3K4me3_pca_model <- prcomp(df_H3K4me3_pca, center = T, scale=T)

summary(H3K4me3_pca_model)
xmin <- -3.3
xmax <- 4.5
ymin <- -3.3
ymax <- 3
gpca_K4 <- ggbiplot(H3K4me3_pca_model, obs.scale=1, var.scale=1,
                 labels=name_tf2, ellipse=F, circle=T, alpha=0.5, var.axes=F, 
                 labels.size=5/ggplot2::.pt) +
  coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax)) +
  # scale_colour_manual(name="Group", values= c(rgb(1,0.65,0), rgb(0.4,0.8,1), rgb(0,0,0.5), rgb(0.55,0,0)))+
  theme_bw()+
  theme(legend.title=element_text(size=6, colour = "black"),
        legend.text=element_text(size=6, colour = "black"),
        plot.title=element_text(size=6, colour = "black", face="bold", vjust = -0.7),
        axis.title=element_text(size=6, colour = "black"), 
        axis.text=element_text(size=6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1,1,1,1), "mm"),
        plot.tag = element_text(size = 8, face = "bold")) + 
  labs(title="PCA of TF families based on H3K4me3")
ggsave(paste0(out, "PCA_H3K4me3_251014.pdf"), 
       gpca_K4, height = 60, width = 60, units = "mm")





### Integration of plots
void <- ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm"))

gall <- {(gpca_RNA + labs(tag = "A")) | (gpca_K27 + labs(tag = "B")) | (gpca_K4 + labs(tag = "C"))} /
  {(glist_all[[34]] + labs(tag = "D")) +
      (glist_K27_all[[34]] + theme(axis.title.y.right =element_blank()) + labs(tag = "E")) +
      (glist_K4_all[[34]] + theme(axis.title.y.left =element_blank())) +
      plot_layout(widths = c(2, 1, 1))} +
  plot_layout(heights = c(1, 1))

ggsave(paste0(out, "Fig.4_MADS_seasonal.pdf"),
       gall, width = 170, height = 80, units = "mm")

