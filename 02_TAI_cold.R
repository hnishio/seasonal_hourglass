
# Load packages
library(tidyverse)
library(patchwork)
library(data.table)

# Create the output directory
out <- "02_TAI_cold/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}



##### Raw data (rep mean, maxRNA > 2, including warm after cold) #####

# Load data
log2rpkm_s <- as.data.frame(fread("data/RNAseq_4wkcold_rpkm_log2.tsv"))
dim(log2rpkm_s)
names(log2rpkm_s)
names(log2rpkm_s)[1] <- "Ahal_ID"
Ahal_ID <- str_split_fixed(log2rpkm_s$Ahal_ID, "\\.", 2)[,1]
log2rpkm_s$Ahal_ID <- Ahal_ID

# Modify data
log2rpkm_s2 <- log2rpkm_s %>%
  select(Ahal_ID, 
         `4wC0_0829_1`, `4wC0_0829_2`, `4wC0_0829_3`, `4wC0_0829_4`,
         `4wC2_0912_1`, `4wC2_0912_2`, `4wC2_0912_3`, `4wC2_0912_4`,
         `4wC4_0926_1`, `4wC4_0926_2`, `4wC4_0926_3`, `4wC4_0926_4`,
         `C4W1_1003_1`, `C4W1_1003_2`, `C4W1_1003_3`, `C4W1_1003_4`)

# Calculate mean of replicates
log2rpkm_s3 <- log2rpkm_s2 %>%
  mutate(mean_C0 = (`4wC0_0829_1` + `4wC0_0829_2` + `4wC0_0829_3` + `4wC0_0829_4`) / 4,
         mean_C2 = (`4wC2_0912_1` + `4wC2_0912_2` + `4wC2_0912_3` + `4wC2_0912_4`) / 4,
         mean_C4 = (`4wC4_0926_1` + `4wC4_0926_2` + `4wC4_0926_3` + `4wC4_0926_4`) / 4,
         mean_C4W1 = (`C4W1_1003_1` + `C4W1_1003_2` + `C4W1_1003_3` + `C4W1_1003_4`) / 4)

# Extract genes with maxRNA > 2
max_log2rpkm <- apply(log2rpkm_s3[,c("mean_C0", "mean_C2", "mean_C4", "mean_C4W1")], 1, max)
log2rpkm_s4 <- log2rpkm_s3 %>%
  mutate(max_log2rpkm = max_log2rpkm) %>%
  filter(max_log2rpkm > 2)

# Load orthomap data
df_orthomap <- as.data.frame(fread(paste0("data/", "orthomap_main_rev_250416.tsv")))

# Merge rpkm and orthomap
df_merge <- inner_join(log2rpkm_s4, df_orthomap, by = "Ahal_ID") %>%
  select(-PSname, -PScontinuity) %>%
  relocate(Ahal_ID, PSnum)
df_merge2 <- inner_join(log2rpkm_s4, df_orthomap, by = "Ahal_ID") %>%
  select(-PScontinuity) %>%
  relocate(Ahal_ID, PSnum, PSname)
fwrite(df_merge2, file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_repmean2.csv"), row.names=F)


# Calculate TAI
vec_tai <- NULL
for(i in 1:4){
  vec_tai[i] <- sum(df_merge$PSnum*df_merge[,i+18]) / sum(df_merge[,i+18])
}
df_tai <- data.frame(tai = vec_tai, date = c("0 wk", "2 wk", "4 wk", "4 wk + 1 wk"))

# bootstrapping
bootstrap <- matrix(NA,nrow=100,ncol=4)
for(j in 1:4){
  for(i in 1:100){
    sampleNum <- sample(nrow(df_merge), nrow(df_merge), replace=T)
    df_merge_sub <- df_merge[sampleNum,]
    bootstrap[i,j] <- sum(df_merge_sub$PSnum*df_merge_sub[,j+18]) / sum(df_merge_sub[,j+18])
  }
}

boot.quant <- matrix(NA,nrow=2,ncol=4)
for(i in 1:4){
  boot.quant[,i] <- quantile(bootstrap[,i],p=c(0.025,0.975))
}

boot.sd <- matrix(NA,nrow=2,ncol=4)
for(i in 1:4){
  mean_b <- mean(bootstrap[,i])
  sd_b <- sd(bootstrap[,i])
  boot.sd[,i] <- c(mean_b - sd_b, mean_b + sd_b)
}

# Save output
df_tai <- data.frame(tai = vec_tai, 
                     q2.5 = boot.quant[1,], 
                     q97.5 = boot.quant[2,], 
                     mean_minus_sd = boot.sd[1,], 
                     mean_plus_sd = boot.sd[2,], 
                     date = c("0 wk", "2 wk", "4 wk", "4 wk + 1 wk"))
fwrite(df_tai, file=paste0(out, "TAI_cold_rawdata_maxRNAover2_repmean2.csv"), row.names=F)



# Raw data per replicate (maxRNA > 2)
# Calculate TAI
vec_tai <- NULL
for(i in 1:16){
  vec_tai[i] <- sum(df_merge$PSnum*df_merge[,i+2]) / sum(df_merge[,i+2])
}
df_tai2 <- data.frame(tai = vec_tai, date = c(rep("0 wk",4), rep("2 wk",4), rep("4 wk",4), rep("4 wk + 1 wk",4)))

fwrite(df_tai2, file=paste0(out, "TAI_cold_rawdata_maxRNAover2_perrep2.csv"), row.names=F)





#### Plotting rep MEAN +- SE #####
se <- function(x){sd(x) / sqrt(length(x))}
df_tai2 <- as.data.frame(fread(paste0(out, "TAI_cold_rawdata_maxRNAover2_perrep2.csv")))
df_tai2$time <- as.POSIXct(c(rep("2022-08-29 12:00:00 JST",4), rep("2022-09-12 12:00:00 JST",4),
                             rep("2022-09-26 12:00:00 JST",4), rep("2022-10-03 12:00:00 JST",4)))
df_tai2$days <- difftime(df_tai2$time, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")
df_tai3 <- df_tai2 %>%
  group_by(days) %>%
  summarize(mean_tai = mean(tai),
            mean_minus_se_tai = mean(tai) - se(tai),
            mean_plus_se_tai = mean(tai) + se(tai))

# Tukey-Kramer test
df_test <- data.frame(data = df_tai2$tai, group = df_tai2$days)
df_test$group <- factor(df_test$group)

res <- aov(data ~ group, data = df_test)
tuk <- glht(res, linfct = mcp(group = "Tukey"))
tuk_cld <- cld(tuk, decreasing = F)
mark <- as.character(tuk_cld[["mcletters"]][["Letters"]])

# Load temp setting
df_tai <- df_tai2 %>% group_by(time) %>% summarise(mean = mean(tai))
df_temp_warm <- as.data.frame(fread("data/warm_temp.csv"))
df_temp_cold <- as.data.frame(fread("data/cold_temp.csv"))
df_temp_warm$time <- force_tz(mdy_hms(df_temp_warm$time), tzone = "Asia/Tokyo")
df_temp_cold$time <- force_tz(mdy_hms(df_temp_cold$time), tzone = "Asia/Tokyo")
df_temp_warm1 <- df_temp_warm %>%
  filter(time >= (df_tai$time[1] - 60*60*24*7) & time <= df_tai$time[1])
df_temp_cold1 <- df_temp_cold %>%
  filter(time >= df_tai$time[1] & time <= df_tai$time[3])
df_temp_warm2 <- df_temp_warm %>%
  filter(time >= df_tai$time[3] & time <= df_tai$time[4] + 60*60*24*7)
df_temp <- rbind(df_temp_warm1, df_temp_cold1, df_temp_warm2)
df_temp$days <- difftime(df_temp$time, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")

# Plotting
local({
  
  # Setting var_nameiables
  xmax <- max(df_temp$days, na.rm = T)
  xmin <- min(df_temp$days, na.rm = T)
  ymax <- max(df_tai3[,2:4], na.rm = T)
  ymin <- min(df_tai3[,2:4], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_tai <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_tai_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
  df_temp$temp_tai <- scale_to_tai(df_temp$temp)
  
  g_tai <<- ggplot() +
    geom_line(data = df_temp, aes(x = days, y = temp_tai), col = "gray50", linewidth = 0.2) + 
    # geom_hline(yintercept = 0, col='gray50', linetype='dashed') +
    # geom_pointrange(data = df_tai, aes(x=days, y=tai, ymin=mean_minus_se, ymax=mean_plus_se), 
    #                 color='grey5', fill='grey95', size=0.3, fatten = 0.2) +
    geom_point(data = df_tai3, aes(x=days, y=mean_tai), size=0.8) +
    geom_errorbar(data = df_tai3, aes(x=days, ymin=mean_minus_se_tai, ymax=mean_plus_se_tai), width = 1.2, linewidth = 0.5) +
    # geom_jitter(data = df_tai2, aes(x = days, y = tai), color = "black", alpha = 1, size = 0.5, width = 0.25) +
    annotate(geom = "text", x=df_tai3$days, y=rep(yceiling+yrange*0.07,length(mark)), label=mark, size=7/ggplot2::.pt) +
    annotate(geom = "text", label = "Young", x=-8.5, y=yceiling+yrange*0.22, size=7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "Old", x=-8.5, y=yfloor-yrange*0.22, size=7/ggplot2::.pt, vjust=0) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling+yrange*0.15), clip='off') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_tai_rev(.), name = "Temperature (°C)")) +
    theme_test(base_size = 7) +
    theme(plot.title=element_blank(), 
          axis.title=element_text(size = 7, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text=element_text(size = 7, color = "black"),
          plot.margin = unit(c(1, 1, 1, 1),"mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = "",
         x = "Days after transfer to cold",
         y = "Transcriptome age index")
})





##### PS1 and PS18 with temperature #####

### Comparison between PSs
df_merge <- as.data.frame(fread(file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_repmean2.csv")))
df_merge <- df_merge[,1:19]
table(df_merge$PSnum)

# Calculate mean expression for each PS at each date
df_eachPS <- df_merge %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), mean, .names = "{.col}"))

names(df_eachPS) <- c("PSnum", rep("C0",4), rep("C2",4), rep("C4",4), rep("C4W1",4))
df_eachPS_long <- df_eachPS %>%
  pivot_longer(cols = -PSnum, names_to = "condition", values_to = "mean_allgenes")
fwrite(df_eachPS_long, paste0(out, "PS_cold_log2rpkm_maxRNAover2_allgenesmean.csv"))

# Mean of replicates
df_eachPS_meanreps <- df_eachPS_long %>%
  group_by(PSnum, condition) %>%
  summarize(mean_rep = mean(mean_allgenes))

# SE of replicates
se <- function(x){sd(x) / sqrt(length(x))}
df_eachPS_sereps <- df_eachPS_long %>%
  group_by(PSnum, condition) %>%
  summarize(se_rep = se(mean_allgenes))

df_eachPS_minusse <- cbind(df_eachPS_meanreps[,1:2], (df_eachPS_meanreps[,3] - df_eachPS_sereps[,3]))
df_eachPS_plusse <- cbind(df_eachPS_meanreps[,1:2], (df_eachPS_meanreps[,3] + df_eachPS_sereps[,3]))

# Modify data frames
df_plot <- inner_join(df_eachPS_meanreps, df_eachPS_minusse, by = c("PSnum", "condition"))
df_plot <- inner_join(df_plot, df_eachPS_plusse, by = c("PSnum", "condition"))
df_plot$date <- NA
df_plot$date[df_plot$condition == "C0"] <- "2022-08-29 12:00:00 JST"
df_plot$date[df_plot$condition == "C2"] <- "2022-09-12 12:00:00 JST"
df_plot$date[df_plot$condition == "C4"] <- "2022-09-26 12:00:00 JST"
df_plot$date[df_plot$condition == "C4W1"] <- "2022-10-03 12:00:00 JST"
df_plot$date <- as.POSIXct(df_plot$date)
df_plot$days <- difftime(df_plot$date, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")
names(df_plot) <- c("PSnum", "condition", "mean", "mean_minus_se", "mean_plus_se", "date", "days")
fwrite(df_plot, paste0(out, "PS_cold_log2rpkm_maxRNAover2_repmean.csv"))


### Plotting
df_merge <- as.data.frame(fread(file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_repmean2.csv")))
df_eachPS_long <- as.data.frame(fread(paste0(out, "PS_cold_log2rpkm_maxRNAover2_allgenesmean.csv")))
names(df_eachPS_long)
dim(df_eachPS_long)

# Load data for plot
df_plot <- as.data.frame(fread(paste0(out, "PS_cold_log2rpkm_maxRNAover2_repmean.csv")))

# Plotting
glist <- list()

# PS1
local({
  i = 1
  
  # Tukey-Kramer test
  df_rev <- df_eachPS_long %>%
    filter(PSnum == i)
  df_test <- data.frame(data = df_rev$mean_allgenes, group = df_rev$condition)
  df_test$group <- factor(df_test$group)
  res <- aov(data ~ group, data = df_test)
  tuk <- glht(res, linfct = mcp(group = "Tukey"))
  tuk_cld <- cld(tuk, decreasing = F)
  mark <- as.character(tuk_cld[["mcletters"]][["Letters"]])
  mark <- c("c", "a", "b", "c")
  
  df_plot_rev <- df_plot %>% filter(PSnum == i)
  xmax <- max(df_temp$days, na.rm = T)
  xmin <- min(df_temp$days, na.rm = T)
  ymax <- max(df_plot_rev[,3:5], na.rm = T)
  ymin <- min(df_plot_rev[,3:5], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
  df_temp$temp_rna <- scale_to_rna(df_temp$temp)
  
  # Setting labels
  PSname <- df_merge$PSname[min(which(df_merge$PSnum==i))]
  
  # Plotting
  glist[[i]] <<- ggplot() +
    geom_line(data = df_temp, aes(x = days, y = temp_rna), col = "gray50", linewidth = 0.2) + 
    geom_point(data = df_plot_rev, aes(x = days, y = mean), size=0.8) +
    geom_errorbar(data = df_plot_rev, aes(x=days, ymin=mean_minus_se, ymax=mean_plus_se), width = 2, linewidth = 0.5) +
    annotate(geom = "text", x=df_plot_rev$days, y=rep(yceiling+yrange*0.07,length(mark)), label=mark, size=7/ggplot2::.pt) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling+yrange*0.15), clip='off') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_rna_rev(.), name = "")) +
    theme_test(base_size = 7) +
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.title=element_text(size = 7, colour = "black"), 
          axis.title.y.right=element_blank(), 
          axis.text=element_text(size = 7, colour = "black"),
          axis.ticks = element_line(color = "black"),
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "Days after transfer to cold",
         y = "Mean expression")
})

# PS18
local({
  i = 18
  
  # Tukey-Kramer test
  df_rev <- df_eachPS_long %>%
    filter(PSnum == i)
  df_test <- data.frame(data = df_rev$mean_allgenes, group = df_rev$condition)
  df_test$group <- factor(df_test$group)
  res <- aov(data ~ group, data = df_test)
  tuk <- glht(res, linfct = mcp(group = "Tukey"))
  tuk_cld <- cld(tuk, decreasing = F)
  mark <- as.character(tuk_cld[["mcletters"]][["Letters"]])
  
  df_plot_rev <- df_plot %>% filter(PSnum == i)
  xmax <- max(df_temp$days, na.rm = T)
  xmin <- min(df_temp$days, na.rm = T)
  ymax <- max(df_plot_rev[,3:5], na.rm = T)
  ymin <- min(df_plot_rev[,3:5], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
  df_temp$temp_rna <- scale_to_rna(df_temp$temp)
  
  # Setting labels
  PSname <- df_merge$PSname[min(which(df_merge$PSnum==i))]
  
  # Plotting
  glist[[i]] <<- ggplot() +
    geom_line(data = df_temp, aes(x = days, y = temp_rna), col = "gray50", linewidth = 0.2) + 
    geom_point(data = df_plot_rev, aes(x = days, y = mean), size=0.8) +
    geom_errorbar(data = df_plot_rev, aes(x=days, ymin=mean_minus_se, ymax=mean_plus_se), width = 2, linewidth = 0.5) +
    annotate(geom = "text", x=df_plot_rev$days, y=rep(yceiling+yrange*0.07,length(mark)), label=mark, size=7/ggplot2::.pt) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling+yrange*0.15), clip='off') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_rna_rev(.), name = "Temperature (°C)")) +
    theme_test(base_size = 7) +
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.title=element_text(size = 7, colour = "black"), 
          axis.title.y.left=element_blank(), 
          axis.text=element_text(size = 7, colour = "black"),
          axis.ticks = element_line(color = "black"),
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "Days after transfer to cold",
         y = "")
})



### Integration of plots
void <- ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm"))
gall <- (g_tai + labs(tag = "A")) / 
  {(glist[[1]] + labs(tag = "B")) + glist[[18]]} +
  plot_layout(heights = c(1, 1))

ggsave(paste0(out, "Fig.2_TAI_cold_250503.pdf"),
       gall, width = 100, height = 80, units = "mm")

