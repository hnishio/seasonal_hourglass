
# Load packages
library(tidyverse)
library(patchwork)
library(data.table)
library(multcomp)
library(viridis)

# Create the output directory
out <- "05_TAI_HAI_cold/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}



##### mRNA #####
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

# Load orthomap data
df_orthomap <- as.data.frame(fread(paste0("01_TAI_250403/", "orthomap_main_rev_250416.tsv")))

log2rpkm_s3 <- log2rpkm_s2 %>%
  mutate(mean_C0 = (`4wC0_0829_1` + `4wC0_0829_2` + `4wC0_0829_3` + `4wC0_0829_4`) / 4,
         mean_C2 = (`4wC2_0912_1` + `4wC2_0912_2` + `4wC2_0912_3` + `4wC2_0912_4`) / 4,
         mean_C4 = (`4wC4_0926_1` + `4wC4_0926_2` + `4wC4_0926_3` + `4wC4_0926_4`) / 4,
         mean_C4W1 = (`C4W1_1003_1` + `C4W1_1003_2` + `C4W1_1003_3` + `C4W1_1003_4`) / 4)

# Extract genes with maxRNA > 2
df_spline <- as.data.frame(fread(paste0("01_TAI_seasonal/", "log2rpkm_sp_1wk.csv"), header = T))
df_spline_maxRNAover2 <- df_spline %>%
  mutate(max_log2rpkm = do.call(pmax, c(across(where(is.numeric)), na.rm = TRUE))) %>%
  filter(max_log2rpkm > 2) %>%
  select(-max_log2rpkm)

log2rpkm_s4 <- inner_join(log2rpkm_s3, data.frame(Ahal_ID = df_spline_maxRNAover2$Ahal_ID), by = "Ahal_ID")


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



### Plotting3, rep MEAN +- SE
se <- function(x){sd(x) / sqrt(length(x))}
df_tai2 <- as.data.frame(fread(paste0(out, "TAI_cold_rawdata_maxRNAover2_perrep2.csv")))
df_tai2$time <- as.POSIXct(c(rep("2022-08-29 12:00:00 JST",4), rep("2022-09-12 12:00:00 JST",4),
                             rep("2022-09-26 12:00:00 JST",4), rep("2022-10-03 12:00:00 JST",4)))
df_tai2$days <- difftime(df_tai2$time, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")
df_tai3 <- df_tai2 %>%
  group_by(days) %>%
  summarize(mean_tai = mean(tai),
            mean_minus_se_tai = mean(tai) - se(tai),
            mean_plus_se_tai = mean(tai) + se(tai),
            se_tai = se(tai))

# Tukey-Kramer test
df_test <- data.frame(data = df_tai2$tai, group = df_tai2$days)
df_test$group <- factor(df_test$group)

res <- aov(data ~ group, data = df_test)
tuk <- glht(res, linfct = mcp(group = "Tukey"))
tuk_cld <- cld(tuk, decreasing = F)
mark <- as.character(tuk_cld[["mcletters"]][["Letters"]])

# Load temp setting
df_tai <- df_tai2 %>% group_by(time) %>% summarise(mean = mean(tai))
df_temp_warm <- as.data.frame(fread("data/room1-2_lower_rev.csv"))
df_temp_cold <- as.data.frame(fread("data/cold_1_left_rev.csv"))
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
  
  ymean <- (ymax + ymin) / 2
  ymax <- ymean + max(df_tai3[,5], na.rm = T) * 3.3
  ymin <- ymean - max(df_tai3[,5], na.rm = T) * 3.3
  
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
    annotate(geom = "text", x=df_tai3$days, y=rep(yceiling+yrange*0.07,length(mark)), label=mark, size = 6/ggplot2::.pt) +
    annotate(geom = "text", label = "Young", x=-12, y=yceiling+yrange*0.22, size = 6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "Old", x=-12, y=yfloor-yrange*0.22, size = 6/ggplot2::.pt, vjust=0) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling+yrange*0.15), clip='off') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_tai_rev(.), name = "Temperature (°C)")) +
    theme_test(base_size = 6) +
    theme(plot.title=element_text(size = 6, color = "black", face = "bold"),
          axis.title=element_text(size = 6, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text=element_text(size = 6, color = "black"),
          plot.margin = unit(c(1, 1, 1, 1),"mm"),
          plot.tag = element_text(size = 8, face = "bold")) +
    labs(title = "mRNA",
         x = "Days after transfer to cold",
         y = "Age index")
})




### All PS
df_merge <- as.data.frame(fread(file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_repmean2.csv")))
names(df_merge)
dim(df_merge)
df_merge[1:5,1:5]

df_merge2 <- df_merge[,c(1:3,20:23)]
dates <- c("2022-08-29 12:00:00 JST", "2022-09-12 12:00:00 JST",
           "2022-09-26 12:00:00 JST", "2022-10-03 12:00:00 JST")
names(df_merge2) <- c(names(df_merge)[1:3], dates)

# Calculate mean expression for each PS at each date
df_eachPS <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), mean, .names = "{.col}"))

df_plot <- df_eachPS %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean"
  ) %>%
  mutate(date = as.POSIXct(date))  
df_plot$days <- as.numeric(difftime(df_plot$date, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days"))


# Plotting
xmax <- max(df_temp$days, na.rm = T)
xmin <- min(df_temp$days, na.rm = T)
ymax <- max(df_plot[,3], na.rm = T)
ymin <- min(df_plot[,3], na.rm = T)

ymean <- (ymax + ymin) / 2
ymax <- ymean + 2
ymin <- ymean - 2

yrange <- (ymax - ymin)
yceiling <-  ymax + yrange * 0.05
yfloor <- ymin - yrange * 0.05

# Scale adjustment
scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
df_temp$temp_rna <- scale_to_rna(df_temp$temp)

# Setting labels
g_RNA <- ggplot(df_plot, aes(x = days, y = mean, group = PSnum, color = PSnum)) +
  geom_line() +
  scale_color_viridis(name="PS", option="turbo") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
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
  labs(title = "mRNA",
       x = "Days after transfer to cold",
       y = "Mean level")





##### Load ChIP-seq data #####
var_name <- c("H3K4me3", "H3K27me3")
date <- c("220829_4wC0", "220912_4wC2", "220926_4wC4", "221003_C4W1")
rep <- 1:3
inp <- "data/20250725_Ratchet_ChIP-seq_4W/"

j=1;t=1;i=1
for(j in 1:3){ #var_name
  for(t in 1:4){ #date
    for(i in 1:3){ #rep
      region <- ifelse(var_name[j] == "H3K4me3", "TSS", "GENE")
      file_name <- paste0(inp, "consensus_unmapped_", var_name[j], "_", date[t], "_rep", rep[i], "_mult_mapq30_properpair_", region, ".txt")
      if(!file.exists(file_name)){
        next
      }
      df <- as.data.frame(fread(file_name)) 
      df2 <- df[order(df$V7),] %>%
        select(V7, V6)
      names(df2) <- c("Ahal_ID", "log2rpkm")
      eval(parse(text = paste0(
        var_name[j], "_", date[t], "_rep", rep[i], "<- df2"
      )))
    } #rep
  } #date
} #var_name

# Integrate data
col_names <- c("4wC0_0829_1", "4wC0_0829_2", "4wC0_0829_3",
               "4wC2_0912_1", "4wC2_0912_2", "4wC2_0912_3",
               "4wC4_0926_1", "4wC4_0926_2", "4wC4_0926_3",
               "C4W1_1003_1", "C4W1_1003_2", "C4W1_1003_3")

df_H3K4me3 <- cbind(H3K4me3_220829_4wC0_rep1, H3K4me3_220829_4wC0_rep2[,2], H3K4me3_220829_4wC0_rep3[,2],
                    H3K4me3_220912_4wC2_rep1[,2], H3K4me3_220912_4wC2_rep2[,2], H3K4me3_220912_4wC2_rep3[,2],
                    H3K4me3_220926_4wC4_rep1[,2], H3K4me3_220926_4wC4_rep2[,2], H3K4me3_220926_4wC4_rep3[,2],
                    H3K4me3_221003_C4W1_rep1[,2], H3K4me3_221003_C4W1_rep2[,2], H3K4me3_221003_C4W1_rep3[,2])
names(df_H3K4me3) <- c("Ahal_ID", col_names)

df_H3K27me3 <- cbind(H3K27me3_220829_4wC0_rep1, H3K27me3_220829_4wC0_rep2[,2], H3K27me3_220829_4wC0_rep3[,2],
                     H3K27me3_220912_4wC2_rep1[,2], H3K27me3_220912_4wC2_rep2[,2], H3K27me3_220912_4wC2_rep3[,2],
                     H3K27me3_220926_4wC4_rep1[,2], H3K27me3_220926_4wC4_rep2[,2], H3K27me3_220926_4wC4_rep3[,2],
                     H3K27me3_221003_C4W1_rep1[,2], H3K27me3_221003_C4W1_rep2[,2], H3K27me3_221003_C4W1_rep3[,2])
names(df_H3K27me3) <- c("Ahal_ID", col_names)

# Load orthomap data
df_orthomap <- as.data.frame(fread(paste0("01_TAI_250403/", "orthomap_main_rev_250416.tsv")))

# Extract genes with maxRNA > 2
df_spline <- as.data.frame(fread(paste0("07_TAI_seasonal_spline_over4M_250430/", "log2rpkm_sp_1wk.csv"), header = T))
df_spline_maxRNAover2 <- df_spline %>%
  mutate(max_log2rpkm = do.call(pmax, c(across(where(is.numeric)), na.rm = TRUE))) %>%
  filter(max_log2rpkm > 2) %>%
  select(-max_log2rpkm)

df_H3K4me3 <- inner_join(df_H3K4me3, df_spline_maxRNAover2["Ahal_ID"], by = "Ahal_ID")
df_H3K27me3 <- inner_join(df_H3K27me3, df_spline_maxRNAover2["Ahal_ID"], by = "Ahal_ID")





##### H3K4me3 #####

# Merge rpkm and orthomap
df_merge <- inner_join(df_H3K4me3, df_orthomap, by = "Ahal_ID") %>%
  select(-PSname, -PScontinuity) %>%
  relocate(Ahal_ID, PSnum)
df_merge2 <- inner_join(df_H3K4me3, df_orthomap, by = "Ahal_ID") %>%
  select(-PScontinuity) %>%
  relocate(Ahal_ID, PSnum, PSname)
fwrite(df_merge2, file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_K4repmean.csv"), row.names=F)

# Raw data per replicate (maxRNA > 2)
# Calculate TAI
vec_tai <- NULL
for(i in 1:12){
  vec_tai[i] <- sum(df_merge$PSnum*df_merge[,i+2]) / sum(df_merge[,i+2])
}
df_tai <- data.frame(tai = vec_tai, date = c(rep("0 wk",3), rep("2 wk",3), rep("4 wk",3), rep("4 wk + 1 wk",3)))
fwrite(df_tai, file=paste0(out, "K4AI_cold_rawdata_maxRNAover2_perrep.csv"), row.names=F)


### Plotting3, rep MEAN +- SE
se <- function(x){sd(x) / sqrt(length(x))}
df_tai <- as.data.frame(fread(paste0(out, "K4AI_cold_rawdata_maxRNAover2_perrep.csv")))
df_tai$time <- as.POSIXct(c(rep("2022-08-29 12:00:00 JST",3), rep("2022-09-12 12:00:00 JST",3),
                            rep("2022-09-26 12:00:00 JST",3), rep("2022-10-03 12:00:00 JST",3)))
df_tai$days <- difftime(df_tai$time, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")
df_tai2 <- df_tai %>%
  group_by(days) %>%
  summarize(mean_tai = mean(tai),
            mean_minus_se_tai = mean(tai) - se(tai),
            mean_plus_se_tai = mean(tai) + se(tai),
            se_tai = se(tai))

# Tukey-Kramer test
df_test <- data.frame(data = df_tai$tai, group = df_tai$days)
df_test$group <- factor(df_test$group)

res <- aov(data ~ group, data = df_test)
tuk <- glht(res, linfct = mcp(group = "Tukey"))
tuk_cld <- cld(tuk, decreasing = F)
mark <- as.character(tuk_cld[["mcletters"]][["Letters"]])


# Load temp setting
df_tai3 <- df_tai %>% group_by(time) %>% summarise(mean = mean(tai))
df_temp_warm <- as.data.frame(fread("data/room1-2_lower_rev.csv"))
df_temp_cold <- as.data.frame(fread("data/cold_1_left_rev.csv"))
df_temp_warm$time <- force_tz(mdy_hms(df_temp_warm$time), tzone = "Asia/Tokyo")
df_temp_cold$time <- force_tz(mdy_hms(df_temp_cold$time), tzone = "Asia/Tokyo")
df_temp_warm1 <- df_temp_warm %>%
  filter(time >= (df_tai3$time[1] - 60*60*24*7) & time <= df_tai3$time[1])
df_temp_cold1 <- df_temp_cold %>%
  filter(time >= df_tai3$time[1] & time <= df_tai3$time[3])
df_temp_warm2 <- df_temp_warm %>%
  filter(time >= df_tai3$time[3] & time <= df_tai3$time[4] + 60*60*24*7)
df_temp <- rbind(df_temp_warm1, df_temp_cold1, df_temp_warm2)
df_temp$days <- difftime(df_temp$time, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")

# Plotting
local({
  # Setting var_nameiables
  xmax <- max(df_temp$days, na.rm = T)
  xmin <- min(df_temp$days, na.rm = T)
  ymax <- max(df_tai2[,2:4], na.rm = T)
  ymin <- min(df_tai2[,2:4], na.rm = T)
  
  ymean <- (ymax + ymin) / 2
  ymax <- ymean + max(df_tai2[,5], na.rm = T) * 3.3
  ymin <- ymean - max(df_tai2[,5], na.rm = T) * 3.3
  
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_tai <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_tai_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
  df_temp$temp_tai <- scale_to_tai(df_temp$temp)
  
  g_K4ai <<- ggplot() +
    geom_line(data = df_temp, aes(x = days, y = temp_tai), col = "gray50", linewidth = 0.2) + 
    # geom_hline(yintercept = 0, col='gray50', linetype='dashed') +
    # geom_pointrange(data = df_tai, aes(x=days, y=tai, ymin=mean_minus_se, ymax=mean_plus_se), 
    #                 color='grey5', fill='grey95', size=0.3, fatten = 0.2) +
    geom_point(data = df_tai2, aes(x=days, y=mean_tai), size=0.8) +
    geom_errorbar(data = df_tai2, aes(x=days, ymin=mean_minus_se_tai, ymax=mean_plus_se_tai), width = 1.2, linewidth = 0.5) +
    # geom_jitter(data = df_tai, aes(x = days, y = tai), color = "black", alpha = 1, size = 0.5, width = 0.25) +
    annotate(geom = "text", x=df_tai2$days, y=rep(yceiling+yrange*0.07,length(mark)), label=mark, size=6/ggplot2::.pt) +
    annotate(geom = "text", label = "Young", x=-12, y=yceiling+yrange*0.22, size=6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "Old", x=-12, y=yfloor-yrange*0.22, size=6/ggplot2::.pt, vjust=0) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling+yrange*0.15), clip='off') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_tai_rev(.), name = "Temperature (°C)")) +
    theme_test(base_size = 6) +
    theme(plot.title=element_text(size = 6, color = "black", face = "bold"), 
          axis.title=element_text(size = 6, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text=element_text(size = 6, color = "black"),
          plot.margin = unit(c(1, 1, 1, 1),"mm"),
          plot.tag = element_text(size = 8, face = "bold")) +
    labs(title = "H3K4me3",
         x = "Days after transfer to cold",
         y = "Age index")
})



### All PS
df_merge <- as.data.frame(fread(file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_K4repmean.csv")))
names(df_merge)
dim(df_merge)
df_merge[1:5,1:5]

df_mean <- as.data.frame(cbind(apply(df_merge[,4:6], 1, mean), apply(df_merge[,7:9], 1, mean), 
                               apply(df_merge[,10:12], 1, mean), apply(df_merge[,13:15], 1, mean)))
df_merge2 <- cbind(df_merge[,1:3], df_mean)
dates <- c("2022-08-29 12:00:00 JST", "2022-09-12 12:00:00 JST",
           "2022-09-26 12:00:00 JST", "2022-10-03 12:00:00 JST")
names(df_merge2) <- c(names(df_merge)[1:3], dates)

# Calculate mean expression for each PS at each date
df_eachPS <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), mean, .names = "{.col}"))

df_plot <- df_eachPS %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean"
  ) %>%
  mutate(date = as.POSIXct(date))  
df_plot$days <- as.numeric(difftime(df_plot$date, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days"))


# Plotting
xmax <- max(df_temp$days, na.rm = T)
xmin <- min(df_temp$days, na.rm = T)
ymax <- max(df_plot[,3], na.rm = T)
ymin <- min(df_plot[,3], na.rm = T)

ymean <- (ymax + ymin) / 2
ymax <- ymean + 2
ymin <- ymean - 2

yrange <- (ymax - ymin)
yceiling <-  ymax + yrange * 0.05
yfloor <- ymin - yrange * 0.05

# Scale adjustment
scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
df_temp$temp_rna <- scale_to_rna(df_temp$temp)

# Setting labels
g_K4 <- ggplot(df_plot, aes(x = days, y = mean, group = PSnum, color = PSnum)) +
  geom_line() +
  scale_color_viridis(name="PS", option="turbo") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
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
       x = "Days after transfer to cold",
       y = "Mean level")





##### H3K27me3 #####

# Merge rpkm and orthomap
df_merge <- inner_join(df_H3K27me3, df_orthomap, by = "Ahal_ID") %>%
  select(-PSname, -PScontinuity) %>%
  relocate(Ahal_ID, PSnum)
df_merge2 <- inner_join(df_H3K27me3, df_orthomap, by = "Ahal_ID") %>%
  select(-PScontinuity) %>%
  relocate(Ahal_ID, PSnum, PSname)
fwrite(df_merge2, file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_K27repmean.csv"), row.names=F)

# Raw data per replicate (maxRNA > 2)
# Calculate TAI
vec_tai <- NULL
for(i in 1:12){
  vec_tai[i] <- sum(df_merge$PSnum*df_merge[,i+2]) / sum(df_merge[,i+2])
}
df_tai <- data.frame(tai = vec_tai, date = c(rep("0 wk",3), rep("2 wk",3), rep("4 wk",3), rep("4 wk + 1 wk",3)))
fwrite(df_tai, file=paste0(out, "K27AI_cold_rawdata_maxRNAover2_perrep.csv"), row.names=F)


### Plotting3, rep MEAN +- SE
se <- function(x){sd(x) / sqrt(length(x))}
df_tai <- as.data.frame(fread(paste0(out, "K27AI_cold_rawdata_maxRNAover2_perrep.csv")))
df_tai$time <- as.POSIXct(c(rep("2022-08-29 12:00:00 JST",3), rep("2022-09-12 12:00:00 JST",3),
                            rep("2022-09-26 12:00:00 JST",3), rep("2022-10-03 12:00:00 JST",3)))
df_tai$days <- difftime(df_tai$time, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")
df_tai2 <- df_tai %>%
  group_by(days) %>%
  summarize(mean_tai = mean(tai),
            mean_minus_se_tai = mean(tai) - se(tai),
            mean_plus_se_tai = mean(tai) + se(tai),
            se_tai = se(tai))

# Tukey-Kramer test
df_test <- data.frame(data = df_tai$tai, group = df_tai$days)
df_test$group <- factor(df_test$group)

res <- aov(data ~ group, data = df_test)
tuk <- glht(res, linfct = mcp(group = "Tukey"))
tuk_cld <- cld(tuk, decreasing = F)
mark <- as.character(tuk_cld[["mcletters"]][["Letters"]])


# Load temp setting
df_tai3 <- df_tai %>% group_by(time) %>% summarise(mean = mean(tai))
df_temp_warm <- as.data.frame(fread("data/room1-2_lower_rev.csv"))
df_temp_cold <- as.data.frame(fread("data/cold_1_left_rev.csv"))
df_temp_warm$time <- force_tz(mdy_hms(df_temp_warm$time), tzone = "Asia/Tokyo")
df_temp_cold$time <- force_tz(mdy_hms(df_temp_cold$time), tzone = "Asia/Tokyo")
df_temp_warm1 <- df_temp_warm %>%
  filter(time >= (df_tai3$time[1] - 60*60*24*7) & time <= df_tai3$time[1])
df_temp_cold1 <- df_temp_cold %>%
  filter(time >= df_tai3$time[1] & time <= df_tai3$time[3])
df_temp_warm2 <- df_temp_warm %>%
  filter(time >= df_tai3$time[3] & time <= df_tai3$time[4] + 60*60*24*7)
df_temp <- rbind(df_temp_warm1, df_temp_cold1, df_temp_warm2)
df_temp$days <- difftime(df_temp$time, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days")

# Plotting
local({
  # Setting var_nameiables
  xmax <- max(df_temp$days, na.rm = T)
  xmin <- min(df_temp$days, na.rm = T)
  ymax <- max(df_tai2[,2:4], na.rm = T)
  ymin <- min(df_tai2[,2:4], na.rm = T)
  
  ymean <- (ymax + ymin) / 2
  ymax <- ymean + max(df_tai2[,5], na.rm = T) * 3.3
  ymin <- ymean - max(df_tai2[,5], na.rm = T) * 3.3
  
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_tai <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_tai_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
  df_temp$temp_tai <- scale_to_tai(df_temp$temp)
  
  g_K27ai <<- ggplot() +
    geom_line(data = df_temp, aes(x = days, y = temp_tai), col = "gray50", linewidth = 0.2) + 
    # geom_hline(yintercept = 0, col='gray50', linetype='dashed') +
    # geom_pointrange(data = df_tai, aes(x=days, y=tai, ymin=mean_minus_se, ymax=mean_plus_se), 
    #                 color='grey5', fill='grey95', size=0.3, fatten = 0.2) +
    geom_point(data = df_tai2, aes(x=days, y=mean_tai), size=0.8) +
    geom_errorbar(data = df_tai2, aes(x=days, ymin=mean_minus_se_tai, ymax=mean_plus_se_tai), width = 1.2, linewidth = 0.5) +
    # geom_jitter(data = df_tai, aes(x = days, y = tai), color = "black", alpha = 1, size = 0.5, width = 0.25) +
    annotate(geom = "text", x=df_tai2$days, y=rep(yceiling+yrange*0.07,length(mark)), label=mark, size=6/ggplot2::.pt) +
    annotate(geom = "text", label = "Young", x=-12, y=yceiling+yrange*0.22, size=6/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "Old", x=-12, y=yfloor-yrange*0.22, size=6/ggplot2::.pt, vjust=0) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling+yrange*0.15), clip='off') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_tai_rev(.), name = "Temperature (°C)")) +
    theme_test(base_size = 6) +
    theme(plot.title=element_text(size = 6, color = "black", face = "bold"), 
          axis.title=element_text(size = 6, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text=element_text(size = 6, color = "black"),
          plot.margin = unit(c(1, 1, 1, 1),"mm"),
          plot.tag = element_text(size = 8, face = "bold")) +
    labs(title = "H3K27me3",
         x = "Days after transfer to cold",
         y = "Age index")
})



### All PS
df_merge <- as.data.frame(fread(file=paste0(out, "PS_cold_log2rpkm_maxRNAover2_K27repmean.csv")))
names(df_merge)
dim(df_merge)
df_merge[1:5,1:5]

df_mean <- as.data.frame(cbind(apply(df_merge[,4:6], 1, mean), apply(df_merge[,7:9], 1, mean), 
                               apply(df_merge[,10:12], 1, mean), apply(df_merge[,13:15], 1, mean)))
df_merge2 <- cbind(df_merge[,1:3], df_mean)
dates <- c("2022-08-29 12:00:00 JST", "2022-09-12 12:00:00 JST",
           "2022-09-26 12:00:00 JST", "2022-10-03 12:00:00 JST")
names(df_merge2) <- c(names(df_merge)[1:3], dates)

# Calculate mean expression for each PS at each date
df_eachPS <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), mean, .names = "{.col}"))

df_plot <- df_eachPS %>%
  pivot_longer(
    cols = -PSnum,
    names_to = "date",
    values_to = "mean"
  ) %>%
  mutate(date = as.POSIXct(date))  
df_plot$days <- as.numeric(difftime(df_plot$date, as.POSIXct("2022-08-29 12:00:00 JST"), units = "days"))


# Plotting
xmax <- max(df_temp$days, na.rm = T)
xmin <- min(df_temp$days, na.rm = T)
ymax <- max(df_plot[,3], na.rm = T)
ymin <- min(df_plot[,3], na.rm = T)

ymean <- (ymax + ymin) / 2
ymax <- ymean + 2
ymin <- ymean - 2

yrange <- (ymax - ymin)
yceiling <-  ymax + yrange * 0.05
yfloor <- ymin - yrange * 0.05

# Scale adjustment
scale_to_rna <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
scale_to_rna_rev <- function(x) ((x - ymin) * (max(df_temp$temp) - min(df_temp$temp))) / (ymax - ymin) + min(df_temp$temp)
df_temp$temp_rna <- scale_to_rna(df_temp$temp)

# Setting labels
g_K27 <- ggplot(df_plot, aes(x = days, y = mean, group = PSnum, color = PSnum)) +
  geom_line() +
  scale_color_viridis(name="PS", option="turbo") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
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
       x = "Days after transfer to cold",
       y = "Mean level")





### Integration of plots
void <- ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm"))

gall <- {(g_tai + labs(tag = "A")) + (g_K27ai + labs(tag = "C")) + (g_K4ai + labs(tag = "E"))} /
  {(g_RNA + labs(tag = "B")) + (g_K27 + labs(tag = "D")) + (g_K4 + labs(tag = "F"))} +
  plot_layout(heights = c(1, 1.5))

ggsave(paste0(out, "Fig.5_TAI_HAI_cold.pdf"),
       gall, width = 150, height = 70, units = "mm")

