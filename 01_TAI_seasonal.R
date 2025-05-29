
# Load packages
library(tidyverse)
library(patchwork)
library(data.table)

# Create the output directory
out <- "01_TAI_seasonal/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}



##### Seasonal TAI #####
# Load data from https://sohi.ecology.kyoto-u.ac.jp/AhgRNAseq/Nishio_Nat.Plants_script_data.zip
# and place 20181120_bowtie2_rpm.Rdata and 140519_SampleAttribute.txt in data/ folder
load("data/20181120_bowtie2_rpm.Rdata")
dim(rpkm)

rawcnt_s <- rawcnt[8410:nrow(rpkm), 1:490]
vec_mappedreads <- apply(rawcnt_s, 2, sum)
hist(vec_mappedreads, breaks = 48)
idx_4M <- as.numeric(which(vec_mappedreads > 4*10^6))

log2rpkm_s <- log2(rpkm[8410:nrow(rpkm), 1:490] + 1)
log2rpkm_s <- log2rpkm_s[-24365, idx_4M]
dim(log2rpkm_s)

attribute <- read.table("data/140519_SampleAttribute.txt", header=T,sep="\t")
dim(attribute)
attribute_s <- attribute[1:490, 1:5]
attribute_s <- attribute_s[idx_4M,]
dates <- as.POSIXct(paste(attribute_s[,'year'], 
                          attribute_s[,'month'], 
                          attribute_s[,'day'],
                          attribute_s[,'hour'],
                          sep='/'), 
                    format='%Y/%m/%d/%H')
days <- as.numeric(dates - as.POSIXct("2011-06-30 12:00:00"))
days2 <- c(days,730+days,1460+days)

# Change colnames of log2rpkm_s
first <- str_split_fixed(colnames(log2rpkm_s), "_", 2)[,1]
second <- as.character(as.numeric(str_split_fixed(colnames(log2rpkm_s), "_", 2)[,2]))
colnames_rev <- paste(first, second, sep = "-")
colnames(log2rpkm_s) <- colnames_rev

# Prepare a data frame
df_spline <- as.data.frame(log2rpkm_s)
Ahal_ID <- str_split_fixed(row.names(log2rpkm_s), "\\.", 2)[,1]
df_spline <- df_spline %>%
  mutate(Ahal_ID = Ahal_ID) %>%
  relocate(Ahal_ID)





##### Spline #####
log2rpkm_s2 <- cbind(log2rpkm_s, log2rpkm_s, log2rpkm_s)
dim(log2rpkm_s2)

# Container of results
start_time <- as.POSIXct("2011-07-07 12:00:00", tz = "Asia/Tokyo")
end_time <- start_time + lubridate::years(2)
time_seq <- seq(from = start_time, to = end_time, by = "7 days")
mat_spline <- matrix(nrow = nrow(log2rpkm_s2), ncol = length(time_seq))

for(i in 1:nrow(log2rpkm_s2)){
  spRNA <- smooth.spline(days2, log2rpkm_s2[i,], spar=0.3)
  length <- 10000
  x <- seq(730, 1470, length=length)
  predRNA <- predict(spRNA, x)
  
  # 1 wk interval
  # 0, 730, 1460 days: "2011-06-30 12:00:00"
  for(j in 1:length(time_seq)){
    mat_spline[i,j] <- predRNA$y[which.min(abs(predRNA$x - (730 + 7*j)))]
  }
}

# Save result
df_spline <- as.data.frame(mat_spline)
names(df_spline) <- time_seq
df_spline2 <- cbind(data.frame(Ahal_ID = Ahal_ID),
                    df_spline)
fwrite(df_spline2, paste0(out, "log2rpkm_sp_1wk.csv"))





##### Raw data (rep mean, maxRNA > 2) #####

# Load orthomap data
df_orthomap <- as.data.frame(fread(file = paste0("data/", "orthomap_main_rev_250416.tsv")))
df_spline <- as.data.frame(fread(paste0(out, "log2rpkm_sp_1wk.csv"), header = T))

# Extract genes with maxRNA > 2
df_spline_maxRNAover2 <- df_spline %>%
  mutate(max_log2rpkm = do.call(pmax, c(across(where(is.numeric)), na.rm = TRUE))) %>%
  filter(max_log2rpkm > 2) %>%
  select(-max_log2rpkm)

# Merge rpkm and orthomap
df_merge <- inner_join(df_spline_maxRNAover2, df_orthomap, by = "Ahal_ID") %>%
  select(-PSname, -PScontinuity) %>%
  relocate(Ahal_ID, PSnum)
df_merge2 <- inner_join(df_spline_maxRNAover2, df_orthomap, by = "Ahal_ID") %>%
  select(-PScontinuity) %>%
  relocate(Ahal_ID, PSnum, PSname)
fwrite(df_merge2, file=paste0(out, "PS_log2rpkm_maxRNAover2_sp_1wk.csv"), row.names=F)

# Calculate TAI
vec_tai <- NULL
for(i in 1:(ncol(df_merge) - 2)){
  vec_tai[i] <- sum(df_merge$PSnum*df_merge[,i+2]) / sum(df_merge[,i+2])
}
df_tai <- data.frame(tai = vec_tai, date = names(df_merge)[-(1:2)])

# bootstrapping
n_iter <- 10000
bootstrap <- matrix(NA, nrow=n_iter, ncol=ncol(df_merge) - 2)
set.seed(42)
for(i in 1:n_iter){
  sampleNum <- sample(nrow(df_merge), nrow(df_merge), replace=T)
  df_merge_sub <- df_merge[sampleNum,]
  for(j in 1:(ncol(df_merge) - 2)){
    bootstrap[i,j] <- sum(df_merge_sub$PSnum*df_merge_sub[,j+2]) / sum(df_merge_sub[,j+2])
  }
}

boot.quant <- matrix(NA,nrow=2,ncol=ncol(df_merge) - 2)
for(i in 1:(ncol(df_merge) - 2)){
  boot.quant[,i] <- quantile(bootstrap[,i],p=c(0.025,0.975))
}

boot.sd <- matrix(NA,nrow=2,ncol=ncol(df_merge) - 2)
for(i in 1:(ncol(df_merge) - 2)){
  mean_b <- mean(bootstrap[,i])
  sd_b <- sd(bootstrap[,i])
  boot.sd[,i] <- c(mean_b - sd_b, mean_b + sd_b)
}

se <- function(x){sd(x) / sqrt(length(x))}
boot.se <- matrix(NA,nrow=2,ncol=ncol(df_merge) - 2)
for(i in 1:(ncol(df_merge) - 2)){
  mean_b <- mean(bootstrap[,i])
  se_b <- se(bootstrap[,i])
  boot.se[,i] <- c(mean_b - se_b, mean_b + se_b)
}

# Save output
df_tai <- data.frame(tai = vec_tai, 
                     q2.5 = boot.quant[1,], 
                     q97.5 = boot.quant[2,], 
                     mean_minus_sd = boot.sd[1,], 
                     mean_plus_sd = boot.sd[2,], 
                     mean_minus_se = boot.se[1,], 
                     mean_plus_se = boot.se[2,], 
                     date = names(df_merge)[-(1:2)])
fwrite(df_tai, file=paste0(out, "TAI_rawdata_maxRNAover2_sp_1wk.csv"), row.names=F)



### Permutation test of TAI
df_tai <- as.data.frame(fread(paste0(out, "TAI_rawdata_maxRNAover2_sp_1wk.csv")))
obs_var <- var(df_tai$tai)

# permutation
n_iter <- 10000
permuted_TAI <- matrix(NA, nrow=n_iter, ncol=(ncol(df_merge) - 2))
set.seed(42)
for(i in 1:n_iter){
  permuted_PSnum <- sample(df_merge$PSnum)
  for(j in 1:(ncol(df_merge) - 2)){
    permuted_TAI[i,j] <- sum(permuted_PSnum*df_merge[,j+2]) / sum(df_merge[,j+2])
  }
}

# Calculate variance and p-values
permuted_var <- apply(permuted_TAI, 1, var)
p_val <- sum(permuted_var >= obs_var) / length(permuted_var)
mean(permuted_var > obs_var)

# Save results
df_permuted_TAI <- as.data.frame(permuted_TAI)
names(df_permuted_TAI) <- df_tai$date
fwrite(df_permuted_TAI, file=paste0(out, "permutedTAI_rawdata_maxRNAover2_sp_1wk.csv"), row.names=F)
df_p <- data.frame(p_value = p_val)
fwrite(df_p, file=paste0(out, "permutionP_rawdata_maxRNAover2_sp_1wk.csv"), row.names=F)





#### Plotting (daily temp)
df_tai <- as.data.frame(fread(paste0(out, "TAI_rawdata_maxRNAover2_sp_1wk.csv")))
df_tai$date <- force_tz(df_tai$date, tzone = "Asia/Tokyo")
obs_var <- var(df_tai$tai)

df_permuted_TAI <- as.data.frame(fread(file=paste0(out, "permutedTAI_rawdata_maxRNAover2_sp_1wk.csv"), header = T))
permuted_var <- apply(df_permuted_TAI, 1, var)
p_val <- sum(permuted_var >= obs_var) / length(permuted_var)
if(p_val == 0){
  p_val <- 1 / length(permuted_var)
}

# Load Nishiwaki dailytemp
df_temp <- as.data.frame(fread("data/Nishiwaki_dailytemp2011-2013.tsv"))
df_temp$date <- as.POSIXct(paste0(df_temp$year, "-", df_temp$month, "-", df_temp$day, " 12:00 JST"))
df_temp2 <- df_temp %>%
  filter(date >= df_tai$date[1] & date <= df_tai$date[nrow(df_tai)])

# Setting var_nameiables
labelx<-c("","7","","9","","11","",
          "1","","3","","5","","7","","9","","11","",
          "1","","3","","5","","7","")
length(labelx)

local({
  xmax <- max(df_tai$date, na.rm = T)
  xmin <- min(df_tai$date, na.rm = T)
  ymax <- max(df_tai[,2:3], na.rm = T)
  ymin <- min(df_tai[,2:3], na.rm = T)
  yrange <- (ymax - ymin)
  yceiling <-  ymax + yrange * 0.05
  yfloor <- ymin - yrange * 0.05
  
  # Scale adjustment
  scale_to_tai <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
  scale_to_tai_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
  df_temp2$temp_tai <- scale_to_tai(df_temp2$dailytemp)
  
  g_ts <<- ggplot() +
    geom_line(data = df_temp2, aes(x = date, y = temp_tai), col = "gray50", linewidth = 0.2) + 
    geom_ribbon(data = df_tai, aes(x = date, ymin = q2.5, ymax = q97.5),
                alpha = 0.3, fill = "gray30") +
    geom_line(data = df_tai, aes(x = date, y = tai), col = "black", linewidth = 0.6) +
    # geom_point(data = df_tai2, aes(x = date, y = tai), col = "black", alpha = 0.6, size = 0.2) +
    # scale_x_continuous(breaks = 1:12) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    annotate("text", 
             x=as.POSIXct("2011-12-01"), y=yceiling-yrange*0.08,
             label = bquote(italic(p) < .(format(p_val, scientific = FALSE))), 
             size=7/ggplot2::.pt) +
    theme_test(base_size=7) +
    theme(plot.title = element_blank(),
          axis.title=element_text(size=7, colour = "black"), 
          axis.text=element_text(size=7, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(1.5,0.2,1,0.2), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = "", x = "", y = "Transcriptome age index") +
    scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                       expand = c(0, 0),
                       sec.axis = sec_axis(~ scale_to_tai_rev(.), name = "Temperature (°C)")) +
    annotate(geom = "text", label = "Young", 
             x=as.POSIXct("2011-05-15"), y=yceiling+yrange*0.05,
             size=7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "Old", 
             x=as.POSIXct("2011-05-15"), y=yfloor-yrange*0.05,
             size=7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2011", 
             x=as.POSIXct("2011-09-25"), y=yfloor-yrange*0.25,
             size=7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2012", 
             x=as.POSIXct("2012-07-01"), y=yfloor-yrange*0.25,
             size=7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013", 
             x=as.POSIXct("2013-04-10"), y=yfloor-yrange*0.25,
             size=7/ggplot2::.pt, vjust=0) +
    geom_segment(data = df_tai, mapping = aes(x = date, y = tai),
                 x=as.POSIXct("2011-6-26"), y=yfloor-yrange*0.15, 
                 xend=as.POSIXct("2011-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
    geom_segment(data = df_tai, mapping = aes(x = date, y = tai),
                 x=as.POSIXct("2012-01-07"), y=yfloor-yrange*0.15, 
                 xend=as.POSIXct("2012-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
    geom_segment(data = df_tai, mapping = aes(x = date, y = tai),
                 x=as.POSIXct("2013-01-07"), y=yfloor-yrange*0.15, 
                 xend=as.POSIXct("2013-07-06"), yend=yfloor-yrange*0.15, linewidth = 0.2)
})


### Correlation between TAI and temperature
df_tai_temp2 <- inner_join(df_tai, df_temp2, by = "date")

# Linear regression
# model <- RobustLinearReg::siegel_regression(tai ~ dailytemp, data = df_tai_temp2)
model <- stats::lm(tai ~ dailytemp, data = df_tai_temp2)
suppressWarnings(conf_interval <- stats::predict(model, interval="confidence", level = 0.95))
conf_interval2 <- as.data.frame(cbind(df_tai_temp2$dailytemp, conf_interval)[order(df_tai_temp2$dailytemp, decreasing = F),])
names(conf_interval2)[1] <- "dailytemp"

# Labels
cor_res <- cor.test(df_tai_temp2$tai, df_tai_temp2$dailytemp, method = "pearson")
mark <- ifelse(cor_res$p.value < 0.0001, "****", 
               ifelse(cor_res$p.value < 0.001, "***", 
                      ifelse(cor_res$p.value < 0.01, "**",
                             ifelse(cor_res$p.value < 0.05, "*", "ns"))))
r <- format(cor_res$estimate, digits=2, nsmall = 2)
rlab <- bquote(paste(italic(r), " = ", .(r), .(mark), sep=""))
intercept <- format(summary(model)$coefficients[1,1], digits=2, nsmall = 2)
coefficient <- format(summary(model)$coefficients[2,1], digits=2, nsmall = 2)
equationlab <- bquote(paste(italic(y), " = ", .(intercept), " + ", .(coefficient), " ", italic(x), sep=""))
if(as.numeric(coefficient) < 0){
  coefficient <- format(-summary(model)$coefficients[2,1], digits=2, nsmall = 2)
  equationlab <- bquote(paste(italic(y), " = ", .(intercept), " - ", .(coefficient), " ", italic(x), sep=""))
}

# Plotting
g_cor <- ggplot() +
  geom_ribbon(data = conf_interval2, aes(x=dailytemp, ymin = lwr, ymax = upr), alpha = 0.4, fill = "steelblue") +
  geom_line(data = conf_interval2, aes(x=dailytemp, y=fit), color = "steelblue", linewidth = 1) +
  geom_point(data = df_tai_temp2, aes(x = dailytemp, y = tai), color = "black", alpha = 1, size = 0.5) +
  # coord_cartesian(xlim=c(xmin, xmax), clip='on') +
  # scale_x_continuous(breaks = seq(0,1,0.2)) +
  # scale_y_continuous(breaks = seq(0,1,0.2)) +
  annotate("text", x=-2, y=8.9,
           label=rlab, size=7/ggplot2::.pt, hjust = 0) +
  theme_test(base_size=7) +
  theme(plot.title=element_blank(),
        axis.title=element_text(size=7, colour = "black"), 
        axis.text=element_text(size=7, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1,2,1,1), "mm"),
        plot.tag = element_text(size = 12)) + 
  labs(x = expression(paste("Temperature (°C)", sep="")),
       y = expression(paste("Transcriptome age index", sep="")))



### Simple Moving Average
temp_ave <- matrix(ncol=200, nrow=nrow(df_tai))
samplepoint <- difftime(df_tai$date, as.POSIXct("2011-01-01 12:00:00 JST"), units = "days")
samplepoint2 <- samplepoint# - 1
length(samplepoint2)

for(i in 1:200){
  mvtemp <- TTR::SMA(df_temp$dailytemp, i)
  for(j in 1:nrow(df_tai)){
    temp_ave[j,i] <- mvtemp[samplepoint2[j]]
  }
}
df_sma <- as.data.frame(temp_ave)
names(df_sma) <- paste0("SMA", 1:200)
df_sma <- cbind(data.frame(date = df_tai$date),
                df_sma)
write.csv(df_sma, paste0(out, "Nishiwaki_dailytemp2011-2013_SMA_sp_1wk.csv"),
          quote=F, row.names=F)

# Calculate correlation
df_sma <- as.data.frame(fread(paste0(out, "Nishiwaki_dailytemp2011-2013_SMA_sp_1wk.csv")))
df_sma$date <- force_tz(df_sma$date, tzone = "Asia/Tokyo")
df_tai_sma <- full_join(df_tai, df_sma, by = "date") %>%
  arrange(date)

df_cors <- data.frame(sma = 1:200,
                      cors = NA)
for(i in 1:nrow(df_cors)){
  df_cors$cors[i] <- cor(df_tai_sma$tai, df_tai_sma[,i+8], use = "pairwise.complete.obs")
}

max_sma <- df_cors$sma[which.max(df_cors$cors)]

# Plotting
g_cors <- ggplot() +
  geom_line(data = df_cors, aes(x=sma, y=cors), color = "black", linewidth = 0.5) +
  geom_vline(xintercept = max_sma, col='gray50', linetype='dashed') +
  # coord_cartesian(xlim=c(xmin, xmax), clip='on') +
  # scale_x_continuous(breaks = seq(0,1,0.2)) +
  # scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme_test(base_size=7) +
  theme(plot.title=element_blank(),
        axis.title=element_text(size=7, colour = "black"), 
        axis.text=element_text(size=7, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1,2,1,1), "mm"),
        plot.tag = element_text(size = 12)) + 
  labs(x = "SMA period (days)",
       y = "Correlation coefficient")





##### Each PS mean expression #####

### Comparison between PSs
df_merge2 <- as.data.frame(fread(file=paste0(out, "PS_log2rpkm_maxRNAover2_sp_1wk.csv")))
df_merge2[1:5,1:5]

# Calculate mean expression for each PS at each date
df_eachPS <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), mean, .names = "{.col}"))

# Calculate SE expression for each PS at each date
se <- function(x){sd(x) / sqrt(length(x))}
df_eachPS_se <- df_merge2 %>%
  group_by(PSnum) %>%
  summarize(across(where(is.numeric), se, .names = "{.col}"))
df_eachPS_minusse <- cbind(df_eachPS_se[,1], (df_eachPS[,2:ncol(df_eachPS)] - df_eachPS_se[,2:ncol(df_eachPS_se)]))
df_eachPS_plusse <- cbind(df_eachPS_se[,1], (df_eachPS[,2:ncol(df_eachPS)] + df_eachPS_se[,2:ncol(df_eachPS_se)]))

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
# Correlation
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
  labelx<-c("","7","","9","","11","",
            "1","","3","","5","","7","","9","","11","",
            "1","","3","","5","","7","")
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
    theme_test(base_size = 7) +
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.title=element_text(size=7, colour = "black"), 
          axis.title.y.right = element_blank(), 
          axis.text=element_text(size=7, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "",
         y = "Mean expression") +
    annotate("text", label=rlab,
             x=as.POSIXct("2012-08-15"), y=yceiling-yrange*0.08, 
             size=7/ggplot2::.pt, hjust = 0) +
    annotate(geom = "text", label = "2011", 
             x=as.POSIXct("2011-09-25"), y=yfloor-yrange*0.34,
             size = 7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2012", 
             x=as.POSIXct("2012-07-01"), y=yfloor-yrange*0.34,
             size = 7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013", 
             x=as.POSIXct("2013-04-10"), y=yfloor-yrange*0.34,
             size = 7/ggplot2::.pt, vjust=0) +
    geom_segment(data = df_plot_rev, mapping = aes(x = date, y = mean),
                 x=as.POSIXct("2011-6-26"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2011-12-25"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(data = df_plot_rev, mapping = aes(x = date, y = mean),
                 x=as.POSIXct("2012-01-07"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2012-12-25"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(data = df_plot_rev, mapping = aes(x = date, y = mean),
                 x=as.POSIXct("2013-01-07"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2013-07-06"), yend=yfloor-yrange*0.2, linewidth = 0.2)
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
  labelx<-c("","7","","9","","11","",
            "1","","3","","5","","7","","9","","11","",
            "1","","3","","5","","7","")
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
                       sec.axis = sec_axis(~ scale_to_rna_rev(.), name = "Temperature (°C)")) +
    scale_x_datetime(expand = c(0.03,0.03), date_breaks="1 month", labels = labelx) +
    coord_cartesian(xlim = c(xmin, xmax), ylim=c(yfloor, yceiling), clip='off') +
    theme_test(base_size = 7) +
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.title=element_text(size=7, colour = "black"), 
          axis.title.y.left=element_blank(), 
          axis.text=element_text(size=7, colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.margin = unit(c(1.5,0.5,1,0.5), "mm"),
          plot.tag = element_text(size = 12)) +
    labs(title = paste0("PS", i, ": ", PSname),
         x = "",
         y = "") +
    annotate("text", label=rlab,
             x=as.POSIXct("2012-09-01"), y=yceiling-yrange*0.08, 
             size=7/ggplot2::.pt, hjust = 0) +
    annotate(geom = "text", label = "2011", 
             x=as.POSIXct("2011-09-25"), y=yfloor-yrange*0.34,
             size = 7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2012", 
             x=as.POSIXct("2012-07-01"), y=yfloor-yrange*0.34,
             size = 7/ggplot2::.pt, vjust=0) +
    annotate(geom = "text", label = "2013", 
             x=as.POSIXct("2013-04-10"), y=yfloor-yrange*0.34,
             size = 7/ggplot2::.pt, vjust=0) +
    geom_segment(data = df_plot_rev, mapping = aes(x = date, y = mean),
                 x=as.POSIXct("2011-6-26"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2011-12-25"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(data = df_plot_rev, mapping = aes(x = date, y = mean),
                 x=as.POSIXct("2012-01-07"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2012-12-25"), yend=yfloor-yrange*0.2, linewidth = 0.2) +
    geom_segment(data = df_plot_rev, mapping = aes(x = date, y = mean),
                 x=as.POSIXct("2013-01-07"), y=yfloor-yrange*0.2, 
                 xend=as.POSIXct("2013-07-06"), yend=yfloor-yrange*0.2, linewidth = 0.2)
})



### Integration of plots
void <- ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm"))
gall <- (g_ts + labs(tag = "A")) / 
  {(g_cor + labs(tag = "B")) + (g_cors + labs(tag = "C"))} /
  {(glist[[1]] + labs(tag = "D")) + glist[[18]]} +
  plot_layout(heights = c(1, 1, 0.8))

ggsave(paste0(out, "Fig.1_TAI_seasonal.pdf"),
       gall, width = 100, height = 125, units = "mm")

