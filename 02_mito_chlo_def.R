
# Load packages
library(tidyverse)
library(patchwork)
library(data.table)
library(viridis)

# Create output directory
out <- "02_mito_chlo_def/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}



### Gene ratio mito
df_ps <- as.data.frame(fread(paste0("16_Mito_genes_251005/csv_each_ps/ratiogenePS_allmitoGO_log2rpkm_maxRNAover2_sp_1wk.csv")))

cols <- viridis::turbo(n = max(df_ps$PS))
cols <- cols[df_ps$PS]
df_ps$PS <- factor(df_ps$PS, levels = df_ps$PS)

g_ratio_mito <- ggplot(data = df_ps, aes(x = PS, y = ratio_gene, fill = PS)) +
  geom_bar(stat = "identity", width = 0.7) +
  # coord_flip(ylim = c(0, ymax)) +
  scale_fill_manual(values = cols) +
  # scale_x_discrete(breaks = PSuniq) +
  scale_y_continuous(labels = scales::comma) +
  theme_test() +
  theme(legend.position="none",
        plot.title = element_text(size = 6, colour = "black", face = "bold", vjust = -1),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm")) +
  labs(title = "Mitochondrion-localized genes",
       y = "Percentage of genes")



### Gene ratio chlo
df_ps <- as.data.frame(fread(paste0("17_Chlo_genes_251005/csv_each_ps/ratiogenePS_allchloGO_log2rpkm_maxRNAover2_sp_1wk.csv")))

cols <- viridis::turbo(n = max(df_ps$PS))
cols <- cols[df_ps$PS]
df_ps$PS <- factor(df_ps$PS, levels = df_ps$PS)

g_ratio_chlo <- ggplot(data = df_ps, aes(x = PS, y = ratio_gene, fill = PS)) +
  geom_bar(stat = "identity", width = 0.7) +
  # coord_flip(ylim = c(0, ymax)) +
  scale_fill_manual(values = cols) +
  # scale_x_discrete(breaks = PSuniq) +
  scale_y_continuous(labels = scales::comma) +
  theme_test() +
  theme(legend.position="none",
        plot.title = element_text(size = 6, colour = "black", face = "bold", vjust = -1),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm")) +
  labs(title = "Chloroplast-localized genes",
       y = "Percentage of genes")



### Gene ratio def
df_ps <- as.data.frame(fread(paste0("14_Defense_genes_250927/csv_each_ps/ratiogenePS_alldefGO_log2rpkm_maxRNAover2_sp_1wk.csv")))

cols <- viridis::turbo(n = max(df_ps$PS))
cols <- cols[df_ps$PS]
df_ps$PS <- factor(df_ps$PS, levels = df_ps$PS)

g_ratio_def <- ggplot(data = df_ps, aes(x = PS, y = ratio_gene, fill = PS)) +
  geom_bar(stat = "identity", width = 0.7) +
  # coord_flip(ylim = c(0, ymax)) +
  scale_fill_manual(values = cols) +
  # scale_x_discrete(breaks = PSuniq) +
  scale_y_continuous(labels = scales::comma) +
  theme_test() +
  theme(legend.position="none",
        plot.title = element_text(size = 6, colour = "black", face = "bold", vjust = -1),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(color = "black"),
        plot.tag = element_text(size = 8, face = "bold"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm")) +
  labs(title = "Defense-related genes",
       y = "Percentage of genes")



### Expression of mito GO
df_mito <- as.data.frame(fread(paste0("16_Mito_genes_251005/mitoGO_gene_num.csv")))
mito_go <- df_mito$GOname
first_upper <- function(cha){
  h <- substr(cha,1,1)
  t <- substr(cha,2,nchar(cha))
  out <- paste0(toupper(h),t)
  return(out)
}

# Merge rpkm and orthomap
glist_mito <- list()
for(j in 1:length(mito_go)){
  
  ### Plotting
  df_stat <- as.data.frame(fread(paste0("16_Mito_genes_251005/csv_each/", mito_go[j], "_log2rpkm_maxRNAover2_sp_1wk_stat.csv")))
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
    if(sum(is.na(df_stat$se))==nrow(df_stat)){
      ymax <- max(df_stat$mean, na.rm = T)
      ymin <- min(df_stat$mean, na.rm = T)
    }else{
      ymax <- max(df_stat$mean+df_stat$se, na.rm = T)
      ymin <- min(df_stat$mean-df_stat$se, na.rm = T)
    }
    yrange <- (ymax - ymin)
    yceiling <-  ymax + yrange * 0.05
    yfloor <- ymin - yrange * 0.05
    
    # Scale adjustment
    scale_to_stat <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
    scale_to_stat_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
    df_temp2$temp_stat <- scale_to_stat(df_temp2$dailytemp)
    
    glist_mito[[j]] <<- ggplot() +
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
            plot.margin = unit(c(0.5,0.5,1,0.5), "mm"),
            plot.tag = element_text(size = 8, face = "bold")) +
      labs(title = first_upper(mito_go[j]), x = "", y = "Mean expression") +
      scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                         expand = c(0, 0),
                         sec.axis = sec_axis(~ scale_to_stat_rev(.), name = "Temperature (°C)")) +
      annotate(geom = "text", label = "2011", 
               x=as.POSIXct("2011-09-25"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2012", 
               x=as.POSIXct("2012-07-01"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2013", 
               x=as.POSIXct("2013-04-10"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2011-6-26"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2011-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2012-01-07"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2012-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2013-01-07"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2013-07-06"), yend=yfloor-yrange*0.15, linewidth = 0.2)
  })
}



### Expression of chlo GO
df_chlo <- as.data.frame(fread(paste0("17_Chlo_genes_251005/chloGO_gene_num.csv")))
chlo_go <- df_chlo$GOname
first_upper <- function(cha){
  h <- substr(cha,1,1)
  t <- substr(cha,2,nchar(cha))
  out <- paste0(toupper(h),t)
  return(out)
}

# Merge rpkm and orthomap
glist_chlo <- list()
for(j in 1:length(chlo_go)){
  
  ### Plotting
  df_stat <- as.data.frame(fread(paste0("17_Chlo_genes_251005/csv_each/", chlo_go[j], "_log2rpkm_maxRNAover2_sp_1wk_stat.csv")))
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
    if(sum(is.na(df_stat$se))==nrow(df_stat)){
      ymax <- max(df_stat$mean, na.rm = T)
      ymin <- min(df_stat$mean, na.rm = T)
    }else{
      ymax <- max(df_stat$mean+df_stat$se, na.rm = T)
      ymin <- min(df_stat$mean-df_stat$se, na.rm = T)
    }
    yrange <- (ymax - ymin)
    yceiling <-  ymax + yrange * 0.05
    yfloor <- ymin - yrange * 0.05
    
    # Scale adjustment
    scale_to_stat <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
    scale_to_stat_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
    df_temp2$temp_stat <- scale_to_stat(df_temp2$dailytemp)
    
    glist_chlo[[j]] <<- ggplot() +
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
            plot.margin = unit(c(0.5,0.5,1,0.5), "mm"),
            plot.tag = element_text(size = 8, face = "bold")) +
      labs(title = first_upper(chlo_go[j]), x = "", y = "Mean expression") +
      scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                         expand = c(0, 0),
                         sec.axis = sec_axis(~ scale_to_stat_rev(.), name = "Temperature (°C)")) +
      annotate(geom = "text", label = "2011", 
               x=as.POSIXct("2011-09-25"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2012", 
               x=as.POSIXct("2012-07-01"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2013", 
               x=as.POSIXct("2013-04-10"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2011-6-26"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2011-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2012-01-07"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2012-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2013-01-07"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2013-07-06"), yend=yfloor-yrange*0.15, linewidth = 0.2)
  })
}



### Expression of def GO
df_def <- as.data.frame(fread(paste0("14_Defense_genes_250927/defenseGO_gene_num.csv")))
def_go <- df_def$GOname
first_upper <- function(cha){
  h <- substr(cha,1,1)
  t <- substr(cha,2,nchar(cha))
  out <- paste0(toupper(h),t)
  return(out)
}

# Merge rpkm and orthomap
glist_def <- list()
for(j in 1:length(def_go)){
  
  ### Plotting
  if(!file.exists(paste0("14_Defense_genes_250927/csv_each/", def_go[j], "_log2rpkm_maxRNAover2_sp_1wk_stat.csv"))){
    next
  }
  df_stat <- as.data.frame(fread(paste0("14_Defense_genes_250927/csv_each/", def_go[j], "_log2rpkm_maxRNAover2_sp_1wk_stat.csv")))
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
    if(sum(is.na(df_stat$se))==nrow(df_stat)){
      ymax <- max(df_stat$mean, na.rm = T)
      ymin <- min(df_stat$mean, na.rm = T)
    }else{
      ymax <- max(df_stat$mean+df_stat$se, na.rm = T)
      ymin <- min(df_stat$mean-df_stat$se, na.rm = T)
    }
    yrange <- (ymax - ymin)
    yceiling <-  ymax + yrange * 0.05
    yfloor <- ymin - yrange * 0.05
    
    # Scale adjustment
    scale_to_stat <- function(x) ((x - min(x)) / (max(x) - min(x))) * (ymax - ymin) + ymin
    scale_to_stat_rev <- function(x) ((x - ymin) * (max(df_temp2$dailytemp) - min(df_temp2$dailytemp))) / (ymax - ymin) + min(df_temp2$dailytemp)
    df_temp2$temp_stat <- scale_to_stat(df_temp2$dailytemp)
    
    glist_def[[j]] <<- ggplot() +
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
            plot.margin = unit(c(0.5,0.5,1,0.5), "mm"),
            plot.tag = element_text(size = 8, face = "bold")) +
      labs(title = first_upper(def_go[j]), x = "", y = "Mean expression") +
      scale_y_continuous(labels = function(x) format(x, nsmall = 1),
                         expand = c(0, 0),
                         sec.axis = sec_axis(~ scale_to_stat_rev(.), name = "Temperature (°C)")) +
      annotate(geom = "text", label = "2011", 
               x=as.POSIXct("2011-09-25"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2012", 
               x=as.POSIXct("2012-07-01"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      annotate(geom = "text", label = "2013", 
               x=as.POSIXct("2013-04-10"), y=yfloor-yrange*0.25,
               size=6/ggplot2::.pt, vjust=0) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2011-6-26"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2011-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2012-01-07"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2012-12-25"), yend=yfloor-yrange*0.15, linewidth = 0.2) +
      geom_segment(data = df_stat, mapping = aes(x = date, y = stat),
                   x=as.POSIXct("2013-01-07"), y=yfloor-yrange*0.15, 
                   xend=as.POSIXct("2013-07-06"), yend=yfloor-yrange*0.15, linewidth = 0.2)
  })
}





### Integration of plots
void <- ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm"))
# gall <- (g_ratio_mito + labs(tag = "a")) +
#   (g_ratio_chlo + labs(tag = "b")) +
#   (g_ratio_def + labs(tag = "c")) +
#   glist_mito[[2]] + labs(tag = "d") + theme(axis.title.y.right = element_blank()) + 
#   glist_mito[[3]] + theme(axis.title = element_blank()) + 
#   glist_mito[[4]] + theme(axis.title.y.left = element_blank()) +
#   glist_chlo[[2]] + labs(tag = "e") + theme(axis.title.y.right = element_blank()) +
#   glist_chlo[[3]] + theme(axis.title = element_blank()) + 
#   glist_chlo[[4]] + theme(axis.title.y.left = element_blank()) +
#   glist_def[[1]] + labs(tag = "f") + theme(axis.title.y.right = element_blank()) + 
#   glist_def[[3]] + theme(axis.title = element_blank()) + 
#   glist_def[[4]] + theme(axis.title.y.left = element_blank()) +
#   plot_layout(ncol = 3)

gall <- {(g_ratio_mito + labs(tag = "A")) +
    (g_ratio_chlo + labs(tag = "B")) +
    (g_ratio_def + labs(tag = "C")) +
    plot_layout(ncol = 3)} /
  {
  {glist_mito[[2]] + labs(tag = "D") + 
      glist_mito[[3]] + 
      glist_mito[[4]] + 
      plot_layout(ncol = 1)} |
  {glist_chlo[[2]] + labs(tag = "E") +
      glist_chlo[[3]] +
      glist_chlo[[4]] +
      plot_layout(ncol = 1)} |
  {glist_def[[1]] + labs(tag = "F") +
      glist_def[[3]] + 
      glist_def[[4]] +
      plot_layout(ncol = 1)}
  } +
  plot_layout(heights = c(1,4))

ggsave(paste0(out, "02_mito_chlo_def.pdf"),
       gall, width = 180, height = 140, units = "mm")

