library("dplyr")
library("readxl")
library("gdata")
library(gsubfn)
library(tidyverse)
library("fs")
library(reshape2)
library(ggplot2)
scaleFUN <- function(x) sprintf("%.3f", x)
colours <- c("red", "green", "blue","darkorange","yellow","hotpink")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

get.CRfitness <- function(x) {
  fit.error <- data.frame()
  opt.cr <- vector()
  output.file <- read.table(x,header = TRUE,fill = TRUE,stringsAsFactors = FALSE)
  fname <- tools::file_path_sans_ext(basename(x))
  element.name <- substr(fname, 1, nchar(fname) - 14)
  # index of last 11 rows
  #*******************************************
  only.error <- data.frame()
  fit.ind <- 9
  iter <- ceiling(dim(output.file)[1]/11)
  for(i in 1:iter){
    only.error[i,1:dim(output.file)[2]] <- as.numeric(output.file[fit.ind,])
    fit.ind = fit.ind + 11
  }
  min.fit.ind <- apply(only.error, 1, which.min)[iter]
  apply(only.error, 1, min)[iter]
  len <- dim(output.file)[1]
  CR.err <- output.file[(len-10):(len-1),min.fit.ind]
  #*******************************************
  # rows.index <- dim(output.file)[1] - 10
  # ord.ind <- order(output.file[nrow(output.file) - 1, ],decreasing = TRUE)
  # sort.output <- output.file[, ord.ind]
  # opt.cr <- sort.output[rows.index:(rows.index + 8), 1]
  # min.fitness <- round(as.double(sort.output[nrow(sort.output) - 1, ])[1],digits = 6)
  #fit.error <- as.data.frame(unlist(list(min.fitness, opt.cr)))
  fit.error <- as.data.frame(unlist(CR.err))
  colnames(fit.error) <- element.name
  return(fit.error)
}
#-------------------------------------------------------------------------------------------
get.clusterSize <- function(x) {
  csv.f <-read.csv(x,header = TRUE,fill = TRUE, stringsAsFactors = FALSE)
  element.name <- tools::file_path_sans_ext(basename(x))
  data.point <- dim(csv.f)[1]
  return(list(element.name, data.point))
}

#-----------------------------------------------------------------------------------------------
#***************************** Get fitness and CR of GA Results ********************************
#-----------------------------------------------------------------------------------------------
#ga.path <-"C:/Users/sara/OneDrive - University of Tulsa/MyPhD/Results/Normal/Brito_Felizola/Res_KM"
res.path = ga.path
ga.path = "C:/Users/sara/OneDrive - University of Tulsa/MyPhD/Results/Dimensionless/Brito_Felizola/LM/Res_HCPC"
#ga.path = "C:/Users/sara/OneDrive - University of Tulsa/MyPhD/Results/Normal/WholeData/Res_HCPC"
# ga.path = "C:/Users/sara/Desktop"
#BF_Author_Data <- readRDS("C:/Users/sara/Desktop/phd/EntireData/BritoFelizola/BF_Author_Data.RDS")
#write.csv(BF_Author_Data,paste(ga.path, "Brito_Felizola.csv",sep = "/"))
create.error.result <- function(res.path){
  setwd(res.path)
  extension <- "txt"
  #ga.res <- Sys.glob(paste("*.", extension, sep = ""))[-1]
  ga.res <- Sys.glob(paste("*.", extension, sep = ""))
  length(ga.res)
  fitness.cr <- data.frame()
  k=1
  for(fn in ga.res){
    text.p <- paste(res.path, fn, sep = "/")
    fitness.cr[1:10,k] <- get.CRfitness(x = text.p)
    k = k+1
  }
  #-------------------------- Get Size of each cluster (csv file)---------------------------------
  extension.csv <- "csv"
  cluste.res <- Sys.glob(paste("*.", extension.csv, sep = ""))
  length(cluste.res)
  data.point <- data.frame()
  k=1
  ls.colnames <- vector()
  for(cn in cluste.res){
    csv.p <- paste(res.path, cn, sep = "/")
    clustInfo <- unlist(get.clusterSize(x = csv.p))
    data.point[1,k] <- clustInfo[2]
    ls.colnames[k] <- clustInfo[1]
    k = k+1
  }
  colnames(data.point) <- ls.colnames
  # fitness.cr
  # data.point
  cluster.size <- data.point[,colnames(fitness.cr)]
  ga.result <- rbind(cluster.size,fitness.cr)
  ga.result.temp <- as.data.frame(lapply(ga.result, function(x) as.numeric(as.character(x))))
  # get cluster size and fitness values: rows 1,11
  ga.fit.size <- ga.result.temp[c(1,11),]
  # ga.fit.size
  # 
  avg.error <- data.frame()
  colname.avg <- vector()
  for (i in 1:10){
    if(i ==1)
    {
      filter.rows <- ga.fit.size[,which(grepl( "Brito_Felizola" , colnames(ga.fit.size)))]
      avg.fitness <- sum(filter.rows[1]*filter.rows[2])/(sum(filter.rows[1]))
    }
    else {
    clname <- paste(paste(".",i,sep = ""),"_",sep = "")
    filter.rows <- ga.fit.size[,which(grepl( clname , colnames(ga.fit.size)))]
    avg.fitness <- sum(filter.rows[1,]*filter.rows[2,])/(sum(filter.rows[1,]))
    }
    # avg.error[1,i-1] <- avg.fitness
    # colname.avg[i-1] <- paste("k",i,sep = "=")
    avg.error[1,i] <- avg.fitness
    colname.avg[i] <- paste("k",i,sep = "=")
   
  }
   avg.error
   colname.avg
  colnames(avg.error) <- colname.avg
  # avg.error
  #********************* second type error ***************************
  mean.error <- data.frame()
  for (i in 2:10)
  {
    # if (i == 1)
    # {
    #   filter.rows <- ga.fit.size[, which(grepl("Brito_Felizola" , colnames(ga.fit.size)))]
    #   mean.fitness <- mean(as.numeric(filter.rows[2]))
    # }
    # else {
      clname <- paste(paste(".", i, sep = ""), "_", sep = "")
      filter.rows <- ga.fit.size[, which(grepl(clname , colnames(ga.fit.size)))]
      mean.fitness <- mean(as.numeric(filter.rows[2, ]))
    #}
    mean.error[1,i-1] <- mean.fitness
    #mean.error[1, i] <- mean.fitness
  }
  
  #***************************** Median ******************************
  median.error <- data.frame()
  for (i in 2:10)
  {
    # if (i == 1)
    # {
    #   filter.rows <- ga.fit.size[, which(grepl("Brito_Felizola" , colnames(ga.fit.size)))]
    #   median.fitness <- median(as.numeric(filter.rows[2]))
    # }
    # else {
      clname <- paste(paste(".", i, sep = ""), "_", sep = "")
      filter.rows <- ga.fit.size[, which(grepl(clname , colnames(ga.fit.size)))]
      median.fitness <- median(as.numeric(filter.rows[2, ]))
    #}
    median.error[1,i-1] <- median.fitness
  }
  
  colnames(mean.error) <- colname.avg
  colnames(median.error) <- colname.avg
  resName <- paste(path_file(path_dir(res.path)),unlist(strsplit(path_file(res.path),split = "_"))[2],sep = "_")
  avg.error$ResultType <- paste("WeightedAvg",resName,sep = "_")
  mean.error$ResultType <- paste("Mean",resName,sep="_")
  median.error$ResultType = paste("Median",resName,sep="_")
  saveRDS(rbind(avg.error, mean.error,median.error),
          file = paste(paste(path_dir(path_dir(res.path)), resName, sep = "/"),".rds",sep = ""))
  return(fitness.cr)
}
# Brito_Felizola_HCPC <- readRDS("C:/Users/sara/OneDrive - University of Tulsa/MyPhD/Results/Normal/Brito_Felizola_HCPC.rds")
# Brito_Felizola_HCPC$K = 0.06431895
#fitness.tuffp <- create.error.result(res.path = ga.path)
fitness.lm.tuffp <- create.error.result(res.path = ga.path)
write.csv(file = "C:/Users/sara/Desktop/res/fitness.lm.tuffp.csv",x = fitness.lm.tuffp)

#########################################################################
dim(fitness.lm.tuffp)
dat <- as.data.frame(rbind(HC = as.numeric(fitness.lm.tuffp[10,])#, HC_LM = as.numeric(fitness.tuffp.hc.lm[10,])
                           ))
# dat <- as.data.frame( rbind(HC = as.numeric(fitness.tuffp.hc[1,]),HC_LM = as.numeric(fitness.tuffp.hc[1,]),
#                             EM = as.numeric(fitness.bf.em[1,]),Fuzzy = as.numeric(fitness.bf.fuzzy[1,])))
colnames(dat) = colnames(fitness.lm.tuffp)
dat$method <- rownames(dat)


###########################################################################
dat.melt = melt(data = dat)
data.dens <- dat.melt[order(dat.melt$method),]
ggplot(data.dens, aes(x = value, fill = method)) + geom_density(alpha = 0.3)+
  scale_x_continuous(breaks = seq(-.5, 400, ), limits = c(-.5, 1)
                     ) +
  scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 20)
                     )
  
# plot(density(as.numeric(fitness.bf.km[1,])))
# normalmixEM(as.numeric(fitness.bf.km[1,]), lambda = .5)

#************************************************************************************************************************
Brito_Felizola_KM$`k=1`
Brito_Felizola_HCPC$`K` = Brito_Felizola_KM$`k=1`
Brito_Felizola_EM$`K` = Brito_Felizola_KM$`k=1`
Brito_Felizola_Fuzzy$`K` = Brito_Felizola_KM$`k=1`

colnames(Brito_Felizola_HCPC)[11] = "k=1"
colnames(Brito_Felizola_EM)[11] = "k=1"
colnames(Brito_Felizola_Fuzzy)[11] = "k=1"

Brito_Felizola_HCPC = Brito_Felizola_HCPC[,colnames(Brito_Felizola_KM)]
Brito_Felizola_EM = Brito_Felizola_EM[,colnames(Brito_Felizola_KM)]

################################################################################
 
#data.plot <- rbind(avg.error, mean.error,median.error)
data.plot  = Brito_Felizola_KM
data.anlysis <- data.plot
colnames(data.anlysis) <- c(1:10,"ResultType")
data.anlysis$ResultType = c("Weighted Mean","Mean","Median")
data.melt <- melt(data=data.anlysis)
colnames(dat.melt) <- c("ResultType","k","FitnessError")

#################################################################################
GA.bf.mean.w.median <- dat.melt %>% 
  ggplot(data = dat.melt, mapping = aes(group  = ResultType)) +
  labs(x = "Number of Clusters", y = "Fitness Error")+
  geom_point(aes(x =k, y = FitnessError, colour = ResultType)) +
  geom_line(mapping = aes(x =k, y = FitnessError, colour = ResultType)) +
  scale_y_continuous(#trans = "log2",
                     labels=scaleFUN)+
  # geom_label(data = dat.melt %>%
  #              filter(FitnessError %in% min(FitnessError[which(grepl("Mean" , ResultType))]) |
  #                       FitnessError %in% min(FitnessError[which(grepl("Weighted" , ResultType))]) |
  #                       FitnessError %in% min(FitnessError[which(grepl("Median" , ResultType))]))
  #            ,size=3,label.size = .1,fill = NA,mapping = aes(x =k, y = FitnessError, label = "min"))+
  # scale_color_manual(#labels = c("Mean", "Weighted","Median"), 
  #                    values = colours) +
  theme_light() + theme(legend.title = element_blank(),legend.position =  c(0.86, 0.89),
                        legend.key.size = unit(0.2, "cm"),legend.key.height = unit(0.2, "cm"),legend.key.width = unit(0.2, "cm"))

library(gridExtra)
library(ggpubr)
png("C:/Users/sara/Desktop/Fig_4_15(2).png", units="in", width=4, height=4,res = 300)
#GA.mean.bf
grid.arrange(arrangeGrob(
  GA.bf.mean.w.median + theme(legend.position = c(0.8, 0.89))
  ,nrow = 1), heights = c(10, 1))

dev.off()

################################################################################
# load rds files

LM_KM$`K` = Brito_Felizola_KM$`k=1`
LM_Fuzzy$`K` = Brito_Felizola_KM$`k=1`
LM_EM$`K` = Brito_Felizola_KM$`k=1`
LM_HCPC$`K` = Brito_Felizola_KM$`k=1`
colnames(LM_KM)[11] = "k=1"
colnames(LM_Fuzzy)[11] = "k=1"
colnames(LM_EM)[11] = "k=1"
colnames(LM_HCPC)[11] = "k=1"
LM.KM = LM_KM[,colnames(Brito_Felizola_KM)]
LM.EM = LM_EM[,colnames(Brito_Felizola_KM)]
LM.Fuzzy = LM_Fuzzy[,colnames(Brito_Felizola_KM)]
LM.HCPC = LM_HCPC[,colnames(Brito_Felizola_KM)]
#-------------------------------------------------------------------------------
Barnea_KM$`K` = Brito_Felizola_KM$`k=1`
Barnea_Fuzzy$`K` = Brito_Felizola_KM$`k=1`
Barnea_EM$`K` = Brito_Felizola_KM$`k=1`
Barnea_HCPC$`K` = Brito_Felizola_KM$`k=1`
colnames(Barnea_KM)[11] = "k=1"
colnames(Barnea_Fuzzy)[11] = "k=1"
colnames(Barnea_EM)[11] = "k=1"
colnames(Barnea_HCPC)[11] = "k=1"

Barnea_KM = Barnea_KM[,colnames(Brito_Felizola_KM)]
Barnea_EM = Barnea_EM[,colnames(Brito_Felizola_KM)]
Barnea_Fuzzy = Barnea_Fuzzy[,colnames(Brito_Felizola_KM)]
Barnea_HCPC = Barnea_HCPC[,colnames(Brito_Felizola_KM)]

#-------------------------------------------------------------------------------
DR_KM$`K` = Brito_Felizola_KM$`k=1`
DR_Fuzzy$`K` = Brito_Felizola_KM$`k=1`
DR_EM$`K` = Brito_Felizola_KM$`k=1`
DR_HCPC$`K` = Brito_Felizola_KM$`k=1`
colnames(DR_KM)[11] = "k=1"
colnames(DR_Fuzzy)[11] = "k=1"
colnames(DR_EM)[11] = "k=1"
colnames(DR_HCPC)[11] = "k=1"
DR_KM = DR_KM[,colnames(Brito_Felizola_KM)]
DR_EM = DR_EM[,colnames(Brito_Felizola_KM)]
DR_Fuzzy = DR_Fuzzy[,colnames(Brito_Felizola_KM)]
DR_HCPC = DR_HCPC[,colnames(Brito_Felizola_KM)]

#---------------------------------------------------------------------------------
plot.4clust.mean.weight.median <- function(data.ploting, plot.title, transformation= "log" ){
  
  plot.res <- data.ploting %>% 
    ggplot(data = data.ploting, mapping = aes(group  = ResultType)) +
    labs(x = "Number of Clusters", y = "Average Error (log transformation)",title = plot.title)+
    geom_point(aes(x =k, y = FitnessError, colour = ResultType)) +
    geom_line(mapping = aes(x =k, y = FitnessError, colour = ResultType)) +
    scale_y_continuous(trans = transformation,labels=scaleFUN)+
    geom_label(data = data.ploting %>%
                 filter(FitnessError %in% min(FitnessError[which(grepl("Fuzzy" , ResultType))]) |
                          FitnessError %in% min(FitnessError[which(grepl("EM" , ResultType))]) |
                          FitnessError %in% min(FitnessError[which(grepl("KM" , ResultType))]) |
                          FitnessError %in% min(FitnessError[which(grepl("HCPC" , ResultType))]))
               ,size=3,label.size = .1,fill = NA,mapping = aes(x =k, y = FitnessError, label = "min"))+
    scale_color_manual(labels = c("EM", "Fuzyy","Hierarchical","K-means"), 
                       values = colours) +
    theme_light() + theme(legend.title = element_blank(),legend.position = "bottom",
                          legend.key.size = unit(0.2, "cm"),legend.key.height = unit(0.2, "cm"),legend.key.width = unit(0.2, "cm"))
  return(plot.res)
}

#################################################################################
rm(data.plot)
#data.plot <- rbind(Brito_Felizola_KM,Brito_Felizola_Fuzzy,Brito_Felizola_EM,Brito_Felizola_HCPC)

# data.plot <- rbind(LM.KM,LM.Fuzzy,LM.EM,LM.HCPC)

# data.plot <- rbind(Barnea_KM,Barnea_Fuzzy,Barnea_EM,Barnea_HCPC)

data.plot <- rbind(DR_KM,DR_Fuzzy,DR_EM,DR_HCPC)

rm(mean.data.bf,w.avg.data.bf,median.data.bf,mean.melt.data,median.melt.data,weightedAvg.melt.data)

#*************************************************************************************************

mean.data.bf <- data.plot[grep(pattern = "Mean", data.plot$ResultType),]
w.avg.data.bf <- data.plot[grep(pattern = "WeightedAvg", data.plot$ResultType),]
median.data.bf <- data.plot[grep(pattern = "Median", data.plot$ResultType),]

colnames(median.data.bf) = c(1:10, "ResultType")
median.melt.data <- melt(data=median.data.bf)
colnames(median.melt.data) <- c("ResultType","k","FitnessError")

# median.bf.dimensional <- plot.4clust.mean.weight.median(data.ploting = median.melt.data,
#                                plot.title = "Median of GA Error on Dimensional",transformation ="log" )

# median.bf.lm <- plot.4clust.mean.weight.median(data.ploting = median.melt.data,
#                                                          plot.title = "Median of GA Error on LM",transformation ="log2" )

# median.bf.barnea <- plot.4clust.mean.weight.median(data.ploting = median.melt.data,
#                                                plot.title = "Median of GA Error on Barnea",transformation ="log2" )
                                               
median.bf.dr <- plot.4clust.mean.weight.median(data.ploting = median.melt.data,
                                                   plot.title = "Median of GA Error on Duns_Ros",transformation ="log2" )

#***********************************************************************************************************
colnames(mean.data.bf) = c(1:10, "ResultType")
mean.melt.data <- melt(data=mean.data.bf)
colnames(mean.melt.data) <- c("ResultType","k","FitnessError")
# mean.bf.dimensional <- plot.4clust.mean.weight.median(data.ploting = mean.melt.data,
#                                                         plot.title = "Mean of GA Error on Dimensional",transformation ="log" )

# mean.bf.lm <- plot.4clust.mean.weight.median(data.ploting = mean.melt.data,
#                                                plot.title = "Mean of GA Error on LM",transformation ="log" )

# mean.bf.barnea <- plot.4clust.mean.weight.median(data.ploting = mean.melt.data,
#                                                plot.title = "Mean of GA Error on Barnea",transformation ="log2" )

mean.bf.dr <- plot.4clust.mean.weight.median(data.ploting = mean.melt.data,
                                               plot.title = "Mean of GA Error on Dun_Ros",transformation ="log2" )


#*********************************************************************************************************
colnames(w.avg.data.bf) = c(1:10, "ResultType")
weightedAvg.melt.data <- melt(data=w.avg.data.bf)
colnames(weightedAvg.melt.data) <- c("ResultType","k","FitnessError")
# weightedAvg.bf.dimensional <- plot.4clust.mean.weight.median(data.ploting = weightedAvg.melt.data,
#                                                              plot.title = "Weighted Mean of GA Error on Dimensional",transformation ="log" )

# weightedAvg.bf.lm <- plot.4clust.mean.weight.median(data.ploting = weightedAvg.melt.data,
#                                                plot.title = "Weighted Mean of GA Error on LM",transformation ="log2" )

# weightedAvg.bf.barnea <- plot.4clust.mean.weight.median(data.ploting = weightedAvg.melt.data,
#                                                plot.title = "Weighted Mean of GA Error on Barnea",transformation ="sqrt" )

weightedAvg.bf.dr <- plot.4clust.mean.weight.median(data.ploting = weightedAvg.melt.data,
                                               plot.title = "Weighted Mean of GA Error on Duns_Ros",transformation ="log" )
#****************************************************************************************************
library(gridExtra)
mean.bf.dimensional
median.bf.dimensional
weightedAvg.bf.dimensional

mylegend <- g_legend(mean.bf.dimensional)

png("C:/Users/sara/Desktop/weightedAvg_bf.png", units="in", width=8, height=8,res = 500)
#GA.mean.bf
grid.arrange(arrangeGrob(
  weightedAvg.bf.lm + theme(legend.position = "none"),weightedAvg.bf.barnea + theme(legend.position = "none"),
  weightedAvg.bf.dr + theme(legend.position = "none"), weightedAvg.bf.dimensional + theme(legend.position = "none"),
  nrow = 2), mylegend, heights = c(10, 1))

dev.off()



#****************************************************************************************************
GA.median.bf <- median.melt.data %>% 
  ggplot(data = median.melt.data, mapping = aes(group  = ResultType)) +
  labs(x = "Number of Clusters", y = "Fitness Error",title = "Median of GA Error for Dimensional")+
  geom_point(aes(x =k, y = FitnessError, colour = ResultType)) +
  geom_line(mapping = aes(x =k, y = FitnessError, colour = ResultType)) +
  scale_y_continuous(trans = "log2",
                     labels=scaleFUN)+
  geom_label(data = median.melt.data %>%
               filter(FitnessError %in% min(FitnessError[which(grepl("Fuzzy" , ResultType))]) |
                        FitnessError %in% min(FitnessError[which(grepl("EM" , ResultType))]) |
                        FitnessError %in% min(FitnessError[which(grepl("KM" , ResultType))]) |
                        FitnessError %in% min(FitnessError[which(grepl("HCPC" , ResultType))]))
             ,size=3,label.size = .1,fill = NA,mapping = aes(x =k, y = FitnessError, label = "min"))+
  scale_color_manual(labels = c("EM", "Fuzyy","Hierarchical","K-means"), 
                     values = colours) +
  theme_light() + theme(legend.title = element_blank(),legend.position = "bottom",
                        legend.key.size = unit(0.2, "cm"),legend.key.height = unit(0.2, "cm"),legend.key.width = unit(0.2, "cm"))


################################################################################

GA.weighted.mean.bf.lm
GA.weighted.mean.bf.barnea
GA.weighted.mean.bf.dr
GA.weighted.mean.bf
#**************************
GA.mean.bf.dr
GA.mean.bf.barnea
GA.mean.bf.lm
GA.mean.bf
#**************************
GA.media.bf
GA.median.bf.lm
GA.median.bf.barnea
GA.median.bf.dr
#**************************

library(gridExtra)

mylegend <- g_legend(GA.mean.bf)

png("C:/Users/sara/Desktop/phd/chapter5/Dless_Validation/GA_SimpleMean_bf.png", units="in", width=8, height=8,res = 500)
#GA.mean.bf
grid.arrange(arrangeGrob(
  GA.mean.bf.lm + theme(legend.position = "none"),GA.mean.bf.barnea + theme(legend.position = "none"),
  GA.mean.bf.dr + theme(legend.position = "none"), GA.mean.bf + theme(legend.position = "none"),
  nrow = 2), mylegend, heights = c(10, 1))

dev.off()
####################################################################
# plot mean, weighred avg, median Brito-Felizola (km, em, hc, fuzzy)
mean.melt.data
median.melt.data
weightedAvg.melt.data
mean.ind.val %>% 
  ggplot(aes(k, indexVal, label = type)) + 
  labs(x = "Number of Clusters", y = "Average Validity Index",title ="Hierarchical Method on Duns-Ros") +
  geom_point(aes(x = k, y = indexVal, colour = type)) +
  geom_line(mapping = aes(x = k, y = indexVal, colour = type)) +
  geom_label(data = mean.ind.val %>% 
               filter(indexVal %in% min(indexVal[which(type =="min")]) | indexVal %in% max(indexVal[which(type =="max")])),
             size=3,label.size = .1,fill = NA)+
  scale_x_continuous(breaks = seq(1, 10, 1), limits = c(2, 10)) +
  scale_y_continuous(#trans = "log2",
    labels=sciFun)+
  scale_colour_manual(values = colours[c(1,4)])+
  theme_light() + theme(legend.title = element_blank(),legend.position = "bottom",legend.key.size = unit(0.2, "cm"),legend.key.height = unit(0.2, "cm"),legend.key.width = unit(0.2, "cm"))

####################################################################

#---------------------------------------------------------------------------------
ggplot(data = mean.melt.data, aes(group  = ResultType)) +
  labs(x = "Number of Clusters", y = "Average Error")+
  geom_point(aes(x =k, y = FitnessError, colour = ResultType)) +
  geom_line(mapping = aes(x =k, y = FitnessError, colour = ResultType)) +
  #scale_x_continuous(breaks = seq(1, 10, 1), limits = c(2, 10)) +
  scale_y_continuous(labels=scaleFUN)+
  scale_colour_manual(values = colours)+
  theme_light() + theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm")
  )
#-------------------------------------------------------------------------------
library(magrittr) 
library(dplyr) 
melt.data %>% 
  ggplot(aes(x =k, y = FitnessError, label = "min",group  = ResultType)) + 
  labs(x = "Number of Clusters", y = "Average Error") +
  geom_point(aes(x =k, y = FitnessError, colour = ResultType)) +
  geom_line(mapping = aes(x =k, y = FitnessError, colour = ResultType)) +
  geom_label(data = melt.data %>% 
               filter(FitnessError %in% min(FitnessError)),
             size=3,label.size = .1,fill = NA)+
  #scale_x_continuous(breaks = seq(1, 10, 1), limits = c(2, 10)) +
  scale_y_continuous(labels=scaleFUN)+
  scale_colour_manual(values = colours)+
  theme_light() + theme(legend.title = element_blank(),legend.position = "bottom",
                        legend.key.size = unit(0.2, "cm"),legend.key.height = unit(0.2, "cm"),legend.key.width = unit(0.2, "cm"))

#****************************************************************************
# txt.f <- "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/Brito_Felizola_FlowPattern_DM/DR_clusters"
# setwd(txt.f)
# extension <- "txt"
# fNames <- Sys.glob(paste("*.", extension, sep = ""))
# subName <- substr(fNames,start = 10, stop = nchar(fNames) )
# colN <- vector()
# for(i in seq(subName))
# {
#   splitName <- unlist(strsplit(subName[i], split = "_"))
#   colN[i] <- paste(splitName[1], splitName[2],sep = "_")
# }
#******************************************************************************************************
plot.data <- read.csv("C:/Users/sara/Desktop/My PhD Research/ABM/Time_Error.csv")

library(ggplot2)
ggplot(plot.data, aes(x = K)) + 
  geom_line(aes(y = AvgError)) + 
  geom_line(aes(y = runtime))

ggplot(plot.data, aes(x = K)) +
  # geom_line(aes(y = AvgError, colour = "Average Error")) +
  # geom_line(aes(y = runtime / 500, colour = "Running Time")) +
  geom_smooth(aes(y = AvgError, colour = "Average Error")) +
  geom_smooth(aes(y = runtime / 500, colour = "Running Time"))+
  scale_y_continuous(sec.axis = sec_axis(~ . * 500, name = "Running Time")) +
  scale_colour_manual(values = c("blue", "red")) +
  labs(y = "Average Error", x = "Cluster #") +
  scale_x_discrete(limits = seq(10)) +
  theme_grey()+
  theme(
    legend.position = c(0.186, 0.914),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"))

#theme_light()
# theme_grey()
#theme_get()

#********************** get the fitness values (error) for each text files (result og GA) + the cluster name

txt.files <- "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/Brito_Felizola_FlowPattern_DM/DR_clusters"
files <- list.files(path=txt.files, pattern="*.txt", full.names=TRUE, recursive=FALSE)
fitness.vals <- lapply(files,get.fitt )
substr(files[1],start = 40, stop = nchar(files[1])-40 )
fitnessVal <- as.data.frame(unlist(fitness.vals),stringsAsFactors = FALSE)
fitnessVal$clname <- rownames(fitnessVal)
colnames(fitnessVal) <- c("fitnessValue" ,"clname")
#fitnessVal <- round(fitnessVal,digits = 5)
#***************************************************

csv.files.path <- "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/Brito_Felizola_FlowPattern_DM/DR_clusters"
csv.files <- list.files(path=csv.files.path, pattern="*.csv", full.names=TRUE, recursive=FALSE)

get.invalid.rows <- function(dataf) {
  l = 1
  index <- list()
  for (i in 1:dim(dataf)[1]) {
    colnum <- which(dataf[i, ] == (-999))
    index[[l]] <- c(colnum = ifelse(length(colnum) > 0, colnum, 0),rownum = l)
    l = l + 1
  }
  index.df <- as.data.frame(do.call(rbind, index))
  invalid.rows <- index.df$rownum[index.df$colnum != 0]
  return (invalid.rows)
} # end of function


um.cols <-c("Liq1DensityGL","GasDensityGL","Liq1ViscosityGL","GasViscosityGL",
            "GasLiqSurTensionGL","PipeDiaGL","PipeIncliGL","PipeRoughGL",
            "Liq1SuperficialVelGL","GasSuperficialVelGL","PressureDropGL")


#*********** MAKE A BACKUP COPY FROM YOUR FILES BEFORE THIS *********************

for(i in seq(csv.files)){
  clustcsv <- read.csv(csv.files[i],header = TRUE,fill = TRUE, stringsAsFactors = FALSE)
  clustcsv <- clustcsv[,-1]
  invalR <-  get.invalid.rows(dataf =clustcsv[,um.cols] )
  if (length(invalR)>0) {
    csv.new <- clustcsv[-invalR, ]
  } else {
    csv.new <- clustcsv
  }
  write.csv(csv.new, file = csv.files[i])
}
#***************************************************
datapoints.numbers <- lapply(csv.files,get.datapoint.number )
sum(unlist(datapoints.numbers)[1:10])
unlist(datapoints.numbers)
data.p.no <- as.data.frame(unlist(datapoints.numbers),stringsAsFactors = FALSE)
#data.p.no$clname <- fitnessVal$clname

# *** n = number of data points in each cluster
colnames(data.p.no) <- "n"
cbind(data.p.no,fitnessVal)

datas <- as.data.frame(cbind(data.p.no,fitnessVal),stringsAsFactors = FALSE)
write.csv(datas, "C:/Users/sam805/Desktop/datas.csv")
getwd()

#clust10 <- grep("^10", colnames(datas))
kval <- as.numeric(substr(datas$clname,start = 1,stop =nchar(rownames(datas))-2 ))
kval[2 ] = 10
new.datas <- cbind(kval,datas)
colnames(new.datas) <- c("K","n","FitnessError","clname")
new.datas$multi <- new.datas[,2] * new.datas[,3]
Ks <- unique(new.datas$K)
w.list <- list()
j <- 1
for(i in Ks){
  req.data <- new.datas[which(new.datas$K== i),]
  w.list[[j]] <- c(weightedSum = sum(req.data$multi)/sum(req.data$n),Kvalue = i)
  j <- j+1
}
pre.wdf <- as.data.frame(do.call(rbind, w.list))
w.df <- pre.wdf[order(pre.wdf$Kvalue),]
colnames(w.df) <- c("AvgError", "K")
#Brito_Felizola Error: 0.0521 
df.avgEr <- rbind(c(0.0521, 1),w.df)

#************************* Plot Average Error VS K-Values ******************
library(ggplot2)
ggplot(df.avgEr)
ggplot(df.avgEr, aes(x = K, y = AvgError))+
  #geom_smooth(se = FALSE)+
  geom_point(color = "green")+
  geom_line(color="blue")+
  scale_x_continuous() +
  scale_y_continuous() +
  scale_colour_discrete()

#***************************************************************************
#********* ADD RUNTIME **************************
Wdata <- read.csv("C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/Brito_Felizola_FlowPattern_DM/DRrunTime.csv",header = FALSE)

colnames(new.datas)
#[1] "K"            "n"            "FitnessError" "clname"       "multi"   

Wdata$V1 <- fitnessVal$clname
Wdata$K <- new.datas$K

Ks <- unique(Wdata$K)
t.list <- list()
j <- 1
for(i in Ks){
  req.data <- Wdata[which(Wdata$K== i),]
  t.list[[j]] <- c(totTime = sum(req.data$V2),Kvalue = i)
  j <- j+1
}
pre.t <- as.data.frame(do.call(rbind, t.list))
t.df <- pre.t[order(pre.t$Kvalue),]
colnames(t.df) <- c("time", "K")
t.df
w.df

# adding k=1 (218 datapoints of brito and felizola)
df.avgEr

#DFsumErr <- as.data.frame(cbind(unique(Wdata[,2]),na.omit(Wdata[,6])))
DFsumErr <- as.data.frame(cbind(time = t.df$time,w.df))

# Brito_felizola results
err.df <- rbind( c(500.736,0.0521, 1),DFsumErr)
#df.avgEr[1,]
#colnames(DFsumErr) <- c("K", "ErrorSum")
# pre.time <- apply(output.file,MARGIN = 1,strsplit,split = "/")
# as.data.frame(pre.time,stringsAsFactors = FALSE)

ggplot(DFsumErr[,1:2])
dim(DFsumErr)

ggplot(DFsumErr, aes(x = DFsumErr$K, y = DFsumErr$AvgError))+
  #geom_smooth(se = FALSE)+
  geom_point(color = "green")+
  geom_line(color="blue")+
  scale_x_continuous() +
  scale_y_continuous() +
  scale_colour_discrete()

#******************************************************************************
#*************************** PLOT OF AVERAGE ERROR VS K ***********************
#******************************************************************************

library(reshape2)
df.melt = melt(err.df, id="K")
ggplot(data = df.melt,aes(x = err.df$AvgError, y = err.df$K, colour = variable)) + 
  #geom_point() +
  #geom_line() + 
  xlab("K") + 
  ylab("Average Error") + 
  #xlim(1, 10) + 
  scale_x_discrete(limits =seq(10))+
  scale_y_continuous(breaks =seq(0.025, 0.055, 0.005), limits = c(0.025, 0.055))+
  geom_line(aes(group=2), colour="blue") +  # Blue lines
  geom_point(size=1, colour="red")

#******************************************************************************
#******************************************************************************

RunTimeGA <- read.csv("C:/Users/sara/Desktop/My PhD Research/ClusteringRelated/RunTimeGA.csv",header = FALSE,stringsAsFactors = FALSE)
head(RunTimeGA)
dim(RunTimeGA)
DFtime <- as.data.frame(cbind(1:10,na.omit(RunTimeGA[,3])))
colnames(DFtime) <- c("K", "GATotalTime")
#geom_line(data =df.runT, aes(x=K,y=DFtime$GATotalTime))

plot(DFtime,type = "l")

p <- ggplot() +
  # blue plot: 
  geom_line(data=DFsumErr, aes(x=K, y=ErrorSum)) + 
  geom_smooth(data=DFsumErr, aes(x=K, y=ErrorSum),
              colour="darkblue", size=1) +
  # red plot
  geom_line(data=DFtime, aes(x=K, y=GATotalTime)) + 
  geom_smooth(data=DFtime, aes(x=K, y=GATotalTime),
              colour="red", size=1)

DFtime.log<- as.data.frame(cbind(1:10, log10(log(log(DFtime[,2])))))

# make negative values 0
DFtime.log[1,2] = .0001
DFtime.log[2,2] = .0001
colnames(DFtime.log) <- c("K", "GATotalTime")

p = ggplot() + 
  geom_line(data = DFsumErr, aes(x = K, y = ErrorSum), color = "darkblue") +
  geom_line(data = DFtime.log, aes(x = K, y = GATotalTime), color = "darkred") + 
  scale_y_continuous(trans = "log1p") +
  scale_x_discrete(limits =seq(10))+
  geom_smooth(data=DFsumErr, aes(x=K, y=ErrorSum),colour="blue")+
  geom_smooth(data=DFtime.log, aes(x=K, y=GATotalTime),colour="red")+
  #scale_y_continuous(breaks =seq(0.025, 0.055, 0.005), limits = c(0.025, 0.055))+
  xlab('Number of Clusters') +
  ylab('Error/CPU Time')+
  scale_colour_manual(values=c("green","blue"))+
  theme(
    legend.position = c("left", "top"),
    legend.justification = c("left", "top"))

print(p)


methErr <- as.data.frame(cbind(1:5, c(0.03863165,0.03485083,0.03696972,0.03639312,0.03876284)))
colnames(methErr) <- c("Cluster_Algo", "Avg_Error")
rownames(methErr) <- c("Kmeans","HCPC","EM","Fuzzy","SOM")

ggplot(methErr)
ggplot(methErr, aes(x = Cluster_Algo, y = Avg_Error))+
  geom_bar(stat="identity")


ggplot(methErr, aes(x = Cluster_Algo, y = Avg_Error))+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=rownames(methErr)), vjust=1.6, color="white", size=2.5)+
  theme_minimal()

p<-ggplot(methErr, aes(x=Cluster_Algo, y=Avg_Error, fill=rownames(methErr))) +
  geom_bar(stat="identity")+theme_minimal()
p
#*********************************************************************************************
#x =  "C:/Users/sara/OneDrive - University of Tulsa/MyPhD/Results/Dimensionless/WholeData/LM/Res_HCPC/hcpc.6_1_1000pop_0.03m.txt"
x = "C:/Users/sara/Desktop/hcpc.2_1_1000pop_0.03m.txt"
#x= "C:/Users/sara/OneDrive - University of Tulsa/MyPhD/Results/Normal/tuffp_2000pop_0.03m.txt"
x= "C:/Users/sara/OneDrive - University of Tulsa/MyPhD/hcpc.6_1_1000pop_0.03m.txt"


fit.error <- data.frame()
opt.cr <- vector()

output.file <- read.table(x,header = TRUE,fill = TRUE,stringsAsFactors = FALSE)
iter.n <- ceiling(dim(output.file)[1]/11)
only.error <- data.frame()
fit.ind <- 9
dim(output.file)
for(i in 1:iter.n){
  only.error[i,1:dim(output.file)[2]] <- as.numeric(output.file[fit.ind,])
  fit.ind = fit.ind + 11
} 
min.fit.ind <- apply(only.error, 1, which.min)[iter.n]


#sort(only.error[1,], decreasing = FALSE)[15]


apply(only.error, 1, min)[iter.n]
len <- dim(output.file)[1]
CR.err <- output.file[(len-10):(len-1),min.fit.ind]
CR.err

rows.index <- dim(output.file)[1] - 10
ord.ind <- order(output.file[nrow(output.file) - 1, ],decreasing = TRUE)
sort.output <- output.file[, ord.ind]
opt.cr <- sort.output[rows.index:(rows.index + 8), 1]
min.fitness <- round(as.double(sort.output[nrow(sort.output) - 1, ])[1],digits = 6)
fit.error <- as.data.frame(unlist(list(min.fitness, opt.cr)))
colnames(fit.error) <- element.name
