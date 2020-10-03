
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# calinski_harabasz,dunn,gdi11,gdi12,gdi13,gdi21,gdi22,gdi23,gdi31,gdi32,gdi33,gdi41,gdi42,gdi43,gdi51,
# gdi52,gdi53,gamma,pbm,point_biserial,ratkowsky_lance,silhouette,tau,wemmert_gancarski ---->>>>  MAX
#banfeld_raftery,c_index,davies_bouldin,g_plus,mcclain_rao,ray_turi,scott_symons,sd_scat,sd_dis,s_dbw,xie_beni------>>MIN

#allInd <- c("c_index","calinski_harabasz","davies_bouldin","det_ratio","dunn","s_dbw","silhouette","xie_beni")
# allInd <- c("banfeld_raftery","c_index","davies_bouldin","g_plus","mcclain_rao","ray_turi",
#             "scott_symons","sd_scat","sd_dis" ,"s_dbw","xie_beni",
#             #Max Indices
#             "calinski_harabasz","dunn","gamma","gdi11","gdi12","gdi13","gdi21","gdi22","gdi23","gdi31","gdi32","gdi33","gdi41",
#             "gdi42","gdi43","gdi51","gdi52","gdi53","pbm","point_biserial","ratkowsky_lance","silhouette","tau","wemmert_gancarski")  

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

library(clv)
library(cluster)
library(NbClust)
library(clValid)
library(clustree)
library(dendextend)
library(factoextra)
library(FactoMineR)


intraclust = c("complete","average","centroid")
interclust = c("single", "complete", "average","centroid", "aveToCent", "hausdorff")

# define new Dunn and Davies.Bouldin functions
Dunn <- function(data,clust) 
  clv.Dunn( cls.scatt.data(data,clust),
            intracls = c("complete","average","centroid"), 
            intercls = c("single", "complete", "average","centroid", "aveToCent", "hausdorff"))

Davies.Bouldin <- function(data,clust) 
  clv.Davies.Bouldin( cls.scatt.data(data,clust),
                      intracls = c("complete","average","centroid"),
                      intercls = c("single", "complete", "average","centroid", "aveToCent", "hausdorff"))

my.cluster.stats <- function (d ,clustering ,
                              alt.clustering = NULL,noisecluster = FALSE,silhouette = TRUE,
                              G2 = FALSE,G3 = FALSE,wgap = TRUE,sepindex = TRUE,sepprob = 0.1,
                              sepwithnoise = TRUE,compareonly = FALSE,aggregateonly = FALSE)
{
  if (!is.null(d))
    d <- as.dist(d)
  
  if(class(clustering)=="factor"){
    cn <- max(as.integer(levels(clustering)))
  } else 
    cn <- max(clustering)
  
  clusteringf <- as.factor(clustering)
  clusteringl <- levels(clusteringf)
  cnn <- length(clusteringl)
  if (cn != cnn) {
    warning("clustering renumbered because maximum != number of clusters")
    for (i in 1:cnn)
      clustering[clusteringf == clusteringl[i]] <- i
    cn <- cnn
  }
  n <- length(clustering)
  noisen <- 0
  cwn <- cn
  if (noisecluster) {
    noisen <- sum(clustering == cn)
    cwn <- cn - 1
  }
  diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <- within.dist <- between.dist <- numeric(0)
  for (i in 1:cn)
    cluster.size[i] <- sum(clustering == i)
  pk1 <- cluster.size / n
  pk10 <- pk1[pk1 > 0]
  h1 <- -sum(pk10 * log(pk10))
  corrected.rand <- vi <- NULL
  if (!is.null(alt.clustering)) {
    choose2 <- function(v) {
      out <- numeric(0)
      for (i in 1:length(v))
        out[i] <- ifelse(v[i] >= 2,
                         choose(v[i], 2), 0)
      out
    }
    cn2 <- max(alt.clustering)
    clusteringf <- as.factor(alt.clustering)
    clusteringl <- levels(clusteringf)
    cnn2 <- length(clusteringl)
    if (cn2 != cnn2) {
      warning("alt.clustering renumbered because maximum != number of clusters")
      for (i in 1:cnn2)
        alt.clustering[clusteringf == clusteringl[i]] <- i
      cn2 <- cnn2
    }
    nij <- table(clustering, alt.clustering)
    dsum <- sum(choose2(nij))
    cs2 <- numeric(0)
    for (i in 1:cn2)
      cs2[i] <- sum(alt.clustering == i)
    sum1 <- sum(choose2(cluster.size))
    sum2 <- sum(choose2(cs2))
    pk2 <- cs2 / n
    pk12 <- nij / n
    corrected.rand <- (dsum - sum1 * sum2 / choose2(n)) / ((sum1 +
                                                              sum2) / 2 - sum1 * sum2 /
                                                             choose2(n))
    pk20 <- pk2[pk2 > 0]
    h2 <- -sum(pk20 * log(pk20))
    icc <- 0
    for (i in 1:cn)
      for (j in 1:cn2)
        if (pk12[i, j] > 0)
          icc <- icc + pk12[i, j] * log(pk12[i, j] / (pk1[i] *
                                                        pk2[j]))
    vi <- h1 + h2 - 2 * icc
  }
  if (compareonly) {
    out <- list(corrected.rand = corrected.rand, vi = vi)
  }
  else {
    dmat <- as.matrix(d)
    within.cluster.ss <- 0
    overall.ss <- nonnoise.ss <- sum(d ^ 2) / n
    if (noisecluster)
      nonnoise.ss <- sum(as.dist(dmat[clustering <= cwn,clustering <= cwn]) ^ 2) / sum(clustering <= cwn)
    ave.between.matrix <- separation.matrix <- matrix(0,ncol = cn, nrow = cn)
    for (i in 1:cn) {
      cluster.size[i] <- sum(clustering == i)
      di <- as.dist(dmat[clustering == i, clustering ==
                           i])
      if (i <= cwn) {
        within.cluster.ss <- within.cluster.ss + sum(di ^ 2) / cluster.size[i]
        within.dist <- c(within.dist, di)
      }
      if (length(di) > 0)
        diameter[i] <- max(di)
      else
        diameter[i] <- NA
      average.distance[i] <- mean(di)
      median.distance[i] <- median(di)
      bv <- numeric(0)
      for (j in 1:cn) {
        if (j != i) {
          sij <- dmat[clustering == i, clustering ==
                        j]
          bv <- c(bv, sij)
          if (i < j) {
            separation.matrix[i, j] <- separation.matrix[j,
                                                         i] <-
              min(sij)
            ave.between.matrix[i, j] <- ave.between.matrix[j,
                                                           i] <-
              mean(sij)
            if (i <= cwn & j <= cwn)
              between.dist <- c(between.dist, sij)
          }
        }
      }
      separation[i] <- min(bv)
      average.toother[i] <- mean(bv)
    }
    #*****
    average.between <- mean(between.dist)
    average.within <- weighted.mean(average.distance, cluster.size,
                                    na.rm = TRUE)
    nwithin <- length(within.dist)
    nbetween <- length(between.dist)
    between.cluster.ss <- nonnoise.ss - within.cluster.ss
    ch <- between.cluster.ss * (n - noisen - cwn) / (within.cluster.ss * (cwn - 1))
    clus.avg.widths <- avg.width <- NULL
    if (silhouette) {
      if(class(clustering)=="factor")
        clustering <- as.numeric(clustering)
      sii <- silhouette(clustering, dmatrix = dmat)
      sc <- summary(sii)
      clus.avg.widths <- sc$clus.avg.widths
      if (noisecluster)
        avg.width <- mean(sii[clustering <= cwn, 3])
      else
        avg.width <- sc$avg.width
    }
    g2 <- g3 <- cn2 <- cwidegap <- widestgap <- sindex <- NULL
    if (G2) {
      splus <- sminus <- 0
      for (i in 1:nwithin) {
        splus <- splus + sum(within.dist[i] < between.dist)
        sminus <- sminus + sum(within.dist[i] > between.dist)
      }
      g2 <- (splus - sminus) / (splus + sminus)
    }
    if (G3) {
      sdist <- sort(c(within.dist, between.dist))
      sr <- nwithin + nbetween
      dmin <- sum(sdist[1:nwithin])
      dmax <- sum(sdist[(sr - nwithin + 1):sr])
      g3 <- (sum(within.dist) - dmin) / (dmax - dmin)
    }
    pearsongamma <- cor(c(within.dist, between.dist), c(rep(0,
                                                            nwithin), rep(1, nbetween)))
    dunn <-
      min(separation[1:cwn]) / max(diameter[1:cwn], na.rm = TRUE)
    acwn <- ave.between.matrix[1:cwn, 1:cwn]
    dunn2 <-
      min(acwn[upper.tri(acwn)]) / max(average.distance[1:cwn],
                                       na.rm = TRUE)
    if (wgap) {
      cwidegap <- rep(0, cwn)
      for (i in 1:cwn)
        if (sum(clustering == i) > 1)
          cwidegap[i] <- max(hclust(as.dist(dmat[clustering ==
                                                   i, clustering == i]), method = "single")$height)
        widestgap <- max(cwidegap)
    }
    if (sepindex) {
      psep <- rep(NA, n)
      if (sepwithnoise | !noisecluster) {
        for (i in 1:n)
          psep[i] <- min(dmat[i, clustering !=
                                clustering[i]])
        minsep <- floor(n * sepprob)
      }
      else {
        dmatnn <- dmat[clustering <= cwn, clustering <=
                         cwn]
        clusteringnn <- clustering[clustering <= cwn]
        for (i in 1:(n - noisen))
          psep[i] <- min(dmatnn[i,
                                clusteringnn != clusteringnn[i]])
        minsep <- floor((n - noisen) * sepprob)
      }
      sindex <- mean(sort(psep)[1:minsep])
    }
    if (!aggregateonly)
      out <-
      list(
        n = n,
        cluster.number = cn,
        cluster.size = cluster.size,
        min.cluster.size = min(cluster.size[1:cwn]),
        noisen = noisen,
        diameter = diameter,
        average.distance = average.distance,
        median.distance = median.distance,
        separation = separation,
        average.toother = average.toother,
        separation.matrix = separation.matrix,
        ave.between.matrix = ave.between.matrix,
        average.between = average.between,
        average.within = average.within,
        n.between = nbetween,
        n.within = nwithin,
        max.diameter = max(diameter[1:cwn],
                           na.rm = TRUE),
        min.separation = sepwithnoise *
          min(separation) + (!sepwithnoise) * min(separation[1:cwn]),
        within.cluster.ss = within.cluster.ss,
        clus.avg.silwidths = clus.avg.widths,
        avg.silwidth = avg.width,
        g2 = g2,
        g3 = g3,
        pearsongamma = pearsongamma,
        dunn = dunn,
        dunn2 = dunn2,
        entropy = h1,
        wb.ratio = average.within / average.between,
        ch = ch,
        cwidegap = cwidegap,
        widestgap = widestgap,
        sindex = sindex,
        corrected.rand = corrected.rand,
        vi = vi
      )
    else
      out <-
      list(
        n = n,
        cluster.number = cn,
        min.cluster.size = min(cluster.size[1:cwn]),
        noisen = noisen,
        average.between = average.between,
        average.within = average.within,
        max.diameter = max(diameter[1:cwn],
                           na.rm = TRUE),
        min.separation = sepwithnoise *
          min(separation) + (!sepwithnoise) * min(separation[1:cwn]),
        ave.within.cluster.ss = within.cluster.ss / (n - noisen),
        avg.silwidth = avg.width,
        g2 = g2,
        g3 = g3,
        pearsongamma = pearsongamma,
        dunn = dunn,
        dunn2 = dunn2,
        entropy = h1,
        wb.ratio = average.within / average.between,
        ch = ch,
        widestgap = widestgap,
        sindex = sindex,
        corrected.rand = corrected.rand,
        vi = vi
      )
  }
  out
}


#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#************************************FUNCTION Cluster.analysis**********************************
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*

library(clusterCrit)
cl.val <- function(method = "kmeans",x = df.comp,maxK)
{
  library(clusterCrit)
  hc.pca <- PCA(df.comp,ncp = 2,graph = FALSE,scale.unit = TRUE)
  Clustinfo <- list()
  
  allInd <- c("c_index","davies_bouldin","ray_turi","sd_dis" ,"s_dbw","xie_beni",
            #Max Indices
            "calinski_harabasz","dunn","gamma","ratkowsky_lance","silhouette","tau","wemmert_gancarski")

  for (i in 2:maxK) {
    if (method == "kmeans") {
      clust.res <- kmeans(x ,centers = i,iter.max = 100)
      cluster <- as.integer(clust.res$cluster)
      ls <- intCriteria(as.matrix(x), cluster, allInd)
      val.res <-data.frame(val = matrix(unlist(ls), ncol = length(allInd), byrow = TRUE))
      colnames(val.res) <- allInd
      val.res$k <- i
      Clustinfo[[i - 1]] <- val.res
      #---------------------------------------------------
    } else if (method == "HCPC") {
      clust.res <- HCPC(hc.pca, graph = FALSE, nb.clust = i)
      ls <- intCriteria(as.matrix(x), as.integer(clust.res$data.clust$clust), allInd)
      val.res <-data.frame(val = matrix(unlist(ls), ncol = length(allInd), byrow = TRUE))
      colnames(val.res) <- allInd
      val.res$k <- i
      Clustinfo[[i - 1]] <- val.res
      #---------------------------------------------------
    } else if (method == "Fuzzy") {
      clust.res <- fanny(x, k = i)
      ls <- intCriteria(as.matrix(x), as.integer(clust.res$clustering), allInd)
      val.res <-data.frame(val = matrix(unlist(ls), ncol = length(allInd), byrow = TRUE))
      colnames(val.res) <- allInd
      val.res$k <- i
      Clustinfo[[i - 1]] <- val.res
      #---------------------------------------------------
    } else if (method == "EM") {
      clust.res = Mclust(x, i)
      ls <- intCriteria(as.matrix(x), as.integer(clust.res$classification), allInd)
      val.res <-data.frame(val = matrix(unlist(ls), ncol = length(allInd), byrow = TRUE))
      colnames(val.res) <- allInd
      val.res$k <- i
      Clustinfo[[i - 1]] <- val.res
    }
  } # END OF FOR
  Clust_info = data.frame(val = matrix(unlist(Clustinfo), ncol = length(allInd)+1, byrow = TRUE))
  colnames(Clust_info) <- c(allInd,"K")
  return(Clust_info)
  
} # END OF FUNCTION Cluster.analysis

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*

min.criteria <- c("c_index","davies_bouldin","ray_turi","sd_dis" ,"s_dbw","xie_beni")
max.criteria <- c("dunn","gamma","ratkowsky_lance","silhouette","tau","wemmert_gancarski")

#--------------------------------- Application --------------------------------
#------------------------------------------------------------------------------
dim(Lockhart_Data)
dim(TUFFP_Author_Data)
dim(TUFFP_comp)
dim(TUFFP_Data)
dim(Lockhart_BriFeli)

df.comp <- scale(Lockhart_BriFeli[,2:3])
dim(df.comp)
colours <- c("red", "hotpink",  "green", "blue","purple","orange")

clvals_km <- cl.val(method = "kmeans", x=df.comp,maxK = 10)
clvals_fuzzy <- cl.val(method = "Fuzzy", x=df.comp,maxK = 10)
# clvals_em <- cl.val(method = "EM", x=df.comp,maxK = 10)
 clvals_hc <- cl.val(method = "HCPC", x=df.comp,maxK = 10)
 
 df.comp.lm <- scale(Lockhart_Data[,2:3])
 clvals_km.lm <- cl.val(method = "kmeans", x=df.comp.lm,maxK = 10)
 clvals_hc.lm <- cl.val(method = "HCPC", x=df.comp.lm,maxK = 10)
 
#********************************************************************************

indices <-as.data.frame(cbind(clvals_fuzzy[,min.criteria],clvals_fuzzy[,max.criteria]))
round(indices, digits = 4)
write.csv("C:/Users/sara/Desktop/phd/ValidationRes/bf_LM_FuzzyCriteria.csv",x = round(indices, digits = 5))


#****************************************************************************************************
TUFFP_KM_Criteria
TUFFP_hc_Criteria
TUFFP_hc_lmCriteria
cl.values <- clvals_fuzzy[,c("K",min.criteria)]
colnames(cl.values)

df.val <- data.frame()
for(i in 1:6) {
  temp.val <-as.data.frame(cbind(k = 2:10,index = rep(colnames(cl.values)[i + 1], 9),value = cl.values[, i + 1]),stringsAsFactors = FALSE)
  df.val <- rbind(df.val, temp.val)
}
df.val$type <- "min"
invalid <- which(df.val$value =="NaN" | df.val$value == "Inf" | df.val$value == "-Inf")
invalid
mydata <- df.val[-invalid,]
dataplot <- mydata[order(as.integer(mydata$k),decreasing = FALSE),]
sapply(dataplot, mode)
dataplot$k <- as.numeric(as.character(dataplot$k))
dataplot$value <- as.numeric(as.character(dataplot$value))
which.max(dataplot$value)
min(dataplot$value)
max(dataplot$value)

scaleFUN <- function(x) sprintf("%.2f", x)
sciFun <- function(x) formatC(x, format = "e", digits = 1)

bf.minvalP.km.lm <- ggplot(data = dataplot, aes(group  = index)) +
  labs(x = "Number of Clusters", y = "Min Validity Index (log trans)",title = "K-means Method") +
  geom_point(aes(x = k, y = value, colour = index)) +
  geom_line(mapping = aes(x = k, y = value, colour = index)) +
  scale_x_continuous(breaks = seq(1, 10, 1), limits = c(2, 10)) +
  scale_y_continuous(trans = "log10", labels=sciFun)+
  scale_colour_manual(values = colours)+
  #scale_color_manual(labels = c("EM", "Fuzyy","Hierarchical","K-means"), values = colours)
  theme_light() + theme(legend.title = element_blank(),legend.position = "bottom",
    legend.key.size = unit(0.2, "cm"),legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"))
minvalP.km.lm
minvalP.hc.lm
#**************************************max criteria****************************
#******************************************************************************

cl.values.max <- clvals_km.lm[,c("K",max.criteria)]
colnames(cl.values.max)

df.val.max <- data.frame()
for(i in 1:6) {
  temp.val <-as.data.frame(cbind(k = 2:10,index = rep(colnames(cl.values.max)[i + 1], 9),value = cl.values.max[, i + 1]),stringsAsFactors = FALSE)
  df.val.max <- rbind(df.val.max, temp.val)
}

df.val.max$type <- "max"
invalid <- which(df.val.max$value =="NaN" | df.val.max$value == "Inf" | df.val.max$value == "-Inf")
invalid
mydata <- df.val.max#[-invalid,]
dataplot.max <- mydata[order(as.integer(mydata$k),decreasing = FALSE),]
sapply(dataplot.max, mode)
dataplot.max$k <- as.numeric(as.character(dataplot.max$k))
dataplot.max$value <- as.numeric(as.character(dataplot.max$value))
max(dataplot.max$value)
min(dataplot.max$value)

maxvalP.km.lm <- ggplot(data = dataplot.max, aes(group  = index)) +
  labs(x = "Number of Clusters", y = "Max Validity Index",title = "LM with K-means") +
  geom_point(aes(x = k, y = value, colour = index)) +
  geom_line(mapping = aes(x = k, y = value, colour = index)) +
  scale_x_continuous(breaks = seq(1, 10, 1), limits = c(2, 10)) +
  scale_y_continuous(trans = "log2",labels=scaleFUN)+
  scale_colour_manual(values = colours)+
  theme_light()+ theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"))

# maxvalP.km$theme = theme_light()+ theme(
#   legend.title = element_blank(),
#   legend.position = "bottom",
#   legend.key.size = unit(0.2, "cm"),
#   legend.key.height = unit(0.2, "cm"),
#   legend.key.width = unit(0.2, "cm")
# )
# maxvalP.km$labels$title = "Dimensional Data with K-means"
# minvalP.km$labels$title = "Dimensional Data with K-means"
maxvalP.km

# minvalP.hc$labels$title = "Dimensional Data with Hierarchical"
# maxvalP.hc$labels$title = "Dimensional Data with Hierarchical"

#**************************************************************
# minvalP.km.lm$labels$title = "LM with K-means"
# maxvalP.km.lm$labels$title = "LM with K-means"
# 
# minvalP.hc.lm$labels$title = "LM with Hierarchical"
# maxvalP.hc.lm$labels$title = "LM with Hierarchical"

#*****************************************************************************
library(gridExtra)
library(ggpubr)
mylegend <- g_legend(minvalP.km)
mylegend2 <- g_legend(maxvalP.km)

png("C:/Users/sara/Desktop/phd/chapter5/Dless_Validation/TUFFP_MinMax_hckm.png", units="in", 
      width=8, height=8,res = 500)

ggarrange(minvalP.km, maxvalP.km, minvalP.hc, maxvalP.hc, labels = c("A", "B","C","D"),
          ncol = 2, nrow = 2)
dev.off()

#grid.arrange(maxvalP_exp, maxvalP_log10,maxvalP_log1p,maxvalP_log2,maxvalP_logit,maxvalP_sqrt)
# 
# grid.arrange(arrangeGrob(minvalP ,maxvalP,widths = c(1, 1)),ncol=1,
#   nrow = 2,heights = c(10, 1))




# library("cowplot")
# plot_grid(minvalP, maxvalP , rel_widths = c(1,1),
#           labels = c("1", "2"),
#           ncol = 2, nrow = 1)
# 
#************************************************************

min.index.km <- apply(clvals[,min.criteria],MARGIN = 2,FUN = which.min)
max.index.km <- apply(clvals[,max.criteria],MARGIN = 2,FUN=which.max)

opt.k.km <- as.data.frame(c(unlist(min.index.km)+1,max.index.km+1))
colnames(opt.k.km) <- "IndexVal"
opt.k.km <- t(opt.k.km)
opt.k.km

#-----------------------------------------------------------
clval_hc <- cl.val(method = "HCPC", x=df.comp,maxK = 10)
# saveRDS(clval_hc, file = "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/LM/clval_hc.rds")
# clval_hc <- readRDS(file = "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/LM/clval_hc.rds")
clval_hc[,min.criteria]
clval_hc[,max.criteria]
min.index.hc <- apply(clval_hc[,min.criteria],MARGIN = 2,FUN = which.min)
max.index.hc <- apply(clval_hc[,max.criteria],MARGIN = 2,FUN=which.max)
opt.k.hc <- c(unlist(min.index.hc)+1,max.index.hc+1)
opt.k.hc
#-----------------------------------------------------------

clval_fuzzy <- cl.val(method = "Fuzzy", x=df.comp,maxK = 10)
# saveRDS(clval_fuzzy, file = "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/LM/clval_fuzzy.rds")
# clval_fuzzy <- readRDS(file = "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/LM/clval_fuzzy.rds")

min.index.fz <- unlist(apply(clval_fuzzy[,min.criteria],MARGIN = 2,FUN = which.min))
max.index.fz <- unlist(apply(clval_fuzzy[,max.criteria],MARGIN = 2,FUN=which.max))
opt.k.fz <- c(min.index.fz+1,max.index.fz+1)
opt.k.fz

#-----------------------------------------------------------

clval_EM <- cl.val(method = "EM", x=df.comp,maxK = 10)
# saveRDS(clval_EM, file = "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/LM/clval_EM.rds")
# clval_EM <- readRDS(file = "C:/Users/sam805/Desktop/My PhD Research/ABM_2020SP/LM/clval_EM.rds")

min.index.em <- unlist(apply(clval_EM[,min.criteria],MARGIN = 2,FUN = which.min))
max.index.em <- unlist(apply(clval_EM[,max.criteria],MARGIN = 2,FUN=which.max))
opt.k.em <- c(min.index.em+1,max.index.fz+1)
opt.k.em

#**********************************************************************************************************

#*******************************************************************************
ggplot(data = as.data.frame(ind.val), mapping = aes(x = k, y = davies_bouldin)) +
  geom_line() + 
  labs(x = "K", y = "Indices")

ggplot(ind.val, aes(x = k))+
  geom_line(aes(y = davies_bouldin), color = "darkred") + 
  geom_line(aes(y = s_dbw), color="blue") +
  theme_minimal()


idx <- bestCriterion(vals,"c_index")
cat("Best index value is",vals[idx],"\n")

# Perform the kmeans algorithm
cl <- kmeans(df.comp, centers = 10,iter.max = 30)

# Compute all the internal indices
intCriteria(df.comp,cl$cluster,"all")

# Compute some of them
intCriteria(df.comp,cl$cluster,c("C_index","Calinski_Harabasz","Dunn"))

# The names are case insensitive and can be abbreviated
intCriteria(df.comp,cl$cluster,c("det","cal","dav"))

#******************************************* The Silhouette Method *********************************************

#FUNcluster = c("kmeans", "pam", "clara", "fanny", "hclust", "agnes", "diana")
fviz_nbclust(df.comp, kmeans, method = "silhouette", k.max = 10) + theme_minimal() + ggtitle("The Silhouette Plot")
fviz_nbclust(df.comp, clara, method = "silhouette", k.max = 10) + theme_minimal() + ggtitle("The Silhouette Plot")
fviz_nbclust(df.comp, pam, method = "silhouette", k.max = 10) + theme_minimal() + ggtitle("The Silhouette Plot")
fviz_nbclust(df.comp, hcut, method = "silhouette", k.max = 10) + theme_minimal() + ggtitle("The Silhouette Plot")

fviz_nbclust(df.comp, hcut, method = "wss", k.max = 10) + theme_minimal() + ggtitle("The Silhouette Plot")

#**************************************************************************************************************
#**************************************************************************************************************
cluster.validation <- function(method = "kmeans",x = df.comp,maxK = 10,full.data)
{
  library(clv)
  library(factoextra)
  library(FactoMineR)
  library("fpc")
  intraclust = c("complete","average","centroid")
  interclust = c("single", "complete", "average","centroid", "aveToCent", "hausdorff")
  sdbw.val <- data.frame(K = integer(), sdbwval = integer())
  sdbw.val[1:9,   1] <- seq(2, 10, 1)
  sd.val <- data.frame(K = integer(), sdval = integer())
  sd.val[1:9, 1] <- seq(2, 10, 1)
  dunn1.val <- list()
  davies1.val <- list()
  dunn2.val <- list()
  davies2.val <- list()
  hc.pca <- PCA(full.data, ncp = 2, graph = FALSE,scale.unit = TRUE)
  Clustinfo <- list()
  
  for (i in 2:maxK) {
    if (method == "kmeans") {
      clust.res <- kmeans(x ,centers = i,iter.max = 100)
      clusts<- as.integer(clust.res$cluster)
      cs = cluster.stats(d = dist(x), clustering = clust.res$cluster)
      gc()
      Clustinfo[[i-1]] = unlist(cs[c("cluster.number","within.cluster.ss","avg.silwidth","dunn","dunn2","sindex")])
      
    } else if (method == "HCPC") {
      clust.res <- HCPC(hc.pca, graph = FALSE, nb.clust = i)
      clusts <- as.integer(clust.res$data.clust$clust)
      cs = my.cluster.stats(d = dist(x), clustering = clust.res$data.clust$clust)
      gc()
      Clustinfo[[i-1]] = unlist(cs[c("cluster.number","within.cluster.ss","avg.silwidth","dunn","dunn2","sindex")])
      
    } else if (method == "Fuzzy") {
      clust.res <- fanny(x, k = i)
      clusts <- as.integer(clust.res$clustering)
      cs = cluster.stats(d = dist(x), clustering = clust.res$clustering)
      gc()
      Clustinfo[[i-1]] = unlist(cs[c("cluster.number","within.cluster.ss","avg.silwidth","dunn","dunn2","sindex")])
      
    } else if (method == "EM") {
      clust.res = Mclust(df.comp, i)
      clusts <- as.integer(clust.res$classification)
      cs = cluster.stats(d = dist(x), clustering = clust.res$classification)
      gc()
      Clustinfo[[i-1]] = unlist(cs[c("cluster.number","within.cluster.ss","avg.silwidth","dunn","dunn2","sindex")])
      
    }
    
    scatt <- clv.Scatt(x, clusts)
    dens.bw <- clv.DensBw(x, clusts, scatt)
    dis <- clv.Dis(scatt$cluster.center)
    # SD
    sd.val[i - 1, 2] <-
      clv.SD(scatt$Scatt, dis, alfa = i) # alfa is equal to number of clusters
    # S_Dbw
    sdbw.val[i - 1, 2] <- dens.bw
    # Dunns & Davies-Bouldin
    cls.scatt <- cls.scatt.data(x, clusts, dist = "manhattan")
    dunn1.val[[i - 1]] <- clv.Dunn(cls.scatt, intraclust, interclust)
    dunn2.val[[i - 1]] <- Dunn(x, clusts)
    davies1.val[[i - 1]] <- clv.Davies.Bouldin(cls.scatt, intraclust, interclust)
    davies2.val[[i - 1]] <- Davies.Bouldin(x, clusts)
    
  } # END OF FOR
  
  return(list(sdbw = sdbw.val,sd=sd.val,dunn1 =dunn1.val,davies1= davies1.val,dunn2= dunn2.val ,davies2=davies2.val,Clustinfo ))
  
} # END OF FUNCTION Cluster.analysis

#***************************************************************************

cl_kmeans <- cluster.validation(method = "kmeans", x=df.comp,maxK = 10,full.data = bfdata)
cl_hc <- cluster.validation(method = "HCPC", x=df.comp,maxK = 10,full.data = bfdata)
cl_fuzzy <- cluster.validation(method = "Fuzzy", x=df.comp,maxK = 10,full.data = bfdata)
cl_EM <- cluster.validation(method = "EM", x=df.comp,maxK = 10,full.data = bfdata)
cl_kmeans$sdbw 
cl_kmeans$sdbw[which.min(cl_kmeans$sdbw$sdbwval),] 
cl_kmeans$sd[which.min(cl_kmeans$sd$sdval),] 

# Dunn's: Maximum is the best one
d.index <- lapply(cl_kmeans$dunn1, FUN = max)
dunns.index <- data.frame(K= 2:10,Dunn  = unlist(d.index))
dunns.index[which.max(dunns.index$Dunn),] 
#***************************************

# Davies's: Minimum is the best one
dav.ind <- lapply(cl_kmeans$davies1, FUN = min)
davies.ind <- data.frame(K= 2:10,Davies  = unlist(dav.ind))
davies.ind[which.min(davies.ind$Davies),] 
#***************************************

d2.index <- lapply(cl_kmeans$dunn2, FUN = max)
dunns2.index <- data.frame(K= 2:10,Dunn  = unlist(d2.index))
dunns2.index[which.max(dunns2.index$Dunn),] 
#***************************************
# Davies
dav2.ind <- lapply(cl_kmeans$davies2, FUN = min)
davies2.ind <- data.frame(K= 2:10,Davies  = unlist(dav2.ind))
davies2.ind[which.min(davies2.ind$Davies),] 


#**************************************************************************************
#***************************** Elbow Method *******************************************
#***************************************************************************************
# The Elbow Curve method is helpful because it shows how increasing the number of the clusters contribute separating the clusters 
# in a meaningful way, not in a marginal way. The bend indicates that additional clusters beyond the third have little value 
# The Elbow method is fairly clear, if not a naïve solution based on intra-cluster variance. 
# The gap statistic is more sophisticated method to deal with data that has a distribution with no obvious clustering 
# (can find the correct number of k for globular, Gaussian-distributed, mildly disjoint data distributions).
# function to compute total within-cluster sum of squares

fviz_nbclust(df.comp, kmeans, method = "wss", k.max = 10) + theme_minimal() + ggtitle("the Elbow Method")

#**************************************************************************************
#***************************** The Gap Statistic **************************************
#**************************************************************************************

# The gap statistic compares the total within intra-cluster variation for different values of k with their expected values under null reference
# distribution of the data. The estimate of the optimal clusters will be value that maximize the gap statistic 
# (i.e., that yields the largest gap statistic). This means that the clustering structure is far away from the random uniform 
# distribution of points.

gap_stat <- clusGap(df.comp, FUN = kmeans, nstart = 30, K.max = 10, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

#**************************************************************************************

#**************************************************************************************
#*************************** The Sum of Squares Method ********************************
#**************************************************************************************

# Another clustering validation method would be to choose the optimal number of cluster by minimizing the within-cluster 
# sum of squares (a measure of how tight each cluster is) and maximizing the between-cluster sum of squares 
# (a measure of how seperated each cluster is from the others).

kmean_calc <- function(df, ...) {
  kmeans(df, scaled = ..., nstart = 30)
}
km2 <- kmean_calc(df.comp, 2);km3 <- kmean_calc(df.comp, 3);km4 <- kmeans(df.comp, 4);km5 <- kmeans(df.comp, 5);km6 <- kmeans(df.comp, 6)
km7 <- kmeans(df.comp, 7);km8 <- kmeans(df.comp, 8);km9 <- kmeans(df.comp, 9);km10 <- kmeans(df.comp, 10);km11 <- kmeans(df.comp, 11)
km12 <- kmeans(df.comp, 12);km13 <- kmeans(df.comp, 13);km14 <- kmeans(df.comp, 14);km15 <- kmeans(df.comp, 15);km16 <- kmeans(df.comp, 16)
km17 <- kmeans(df.comp, 17);km18 <- kmeans(df.comp, 18);km19 <- kmeans(df.comp, 19);km20 <- kmeans(df.comp, 20);km21 <- kmeans(df.comp, 21)
km22 <- kmeans(df.comp, 22);km23 <- kmeans(df.comp, 23);km24 <- kmeans(df.comp, 24);km25 <- kmeans(df.comp, 25);km26 <- kmeans(df.comp, 26)
km27 <- kmeans(df.comp, 27);km28 <- kmeans(df.comp, 28)


ssc <- data.frame(
  kmeans = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28),
  within_ss = c(mean(km2$withinss),mean(km3$withinss),mean(km4$withinss),mean(km5$withinss),mean(km6$withinss),
    mean(km7$withinss),mean(km8$withinss),mean(km9$withinss),mean(km10$withinss),mean(km11$withinss),mean(km12$withinss),
    mean(km13$withinss),mean(km14$withinss),mean(km15$withinss),mean(km16$withinss),mean(km17$withinss),mean(km18$withinss),
    mean(km19$withinss),mean(km20$withinss),mean(km21$withinss),mean(km22$withinss),mean(km23$withinss),
    mean(km24$withinss),mean(km25$withinss),mean(km26$withinss),mean(km27$withinss),mean(km28$withinss)
  ),
  between_ss = c(km2$betweenss,km3$betweenss,km4$betweenss,km5$betweenss,km6$betweenss,km7$betweenss,
    km8$betweenss,km9$betweenss,km10$betweenss,km11$betweenss,km12$betweenss,km13$betweenss,km14$betweenss,km15$betweenss,
    km16$betweenss,km17$betweenss,km18$betweenss,km19$betweenss,km20$betweenss,km21$betweenss,km22$betweenss,km23$betweenss,
    km24$betweenss,km25$betweenss,km26$betweenss,km27$betweenss,km28$betweenss
  ))




#**************************************************************************************
#******************************* The NbClust ******************************************
#**************************************************************************************


# The NbClust package provides 30 indices for determining the relevant number of clusters and proposes to users the best clustering
# scheme from the different results obtained by varying all combinations of number of clusters, distance measures, and clustering methods.

res.nbclust <- NbClust(df.comp, distance = "euclidean",min.nc = 2, max.nc = 30, method = "complete", index ="all")
factoextra::fviz_nbclust(res.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")

#**************************************************************************************
#********************* Choosing the appropriate algorithm *****************************
#**************************************************************************************
library(clValid)
intern <- clValid(df.comp, nClust = 2:24, clMethods = c("hierarchical","kmeans","pam"), validation = "internal")# Summary
summary(intern) %>% kable() %>% kable_styling()

# Compute dissimilarity matrix with euclidean distances
d <- dist(df.comp, method = "euclidean")# Hierarchical clustering using Ward's method
res.hc <- hclust(d, method = "ward.D2" )# Cut tree into 5 groups
grp <- cutree(res.hc, k = )# Visualize
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = , border = 2:5) # add rectangle

# Execution of k-means with k=5
final <- kmeans(df.comp, , nstart = 30)
fviz_cluster(final, data = mammals_scaled) + theme_minimal() + ggtitle("k = 5")

#*******************************************************
#install.packages(c("factoextra", "clustertend"))
library(clustertend)
library(factoextra)
fviz_pca_ind(pc, title = "PCA",
             # habillage = iris$Species,
             palette = "jco",
             geom = "point", ggtheme = theme_classic(),
             legend = "bottom")

#*********Hierarchical clustering*********************
#1


#Elbow method
fviz_nbclust(df.comp, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(df.comp, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy.
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(df.comp, kmeans, nstart = 25, method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

nb <- NbClust(df.comp, distance = "euclidean", min.nc = 2,max.nc = 40, method = "kmeans")
library("factoextra")
fviz_nbclust(nb)  



library(ggplot2)
dataplot %>% 
  group_by(index) %>%
  summarise(min = min(value),min2 = sort(dataplot$value)[2]) -> me.2

left_join(dataplot, me.2) %>%
  mutate(color = value == min | value == min2) %>%
  filter(color == TRUE) -> me.3

ggplot(data = dataplot, aes(group  = index,x = k, y = value)) +
  labs(x = "Number of Clusters", y = "Min Validity Index") +
  geom_point(#data=me.3,
    aes(x = k, y = value, colour = index)) +
  facet_wrap(~index, ncol=3, scales = "free_y")+
  geom_line(mapping = aes(x = k, y = value, colour = index)) +
  scale_x_continuous(breaks = seq(1, 10, 1), limits = c(2, 10)) +
  scale_colour_manual(values = colours)+
  theme_light() + theme(
    legend.title = element_blank(),
    legend.position = "none")

