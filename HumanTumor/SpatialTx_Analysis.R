##### Make changes to SPACE code ##############################################
setwd("C:/Users/schromec/Desktop/SPACE/SPACE/R")
devtools::document()
devtools::install()
#restart R session
setwd("C:/Users/schromec/Desktop/SPACE/SpatialTx")
library(SPACE)

# Set working directory
setwd("C:/Users/schromec/Desktop/SPACE/SpatialTx")

# Load data from Visium
d <- Matrix::readMM("matrix.mtx")                                               # Gene x spot matrix  
s <- readr::read_tsv("barcodes.tsv", col_names=F)                               # Spot barcode IDs
colnames(s) <- "ID"
g <- readr::read_tsv("features.tsv", col_names=F)                               # Gene IDs
g <- g[,1:2]
colnames(g) <- c("ID", "Name")
p <- readr::read_csv("all_barcodes.csv", col_names=F)                           # Spot spatial coordinates
colnames(p) <- c("ID", "U1", "U2", "U3", "X", "Y")

# Match spatial coordinates to spots.
s$X <- rep(0, nrow(s))
s$Y <- rep(0, nrow(s))
idx <- unlist(lapply(s$ID, function(x) {which(p[,1] == x)}))                    # Find used spots amid all spots
s$X <- p$X[idx]                                                                 # Record spatial coors of used spots
s$Y <- p$Y[idx]
s$Z <- rep(1, nrow(s))
plot(p$Y, max(p$X) - p$X, pch=1)                                                # Confirm spot coors are correct
points(s$Y, max(p$X) - s$X, pch=16)
remove(p, idx)

# Filter genes to only those that appear in the Visium spatial analyses (17579 instead of 17943)
setwd("C:/Users/schromec/Desktop/SPACE/SpatialTx/Visium_Analysis")
m <- readr::read_csv("MoransI.csv")
colnames(m) <- c("ID", "Name", "I", "P", "P_adj", "Total_Gene_Count", "Median_Normalized_Avg_Count", "Total_Spots")
idx <- unlist(lapply(g$ID, function(x) {sum(m$ID == x)})) == 1                  # Genes that appear in spatial results
g <- g[idx, ]
d <- d[idx, ]
remove(idx)
setwd("C:/Users/schromec/Desktop/SPACE/SpatialTx")

# Filter genes to only those that appear in 50+% of spots
m <- m[match(g$ID, m$ID), ]                                                     # Order spat results same as genes
idx <- m$Total_Spots >= (nrow(s)/2)                                             # Genes in >= 50% of spots
m <- m[idx, ]
g <- g[idx, ]
d <- d[idx, ]
remove(idx)

# Transform gene x spot matrix into an object table, with only one type of object, i.e. a generic spot
TAB <- cbind(s[,2:4], rep(1, nrow(s)), t(as.matrix(d)))
colnames(TAB) <- c("X","Y","Z","Object",g$Name)

# Create a link table linking all transcripts to the lone object type
LNK <- matrix(rep(1,nrow(g)), nrow=1)
colnames(LNK) <- paste("S1.",1:nrow(g),sep="")
rownames(LNK) <- "O1.1"

# Remove old data structures now that they are redundant
remove(d,g,s)

# Determine a radius for censusing
FNN::knn.dist(TAB[,c("X","Y","Z")], k=20)                                       # Spots are 144um apart, hexagonally
RAD <- 300                                                                      # 300um = spot & its 18 neighbors
NUM <- nrow(TAB)                                                                # Just use every spot as a seed point

# Census the data
system.time(CEN <- SPACE::census_table(TAB, list(O1 = LNK), radii = RAD, sample_size = NUM))
PLS <- CEN[[2]]
CEN <- CEN[[1]]
#save(CEN, file = "SpatTx_CEN_150.Rdata")
#save(PLS, file = "SpatTx_PLS_150.Rdata")

# Measure mutual information
system.time(CMI <- SPACE::measure_cisMI(CEN, PLS, depth=1, radii=150, bootstraps = 10, not = "O1.1", max_bins = 10))
#save(CMI, file = "SpatTx_CMI_150.Rdata")

# Panel A: Compare 1D mutual info to Moran's I
# Export as .pdf at 5" x 5" then convert to png.
C1D <- data.frame(ID = m$ID, Name = m$Name, space_i = CMI[[1]]$CisMI, moran_i = m$I)
C1D <- C1D[C1D$moran_i > 0, ]
C1D$space_si <- (C1D$space_i - mean(C1D$space_i)) / sd(C1D$space_i)
C1D$moran_li <- log(C1D$moran_i)
C1D$moran_sli <- (C1D$moran_li - mean(C1D$moran_li)) / sd(C1D$moran_li)
mod <- lm(C1D$moran_sli ~ C1D$space_si)
summary(mod)
ggplot2::ggplot(data=C1D, ggplot2::aes(x=space_si, y=moran_sli)) +
  ggplot2::geom_point(pch=16, fill="black", alpha=0.25, size=3) + 
  ggplot2::geom_smooth(method=lm, color="red", fill="red", alpha=0.25, se=TRUE, level=0.95) +
  ggplot2::xlim(c(-2.75,2.75)) + ggplot2::ylim(c(-3,3)) +
  ggplot2::ylab("Standardized Log Moran's I") + ggplot2::xlab("Standardized CisMI") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
remove(C1D, mod)

# Panel B: Load k-means clustering and differential expression analyses from Visium and plot
# Export as .pdf at 6" x 7", convert to png, and crop.
# Export palette as metafile at 700 x 225
setwd("C:/Users/schromec/Desktop/SPACE/SpatialTx/Visium_Analysis/diffexp/kmeans_8_clusters")
KDE <- read.csv("differential_expression.csv", header = T)
setwd("C:/Users/schromec/Desktop/SPACE/SpatialTx/Visium_Analysis/clustering/kmeans_8_clusters")
KCL <- read.csv("clusters.csv", header = T)
setwd("C:/Users/schromec/Desktop/SPACE/SpatialTx")
KDE <- KDE[match(m$ID, KDE$Feature.ID), ]
KCL <- cbind(KCL, TAB[,c("X","Y","Z")])
KPL <- c("#E31A1C","dodgerblue2","#6A3D9A","#FF7F00","green4","gold1","blue4","palegreen2")
SPACE::plot_palette(KPL)
plot(KCL$Y, max(KCL$X) - KCL$X, pch=16, col=KPL[KCL$Cluster])

# Find genes that best distinguish the Visium k-means clusters and only consider those
pcols <- grepl("p.value", colnames(KDE))
table(rowSums(KDE[,pcols] < 0.05))
vars <- which(rowSums(KDE[,pcols] < 0.05) >= 4)
TABK <- TAB[,c(1:4,vars+4)]
LNKK <- LNK[,vars,drop=F]
colnames(LNKK) <- paste("S1.",1:length(vars),sep="")

# Recensus the data with a focus on just the top distinguishing genes
system.time(CENK <- SPACE::census_table(TABK, list(O1 = LNKK), radii = 150, sample_size = NUM))
PLSK <- CENK[[2]]
CENK <- CENK[[1]]
#save(CENK, file = "SpatTx_CEN_K.Rdata")
#save(PLSK, file = "SpatTx_PLS_K.Rdata")

# Calculate cisMI up to an ensemble depth of 3 for the transcripts that most contribute to k-means
system.time(CMIK <- SPACE::measure_cisMI(CENK, PLSK, depth=3, radii=150, bootstraps = 100, not = "O1.1"))
#save(CMIK, file = "SpatTx_CMI_K.Rdata")

# Make and plot a palette for this select group of transcripts
# Export as metafile at 150 x 700
PALK <- SPACE::make_palette(length(vars))
SPACE::plot_palette(PALK, horizontal=T)
#save(PALK, file="SpatTx_PAL_K.Rdata")

# Panel C: Plot ranked ensembles for this select group of transcripts
# Export as .pdf at 7.5" x 7.5" and convert to .png
CMIKP <- SPACE::plot_MI_rank(CMIK, radius = 150, depth = 1:3, col_pals = list(S1 = PALK), mi_thr = 0.1, 
                              p_thr = exp(-392)) # p_thr just limits to top 25 ensembles
CMIKP[[1]][which.max(CMIKP[[1]]$Zscore), ]
CMIKP[[1]][order(abs(CMIKP[[1]]$CisMI), decreasing=T)[1:40], ]

# Panel D: Investigate the 3-var ensemble with highest cisMI
# Export covariation plot as .pdf at 7" x 5", edit, and export as .png
colnames(TABK)[4+c(3,7,13)]
PTRN <- SPACE::learn_pattern(CENK, c("O1.1_S1.3", "O1.1_S1.7", "O1.1_S1.13"), 150, 
                              list(O1="#000000", S1=PALK), patch_list=PLSK)

# Panel E: Map features of the 3-var ensemble back to the Visium spots
# Export map as .pdf at 6" x 7", convert to png, and crop.
# Export palette as metafile at 225 x 220
PTRNS <- PTRN[,c("X","Y","V")]
PTRNS <- tidyr::pivot_wider(PTRNS, names_from = V, values_from = Y)
TABKS <- TABK[,c(1:4,4+c(3,7,13))]
MAP <- FNN::knnx.index(PTRNS[,(ncol(PTRNS)-2):ncol(PTRNS)], TABKS[,(ncol(TABKS)-2):ncol(TABKS)], k=1)
MAP <- MAP / nrow(TABKS)
MAP <- unlist(lapply(MAP, function(x) {sum(x >= c(0, 0.87))}))
plot(KCL$Y, max(KCL$X) - KCL$X, pch=16, col=c("gray", "green4")[MAP])
SPACE::plot_palette(c("gray", "green4"))

# Panel F: Correspondence of covariation-feature spots with k-means spots
# Export as .pdf at 5.5" x 5", convert to png, and crop
BP <- data.frame(ME = 1:8, Count = unname(table(factor(KCL$Cluster[MAP == 2], levels=c(1:8)))))
BP$Count.Var1 <- NULL
colnames(BP)[2] <- "Count"
ggplot2::ggplot(data=BP, ggplot2::aes(x=as.factor(ME), y=Count, fill=as.factor(ME))) +
  ggplot2::geom_bar(stat="identity") + ggplot2::scale_fill_manual(values = KPL) +
  ggplot2::ylab("Feature 1 Spot Count") + ggplot2::xlab("Original ME") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
tab <- matrix(c(sum(MAP==2 & KCL$Cluster == 5), sum(MAP!=2 & KCL$Cluster == 5),
                sum(MAP==2 & KCL$Cluster != 5), sum(MAP!=2 & KCL$Cluster != 5)), nrow=2)
psych::cohen.kappa(tab)
remove(BP, tab)

# Panels S A-C: correlation among all pairs of variables
# Export each as .pdf at 5" x 5"
MD <- as.data.frame(scale(TABKS[,(ncol(TABKS)-2):ncol(TABKS)]))
summary(lm(GPX2 ~ MGAM2, data=MD))
ggplot2::ggplot(MD, ggplot2::aes(x=MGAM2, y=GPX2)) +
  ggplot2::geom_point(pch=16, fill="black", alpha=0.25, size=3) + 
  ggplot2::geom_smooth(method=lm, color="red", fill="red", alpha=0.25, se=TRUE, level=0.95) +
  ggplot2::xlab("Standardized MGAM2") + ggplot2::ylab("Standardized GPX2") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
summary(lm(MGAM2 ~ PPP1R1B, data=MD))
ggplot2::ggplot(MD, ggplot2::aes(x=PPP1R1B, y=MGAM2)) +
  ggplot2::geom_point(pch=16, fill="black", alpha=0.25, size=3) + 
  ggplot2::geom_smooth(method=lm, color="red", fill="red", alpha=0.25, se=TRUE, level=0.95) +
  ggplot2::xlab("Standardized PPP1R1B") + ggplot2::ylab("Standardized MGAM2") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
summary(lm(PPP1R1B ~ GPX2, data=MD))
ggplot2::ggplot(MD, ggplot2::aes(x=GPX2, y=PPP1R1B)) +
  ggplot2::geom_point(pch=16, fill="black", alpha=0.25, size=3) + 
  ggplot2::geom_smooth(method=lm, color="red", fill="red", alpha=0.25, se=TRUE, level=0.95) +
  ggplot2::xlab("Standardized GPX2") + ggplot2::ylab("Standardized PPP1R1B") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
remove(MD)

# Panel G: Investigate a 2-var ensemble
# Export covariation plot as .pdf at 7" x 5", edit, and export as .png
colnames(TABK)[4+c(2,8)]
PTRN2 <- SPACE::learn_pattern(CENK, c("O1.1_S1.2", "O1.1_S1.8"), 150, list(O1="#000000", S1=PALK), 
                               patch_list=PLSK, toroidal=T, plot_bkgd="W")

# Panel H: Map features of the 2-var ensemble back to the Visium spots
# Export map as .pdf at 6" x 7", convert to png, and crop.
# Export palette as metafile at 225 x 220
PTRN2S <- PTRN2[,c("X","Y","V")]
PTRN2S <- tidyr::pivot_wider(PTRN2S, names_from = V, values_from = Y)
TABK2S <- TABK[,c(1:4,4+c(2,8))]
MAP2 <- FNN::knnx.index(PTRN2S[,(ncol(PTRN2S)-1):ncol(PTRN2S)], TABK2S[,(ncol(TABK2S)-1):ncol(TABK2S)], k=1)
MAP2 <- MAP2 / nrow(TABK2S)
MAP2 <- unlist(lapply(MAP2, function(x) {sum(x >= c(0, 0.65, 0.77, 0.90))}))
plot(KCL$Y, max(KCL$X) - KCL$X, pch=16, col=c("gray", PALK[2], PALK[8],"gray")[MAP2])
SPACE::plot_palette(c(PALK[2], PALK[8]))

# Panel I: Correspondence of covariation-feature spots with k-means spots
# Export both as .pdf at 5.5" x 2.5", convert to png, and crop
BP <- data.frame(ME = 1:8, Count = unname(table(factor(KCL$Cluster[MAP2 == 2], levels=c(1:8)))))
BP$Count.Var1 <- NULL
colnames(BP)[2] <- "Count"
ggplot2::ggplot(data=BP, ggplot2::aes(x=as.factor(ME), y=Count, fill=as.factor(ME))) +
  ggplot2::geom_bar(stat="identity") + ggplot2::scale_fill_manual(values = KPL) +
  ggplot2::ylab("Feature 1 Spots") + ggplot2::xlab("Original ME") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
tab <- matrix(c(sum(MAP2==2 & KCL$Cluster == 3), sum(MAP2!=2 & KCL$Cluster == 3),
                sum(MAP2==2 & KCL$Cluster != 3), sum(MAP2!=2 & KCL$Cluster != 3)), nrow=2)
psych::cohen.kappa(tab)
remove(BP, tab)
BP <- data.frame(ME = 1:8, Count = unname(table(factor(KCL$Cluster[MAP2 == 3], levels=c(1:8)))))
BP$Count.Var1 <- NULL
colnames(BP)[2] <- "Count"
ggplot2::ggplot(data=BP, ggplot2::aes(x=as.factor(ME), y=Count, fill=as.factor(ME))) +
  ggplot2::geom_bar(stat="identity") + ggplot2::scale_fill_manual(values = KPL) +
  ggplot2::ylab("Feature 2 Spots") + ggplot2::xlab("Original ME") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
tab <- matrix(c(sum(MAP2==3 & KCL$Cluster == 3), sum(MAP2!=3 & KCL$Cluster == 3),
                sum(MAP2==3 & KCL$Cluster != 3), sum(MAP2!=3 & KCL$Cluster != 3)), nrow=2)
psych::cohen.kappa(tab)
remove(BP, tab)

# Panel S D: correlation among all pairs of variables
MD <- as.data.frame(scale(TABK2S[,(ncol(TABK2S)-1):ncol(TABK2S)]))
summary(lm(IGHG4 ~ SFRP2, data=MD))
ggplot2::ggplot(MD, ggplot2::aes(x=SFRP2, y=IGHG4)) +
  ggplot2::geom_point(pch=16, fill="black", alpha=0.25, size=3) + 
  ggplot2::geom_smooth(method=lm, color="red", fill="red", alpha=0.25, se=TRUE, level=0.95) +
  ggplot2::xlab("Standardized SFRP2") + ggplot2::ylab("Standardized IGHG4") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
remove(MD)
