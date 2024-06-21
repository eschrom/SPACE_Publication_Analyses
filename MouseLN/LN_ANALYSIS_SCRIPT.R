### Load SPACE package
library(SPACE)   # Note that working directories referenced below will not be the same on your machine!

### FIGURE 1S: Semi-Supervised Categorization of LN Cell Types #################

# PANEL A: Overlay of several markers from the 42-plex mouse LN image
# The image is "LN_OrigImg.tif" and the panel is "LN_OrigPanel.csv".
# The panel is made via color overlays in ImageJ

# PANEL B: CellPose segmentation of the mouse LN image
# Select JOJO & CD45 channels in ImageJ to yield "LN_OrigImg_NucCD45.tif".
# Segment with default CellPose to get "LN_OrigImg_SegMask.tif".
# Extract the centroid of each cell and the MFI for each marker,
# using the "LN_ObjExpr_Script.ijm.ijm" script in ImageJ.

# Assemble data frame with centroid and marker expression per cell
setwd("C:/Users/schromec/Desktop/SPACE/MouseLN")
names <- read.csv("LN_OrigPanel.csv", header=T)$Marker
setwd("C:/Users/schromec/Desktop/SPACE/MouseLN/LN_ObjExpr")
files <- list.files(".")
for (i in 1:length(files)) {
  if (i == 1) {
    cells <- read.csv(files[i], header=T)
    cells <- cells[,2:3]
    colnames(cells) <- c("X","Y")
  } else {
    new_col <- read.csv(files[i], header=T)
    new_col <- new_col[,2,drop=F]
    colnames(new_col) <- names[i-1]
    cells <- cbind(cells, new_col)
    remove(new_col)
  }
}
setwd("C:/Users/schromec/Desktop/SPACE/MouseLN")
remove(names, files, i)

# Load expected expression profiles from "Fig2S_Lineage_Expectations.csv",
# also known as "TabS1.csv", and use only 20 lineage markers 
lin_exp <- read.csv("LN_LinExp.csv", header=T)
cols_to_use <- which(colnames(cells) %in% colnames(lin_exp))
lin_cls <- cells[,c(cols_to_use)]
remove(cols_to_use)

# Define expression thresholds using IsoData for each lineage marker
lin_mrk <- as.data.frame(matrix(NA, nrow=ncol(lin_cls), ncol=2))
colnames(lin_mrk) <- c("Marker", "Thresh")
lin_mrk$Marker <- colnames(lin_cls)
for (i in 1:nrow(lin_mrk)) {
  d <- lin_cls[,i]
  d <- sort(d)
  pos <- min(which(d > 0))
  thr <- d[pos]
  alo <- mean(d[1:pos])
  ahi <- mean(d[(pos+1):length(d)])
  while ((alo + ahi)/2 > thr) {
    thr <- (alo + ahi)/2
    pos <- max(which(d <= thr))
    alo <- mean(d[1:pos])
    ahi <- mean(d[(pos+1):length(d)])
  }
  remove(ahi, alo, pos)
  lin_mrk$Thresh[i] <- thr
  remove(thr,d)
}

# Establish a sigmoid function
sigmoid <- function(x, K, n) {
  # x = input value
  # K = half-maximal constant
  # n = exponent
  return((x^n)/(K^n + x^n))
}

# PANEL C: Binarization of a cellular MFI distribution
# Export at 500 x 500
i <- 2
ggplot2::ggplot() +
  ggplot2::geom_density(data=lin_cls, ggplot2::aes(x=lin_cls[,i], y=ggplot2::after_stat(scaled)),
                        linewidth=3) +
  ggplot2::geom_vline(xintercept = lin_mrk$Thresh[i], linetype="dashed", linewidth=1.5) +
  ggplot2::theme_bw() + ggplot2::xlab(paste("Cellular ", lin_mrk$Marker[i]," MFI",sep="")) +
  ggplot2::scale_y_continuous(name="Normalized Density") +
  ggplot2::theme(text = ggplot2::element_text(size=15)) +
  ggplot2::ggtitle("MFI Distribution, Binarized")

# PANEL D: Translating MFI into probability of expression
# Export at 500 x 500
example_sigmoid <- data.frame(x=0:255, y=sigmoid(0:255, lin_mrk$Thres[i], 4))
ggplot2::ggplot() +
  ggplot2::geom_density(data=lin_cls, ggplot2::aes(x=lin_cls[,i], y=ggplot2::after_stat(scaled)),
                        linewidth=3, color="#808080") +
  ggplot2::geom_line(data=example_sigmoid, ggplot2::aes(x=x, y=y), linewidth=3) +
  ggplot2::geom_vline(xintercept = lin_mrk$Thresh[i], linetype="dashed", linewidth=1.5) +
  ggplot2::theme_bw() + ggplot2::xlab(paste("Cellular ", lin_mrk$Marker[i]," MFI",sep="")) +
  ggplot2::scale_y_continuous(name="Normalized Density",
                              sec.axis=ggplot2::sec_axis(trans=~.*1, name="Expression Probability")) +
  ggplot2::theme(text = ggplot2::element_text(size=15),
                 axis.title.y.left = ggplot2::element_text(color="#808080"),
                 axis.line.y.left = ggplot2::element_line(color="#808080"),
                 axis.text.y.left = ggplot2::element_text(color="#808080"),
                 axis.ticks.y.left = ggplot2::element_line(color="#808080")) +
  ggplot2::ggtitle("Translating MFI to Probability")

# PANEL E: Distribution of expression probability
# Export at 500 x 500
ggplot2::ggplot() +
  ggplot2::geom_density(data=lin_cls, ggplot2::aes(x=sigmoid(lin_cls[,i],lin_mrk$Thres[i],4), 
                                                   y=ggplot2::after_stat(scaled)), linewidth=3) +
  ggplot2::theme_bw() + ggplot2::xlab(paste("Cellular ", lin_mrk$Marker[i]," Expression Probability",sep="")) +
  ggplot2::scale_y_continuous(name="Normalized Density") +
  ggplot2::theme(text = ggplot2::element_text(size=15)) +
  ggplot2::ggtitle("Expression Probability Distribution")
remove(example_sigmoid, i)

# Calculate probability scores of each cell type for each cell
lin_nms <- lin_exp$Name
lin_exp <- lin_exp[,2:ncol(lin_exp)]
lin_typ <- as.data.frame(matrix(0, nrow=nrow(lin_cls), ncol=nrow(lin_exp)))
colnames(lin_typ) <- lin_nms
lin_prb <- data.frame(Map(sigmoid, lin_cls, K = lin_mrk$Thresh, n = 4))
for (i in 1:ncol(lin_typ)) {
  pos <- which(lin_exp[i,] == 1)
  neg <- which(lin_exp[i,] == 0)
  lin_typ[,i] <- (Reduce('*', lin_prb[,pos,drop=F]) * Reduce('*', 1-lin_prb[,neg,drop=F])) ^ 
    (1/(length(c(pos,neg))))
  remove(pos, neg)
}

# PANEL F: Bar graphs depicting calculation of cell x type probability scores
# Export each at 750 x 250
example_cell <- tibble::rownames_to_column(as.data.frame(t(lin_prb[500,])),"Marker")
colnames(example_cell)[2] <- "Prob"
ggplot2::ggplot(example_cell) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Marker), y=Prob), fill="black") + 
  ggplot2::scale_y_continuous(name="Probability", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Marker") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) + 
  ggplot2::ggtitle("Cell 500 Expression Profile")
example_treg <- tibble::rownames_to_column(as.data.frame(t(lin_exp[3,])),"Marker")
colnames(example_treg)[2] <- "Prob"
ggplot2::ggplot(example_treg) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Marker), y=Prob), fill="gray") + 
  ggplot2::scale_y_continuous(name="Probability", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Marker") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
  ggplot2::ggtitle("Canonical T-reg Expression Profile")
example_match <- example_cell
example_match$Prob[example_treg$Prob == 0] <- 1 - example_match$Prob[example_treg$Prob == 0]
ggplot2::ggplot(example_match) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Marker), y=Prob), fill="#808080") + 
  ggplot2::scale_y_continuous(name="Probability", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Marker") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) + 
  ggplot2::ggtitle("Cell 500 Match to Canonical T-reg")
example_nk <- tibble::rownames_to_column(as.data.frame(t(lin_exp[5,])),"Marker")
colnames(example_nk)[2] <- "Prob"
ggplot2::ggplot(example_nk) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Marker), y=Prob), fill="gray") + 
  ggplot2::scale_y_continuous(name="Probability", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Marker") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
  ggplot2::ggtitle("Canonical NK Expression Profile")
example_match <- example_cell
example_nk$Prob[is.na(example_nk$Prob)] <- 1
example_match$Prob[example_nk$Prob == 0] <- 1 - example_match$Prob[example_nk$Prob == 0]
ggplot2::ggplot(example_match) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Marker), y=Prob), fill="#808080") + 
  ggplot2::scale_y_continuous(name="Probability", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Marker") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) + 
  ggplot2::ggtitle("Cell 500 Match to Canonical NK")
remove(example_cell, example_treg, example_nk, example_match)

# Determine a type for each cell, including "unknown"
lin_fnl <- data.frame(Type = apply(lin_typ, 1, which.max),
                      Prob = apply(lin_typ, 1, max))
lin_fnl$Type[lin_fnl$Prob <= 0.5] <- nrow(lin_exp) + 1

# PANEL G: Bar graphs depicting cell type assignment for each cell
# Export each at 750 x 250
example_cell <- tibble::rownames_to_column(as.data.frame(t(lin_typ[499,])),"Type")
colnames(example_cell)[2] <- "Prob"
ggplot2::ggplot(example_cell) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Type), y=Prob), fill="black") + 
  ggplot2::scale_y_continuous(name="Prob Score", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Cell Type") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) + 
  ggplot2::ggtitle("Cell 499 Cell Type Probability Scores")
example_cell <- tibble::rownames_to_column(as.data.frame(t(lin_typ[500,])),"Type")
colnames(example_cell)[2] <- "Prob"
ggplot2::ggplot(example_cell) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Type), y=Prob), fill="black") + 
  ggplot2::scale_y_continuous(name="Prob Score", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Cell Type") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) + 
  ggplot2::ggtitle("Cell 500 Cell Type Probability Scores")
example_cell <- tibble::rownames_to_column(as.data.frame(t(lin_typ[501,])),"Type")
colnames(example_cell)[2] <- "Prob"
ggplot2::ggplot(example_cell) +
  ggplot2::geom_col(ggplot2::aes(x=forcats::fct_inorder(Type), y=Prob), fill="black") + 
  ggplot2::scale_y_continuous(name="Prob Score", limits=c(0,1)) +
  ggplot2::theme_bw() + ggplot2::xlab("Cell Type") +
  ggplot2::theme(text=ggplot2::element_text(size=15), 
                 axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) + 
  ggplot2::ggtitle("Cell 501 Cell Type Probability Scores")
remove(example_cell)

# PANEL H: distribution of maximum lineage probability scores across cells
# Export at 750 x 350
ggplot2::ggplot(data = lin_fnl) +
  ggplot2::geom_density(ggplot2::aes(x=Prob, y=ggplot2::after_stat(scaled)),
                        linewidth=3, color="black") +
  ggplot2::geom_vline(xintercept = 0.5, linetype="dashed", linewidth=1.5) +
  ggplot2::theme_bw() + ggplot2::xlab("Maximum Cell Type Probability Score") +
  ggplot2::scale_y_continuous(name="Normalized Density") +
  #ggplot2::annotate("text", x=c(0.25,0.75), y=0.75,
  #                  label=paste(round(100*as.vector(table(lin_fnl$Prob > 0.5))/nrow(lin_fnl),1),
  #                              "%", sep=""), size=10) +
  ggplot2::theme(text = ggplot2::element_text(size=15))

### FIGURE 2: SPACE Exploration of Mouse LN ####################################

# Make a color palette for the 18 cell types
# Export at 250 x 750
PAL <- SPACE::make_palette(19)
SPACE::plot_palette(PAL, col_labels=paste(c(1:19), rep(". ",19), c(lin_nms, "Unknown"), sep=""), horizontal=T)
#save(PAL, file="LN_Palette.Rdata")

# PANEL A: heatmap of avg MFIs for each cell type
# Export at 750 x 500
tmp_typ <- lin_fnl$Type
TBL <- aggregate(lin_cls, by=list(tmp_typ), FUN=mean)
TBL <- cbind(TBL[,1], as.vector(table(tmp_typ)), 
             TBL[,2:ncol(TBL)])
colnames(TBL)[1:2] <- c("Object", "Count")
SPACE::plot_table(TBL, tile_plots=F)
remove(tmp_typ)
#save(TBL, file="LN_Profile_Table.Rdata")

# PANEL B: Color segmentation mask by cell type and save as a .tif
mask <- cytomapper::loadImages("LN_OrigImage_SegMask.tif")
mask <- mask$LN_OrigImage_SegMask@.Data
mask <- mask * 65535
rgb_pal <- col2rgb(PAL)
img_r <- array(0, dim=dim(mask))
img_g <- array(0, dim=dim(mask))
img_b <- array(0, dim=dim(mask))
for (i in 1:18) {
  img_r[mask %in% which(lin_fnl$Type == i)] <- rgb_pal[1,i]
  img_g[mask %in% which(lin_fnl$Type == i)] <- rgb_pal[2,i]
  img_b[mask %in% which(lin_fnl$Type == i)] <- rgb_pal[3,i]
  print(i)
}
img <- array(c(img_r, img_g, img_b), dim = c(dim(mask),3)) / 255
#tiff::writeTIFF(img, "LN_ObjImg.tif")
remove(img_r, img_g, img_b, rgb_pal, i, img)

# Color segmentation mask by cell type as an array, to avoid color scramble
IMG <- array(0, dim=dim(mask))
for (i in 1:19) {
  IMG[mask %in% which(lin_fnl$Type == i)] <- i
}
IMG <- array(IMG, dim=c(dim(IMG), 1, 1))
SPACE::plot_image(IMG, "O", PAL)
#save("IMG", file="LN_ObjImg.Rdata")
remove(mask, i)

# Choose parameters for censusing the object image
PAR_radii <- SPACE::suggest_radii(c(10,20,30), c(x=0.288,y=0.288,z=1))
PAR_number <- SPACE::suggest_number(5, PAR_radii, list(O1 = IMG))

# Census the image
CEN <- SPACE::census_image(list(O1 = IMG), radii = PAR_radii, sample_size = PAR_number)
PLS <- CEN[[2]]
CEN <- CEN[[1]]
#save(CEN, file="LN_Census.Rdata")
#save(PLS, file="LN_PatchList.Rdata")

# Measure cisMI for ensembles up to 3 at the 10 micron radius
CMI <- SPACE::measure_cisMI(CEN, PLS, 3, 10)
#save(CMI, file="LN_CisMI.Rdata")

# Panel C: Summary of top ensembles
# Export at 5.33" x 5.00" as pdf, edit, and convert to 300dpi png
SPACE::plot_MI_rank(CMI, 10, 1:3, list(O1 = PAL), mi_thr = 0.05) 

# Panel D: Covariation of LZ and DZ B cells. Inset includes CD4
# Export at 5.33" x 5.00" as pdf, edit, and convert to 300dpi png
SPACE::learn_pattern(CEN, c("O1.7","O1.8"), 10, list(O1 = PAL), patch_list = PLS, toroidal = T)
SPACE::learn_pattern(CEN, c("O1.1","O1.7","O1.8"), 10, list(O1 = PAL), toroidal = T)

# Panel E: Visualization of germinal center orientation

# Panel F: Covariation of CD4, CD8, and B.
# Export at 6.33" x 5.00" as pdf, edit, and convert to 300dpi png
SPACE::learn_pattern(CEN, c("O1.1","O1.2","O1.6"), 10, list(O1 = PAL), patch_list = PLS, toroidal = T,
                      smooth_window = 200)

# Panel G: Visualization of opposing T cell gradients

# Panel H: Covariation of CD4, CD8, cDC1, and cDC2
# Export at 5.87" x 5.00" as pdf, edit, and convert to 300dpi png
CMI_hyp <- SPACE::measure_cisMI(CEN, PLS, 4, 10, all=c("O1.1","O1.2","O1.13","O1.14"))
#save(CMI_hyp, file="LN_CisMI_Hypothesis.Rdata")
SPACE::learn_pattern(CEN, c("O1.1","O1.2","O1.13","O1.14"), 10, list(O1 = PAL), patch_list = PLS, toroidal = T,
                      smooth_window = 200)

### FIGURE S2: Comparison to Table-Based Analysis ##############################

# Create and save the member table
MEM <- cbind(cells[,1:2] * 0.288, rep(1,nrow(cells)), lin_fnl[,1])
colnames(MEM) <- c("X","Y","Z","Object")
#write.csv(MEM, file="TabS2.csv")

# Collect census on member table to match census on image
CENT <- SPACE::census_table(MEM, radii = c(10,20,30), sample_size = PAR_number)
PLST <- CENT[[2]]
CENT <- CENT[[1]]
#save(CENT, file="LN_Census_Table.Rdata")
#save(PLST, file="LN_PatchList_Table.Rdata")

# Calculate cisMI for the table-based census
CMIT <- SPACE::measure_cisMI(CENT, PLST, 3, 10)
#save(CMIT, file="LN_CisMI_Table.Rdata")
SPACE::plot_MI_rank(CMIT, 10, 1:3, list(O1 = PAL), mi_thr = 0.05)

# Only Panel: Correlation of Z-scores between image- and table-based approaches
# Export at 5" x 5" as pdf then convert to .png with 300 dpi
MIS <- CMI[[1]]
MIST <- CMIT[[1]]
MIS[is.na(MIS)] <- "none"
MIST[is.na(MIST)] <- "none"
MISS <- plyr::match_df(MIS, MIST, on = c("VA","VB","VC"))
MISTS <- plyr::match_df(MIST, MISS, on = c("VA","VB","VC"))
MISM <- dplyr::inner_join(MISS, MISTS, by = c("VA","VB","VC"))
remove(MISS, MISTS)
MISM <- MISM[,c(1:3,5,9)]
colnames(MISM)[4:5] <- c("Z_Image", "Z_Table")
MISM_mod <- lm(MISM$Z_Table ~ MISM$Z_Image)
ggplot2::ggplot(data=MISM, ggplot2::aes(x=Z_Image, y=Z_Table)) +
  ggplot2::geom_point(pch=16, fill="black", alpha=0.25, size=3) + 
  ggplot2::ylim(c(-75,150)) + 
  ggplot2::xlim(c(-30,120)) +
  ggplot2::geom_smooth(method=lm , color="red", fill="red", alpha=0.25, se=TRUE, level=0.95) +
  ggplot2::ylab("Z Score from Table") + ggplot2::xlab("Z Score from Image") + ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 text=ggplot2::element_text(size=20), legend.position="none")
summary(MISM_mod)
remove(MIS, MIST, MISM, MISM_mod)

### FIGURE 3S: Oblique Cut of LN Obscures Patterns Visually ####################

# All Panels are constructed from screen shots in Imaris

### FIGURE 4S: Existing Methods Are Not as Powerful as SPACE ####################

# Collect a new census using all cells instead of a defined number of neighborhoods
CENTF <- SPACE::census_table(MEM, radii = 10, sample_size = nrow(MEM))
PLSTF <- CENTF[[2]]
CENTF <- CENTF[[1]]
#save(CENTF, file="LN_Census_Table_Benchmark.Rdata")
#save(PLSTF, file="LN_PatchList_Table_Benchmark.Rdata")

# Panel A: Cartoon in Powerpoint

# Panel B: Cartoon in Powerpoint

# Panel C: Pairwise co-occurrence heatmap
# Export at 5.87" x 5.00" as pdf, edit, and convert to 300dpi png
NCZ <- cbind(expand.grid(c(1:19), c(1:19)), matrix(NA, nrow=19^2, ncol=101))
colnames(NCZ) <- c("VA", "VB", "T_Mn", paste("R_Mn",1:100,sep=""))
for (i in 1:101) {
  if (i == 1) {                                # For the first measurement, use true cell phenotypes
    mem <- MEM
  } else {                                     # For subsequent measurements, use permuted cell phenotypes
    mem$Object <- mem$Object[sample(nrow(mem))]
  }
  for (j in 1:nrow(NCZ)) {                     # For all cell type pairs, measure distance
    src <- mem[,1:2][mem$Object == NCZ[j,1], ] # from source cell to 60 nearest target cells
    trg <- mem[,1:2][mem$Object == NCZ[j,2], ] # 60 insures that the full 10 micron neighborhood is included
    dst <- FNN::knnx.dist(trg, src, k=min(60,nrow(trg))) 
    if (NCZ[j,1] == NCZ[j,2]) {                # If source and target are the same, remove smallest distance
      dst <- dst[,2:ncol(dst)]                 # which will be 0, as source and target are the same cell
    } 
    NCZ[j,i+2] <- mean(rowSums(dst <= 10))     # Count number of neighbors within 10um and average
    remove(dst,src,trg)
  } 
  print(paste("Finished iteration ", i, "/101", sep=""))
}        
remove(mem)
NCZ$Zscr <- (NCZ$T_Mn - rowMeans(NCZ[,4:103])) / apply(NCZ[,4:103], 1, sd)
NCZ$Zplt <- NCZ$Zscr
NCZ$Zplt[NCZ$Zplt < -30] <- -30
NCZ$Zplt[NCZ$Zplt > 30] <- 30
ggplot2::ggplot(data=NCZ) +
  ggplot2::geom_tile(ggplot2::aes(x=VA, y=VB, fill=Zplt)) + 
  ggplot2::scale_x_continuous(expand=c(0,0), name="Focal Type", breaks=c(1:19), labels=c(lin_nms, "Unknown")) + 
  ggplot2::scale_y_continuous(expand=c(0,0), name="Neighbor Type", breaks=c(1:19), labels=c(lin_nms, "Unknown")) +
  ggplot2::scale_fill_gradientn(colors=c("#FDE725","#FFFFFF", "#414487"), 
                                values=c(1, -min(NCZ$Zplt)/(max(NCZ$Zplt) - min(NCZ$Zplt)), 0),
                                name="Z Score", breaks=c(-30,-15,0,15,30), 
                                labels=c("<-30", "-15", "0", "15", ">30")) +
  ggplot2::theme_classic() + #ggplot2::ggtitle("Number of Neighbors vs. Randomized Expectations") +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), 
                 legend.title=ggplot2::element_text(size=20),
                 legend.text=ggplot2::element_text(size=18),
                 axis.title=ggplot2::element_text(size=20),
                 axis.text.y=ggplot2::element_text(size=18),
                 axis.text.x=ggplot2::element_text(size=18,angle=45,hjust=1,vjust=1),
                 plot.title=ggplot2::element_text(size=20))
save(NCZ, file="LN_Pairwise_Benchmarking.Rdata")

# Panel D: Cartoon in Powerpoint

# Panel E: Cartoon in Powerpoint

# Panel F: An example k-means clustering (k = 18) ME compositions
# Export at 6.00" x 4.00" as pdf, edit, and convert to 300dpi png
# Export palette as metafile at 720 x 200
INP <- CENTF[,1:19]
kmn_pal <- SPACE::make_palette(18)
kmn <- kmeans(INP, centers=18)
kmn_tab <- as.data.frame(cbind(c(1:18), unname(table(kmn$cluster)), kmn$centers))
colnames(kmn_tab) <- c("Object", "Count", c(lin_nms,"Unknown"))
#kmn_tab <- kmn_tab[c(3,8,16,5,11,18,15,12,14,9,7,10,1,17,6,4,2,13),]
kmn_tab$Object <- c(1:18)
SPACE::plot_table(kmn_tab)
SPACE::plot_palette(kmn_pal, axis_label = "ME")
#save(kmn_pal, file = "LN_MEPalette_Benchmarking.Rdata")
#save(kmn, file = "LN_MEClusters_Benchmarking.Rdata")

# Panel G: An example k-means clustering (k = 18) ME mapping
# Export at 6.5" x 5" as a pdf, convert to png at 300dpi
plot(CENTF$Y, max(CENTF$X) - CENTF$X, col=kmn_pal[c(13,17,1,16,4,15,11,2,10,12,5,8,18,9,7,3,14,6)][kmn$cluster], 
     pch=16, ylab="", yaxt='n', xlab="", xaxt='n') 

# Many K-means runs across different K values to see which key features are captured
DIST <- dist(INP)
#save(DIST, file="LN_DistMatrix_Benchmarking.Rdata")
rpk <- 50 # reps per k (50)
mxk <- 25 # max k to try (25)
n <- rpk * (mxk-1)
KMN <- data.frame(K = rep(2:mxk, each=rpk), S = rep(NA, n), W = rep(NA, n), B = rep(NA, n),
                  F1 = rep(0, n), F2 = rep(0, n), F3 = rep(0, n), F4 = rep(0, n), F5 = rep(0, n), 
                  F6 = rep(0, n), F7 = rep(0, n), F8 = rep(0, n), F9 = rep(0, n), F10 = rep(0, n))
for (i in 1:nrow(KMN)) {
  kmn <- kmeans(INP, centers=KMN$K[i], iter.max=30)
  KMN$S[i] <- mean(cluster::silhouette(kmn$cluster, DIST)[,3])
  KMN$W[i] <- kmn$tot.withinss
  KMN$B[i] <- kmn$betweenss
  z4 <- which(scale(kmn$centers[,1]) >= 0.5)                                    # Where CD4s are enriched
  z8 <- which(scale(kmn$centers[,2]) >= 0.5)                                    # Where CD8s are enriched
  zT <- 1:KMN$K[i] %in% union(z4,z8)                                        # Where CD4s & CD8s are enriched
  if (sum(zT) > 0) {                                                            # F1: T zone where cD4s & CD8s
    dummy <- rep(0, sum(zT))                                                    # are the top two types
    for (j in 1:sum(zT)) {                                                   
      if (sum(order(kmn$centers[which(zT)[j],], decreasing=T)[1:2] %in% c(1,2)) == 2) {
        dummy[j] <- 1
      }
    }
    KMN$F1[i] <- max(dummy)
  }
  if (sum(zT) > 1) {
    if (lm(kmn$centers[zT,1] ~ kmn$centers[zT,2])$coefficients[2] < 
        min(0, lm(kmn$centers[!zT,1] ~ kmn$centers[!zT,2])$coefficients[2], na.rm=T)) {  # F2: Opposing CD4 vs. 8 grad
      KMN$F2[i] <- 1
    }
    if (lm(kmn$centers[zT,1] ~ kmn$centers[zT,6])$coefficients[2] > 
        max(0, lm(kmn$centers[!zT,1] ~ kmn$centers[!zT,6])$coefficients[2], na.rm=T)) {  # F3: Naive B with CD4
      KMN$F3[i] <- 1
    }
    if (lm(kmn$centers[zT,6] ~ kmn$centers[zT,2])$coefficients[2] < 
        min(0, lm(kmn$centers[!zT,6] ~ kmn$centers[!zT,2])$coefficients[2], na.rm=T)) {  # F4: Naive B away from CD8
      KMN$F4[i] <- 1
    }
    if (lm(kmn$centers[zT,13] ~ kmn$centers[zT,14])$coefficients[2] < 
        min(0, lm(kmn$centers[!zT,13] ~ kmn$centers[!zT,14])$coefficients[2], na.rm=T)) {# F5: Opposing DC1 vs. 2 grad
      KMN$F5[i] <- 1
    }
    if (lm(kmn$centers[zT,1] ~ kmn$centers[zT,14])$coefficients[2] > 
        max(0, lm(kmn$centers[!zT,1] ~ kmn$centers[!zT,14])$coefficients[2], na.rm=T)) { # F6: DC2 with CD4
      KMN$F6[i] <- 1
    }
    if (lm(kmn$centers[zT,2] ~ kmn$centers[zT,13])$coefficients[2] > 
        max(0, lm(kmn$centers[!zT,2] ~ kmn$centers[!zT,13])$coefficients[2], na.rm=T)) { # F7: DC1 with CD8
      KMN$F7[i] <- 1
    }
  }
  zL <- which(scale(kmn$centers[,7]) >= 0.5)                                    # Where LZBs are enriched
  zD <- which(scale(kmn$centers[,8]) >= 0.5)                                    # Where DZBs are enriched
  zG <- 1:KMN$K[i] %in% union(zL,zD)                                        # Where L+DZBs are enriched
  if (sum(zG) > 0) {                                                            # F8: GC where LZB & DZB
    dummy <- rep(0, sum(zG))                                                    # are the top 2 types
    for (j in 1:sum(zG)) {                                                   
      if (sum(order(kmn$centers[which(zG)[j],], decreasing=T)[1:2] %in% c(7,8)) == 2) {
        dummy[j] <- 1
      }
    }
    KMN$F8[i] <- max(dummy)
  }
  if (sum(zG) > 1) {
    if (lm(kmn$centers[zG,7] ~ kmn$centers[zG,1])$coefficients[2] <             # F9: CD4 T away from LZB
        min(0, lm(kmn$centers[!zG,7] ~ kmn$centers[!zG,1])$coefficients[2], na.rm=T)) {  
      KMN$F9[i] <- 1
    }
    if (lm(kmn$centers[zG,8] ~ kmn$centers[zG,1])$coefficients[2] >             # F10: CD4 T with DZB
        max(0, lm(kmn$centers[!zG,8] ~ kmn$centers[!zG,1])$coefficients[2], na.rm=T)) {  
      KMN$F10[i] <- 1
    }
  }
  remove(z4, z8, zT, zL, zD, zG, dummy, kmn)
  print(i)
}  
KMN$FT <- rowSums(KMN[,c(5:14)])/10
remove(i,j)
#save(KMN, file="LN_MEFeatures_Benchmarking.Rdata")   

# Panel H: Features recovered by each of 50 clustering runs at k=15
# Export at 650 x 500 as metafile
kmn_plot <- as.data.frame(cbind(1:rpk, rep(1,rpk), (KMN[KMN$K==18, c(5:14)]))) 
colnames(kmn_plot) <- c("Object", "Count", paste("F",1:10,sep=""))
SPACE::plot_table(kmn_plot)

# Panel I: Compare silhouette score with feature recovery across K values
# Export 8.25" x 6.00" as pdf, edit, and convert to 300dpi png
KMN_plot <- data.frame(K = unique(KMN$K), S = aggregate(KMN$S, list(KMN$K), mean)[,2], 
                       S_SE = (aggregate(KMN$S, list(KMN$K), sd)[,2] / sqrt(rpk)),
                       FT = aggregate(KMN$FT, list(KMN$K), mean)[,2],
                       FT_SE = (aggregate(KMN$FT, list(KMN$K), sd)[,2] / sqrt(rpk)))
KMN_plot$S_SE <- KMN_plot$S_SE / KMN_plot$S
KMN_plot$S <- (KMN_plot$S - min(KMN_plot$S)) / (max(KMN_plot$S) - min(KMN_plot$S))
ggplot2::ggplot(data=KMN_plot, ggplot2::aes(x=K)) +
  ggplot2::geom_line(ggplot2::aes(y=FT), col="red", linewidth=3) +
  ggplot2::geom_line(ggplot2::aes(y=S), linewidth=3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=FT-FT_SE, ymax=FT+FT_SE), col="red", width=0.25, linewidth=1.5) + 
  ggplot2::geom_errorbar(ggplot2::aes(ymin=S-S_SE, ymax=S+S_SE), width=0.25, linewidth=1.5) + 
  ggplot2::coord_cartesian(ylim=c(0,1)) + 
  ggplot2::scale_y_continuous(expand=c(0,0), name="Normalized Average Silhouette Score",
                              sec.axis=ggplot2::sec_axis(~ ., name="Fraction of Key Features Detected")) + 
  ggplot2::scale_x_continuous(expand=c(0,0), breaks=c(2:mxk), labels=c(2:mxk), name="Number of MEs") + 
  ggplot2::theme_classic() + ggplot2::ggtitle("Global Clustering Quality vs. Specific Pattern Accuracy") +
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                 axis.title.x=ggplot2::element_text(size=18),
                 axis.title.y=ggplot2::element_text(size=18),
                 axis.title.y.right=ggplot2::element_text(size=18, color="red"), 
                 axis.text.y=ggplot2::element_text(size=15),
                 axis.text.y.right=ggplot2::element_text(size=15, color="red"),
                 axis.text.x=ggplot2::element_text(size=15),
                 plot.title=ggplot2::element_text(size=20))

### FIGURE S5: Details of Microenvironment Simulations #########################

# Panel A: Features recovered by each clustering run across all K
# Export at 18" x 6" as pdf, edit, convert to .png at 300 dpi
kmn_plot <- as.data.frame(cbind(1:n, rep(1,n), (KMN[, c(5:14)])))
colnames(kmn_plot) <- c("Object", "Count", paste("F",1:10,sep=""))
SPACE::plot_table(kmn_plot)

# Panel B: Number of features recovered by each clustering run at each K
# Export at 11" x 18" as portrait pdf, then convert to .png at 300dpi
KMN$Run <- rep(1:rpk, (mxk-1))
KMN$K_text <- paste("K = ", KMN$K, sep="") # 1150 x 1800 as .tif
ggplot2::ggplot(data=KMN, ggplot2::aes(x=Run, y=FT*10)) +
  ggplot2::geom_line(linewidth=3) + ggplot2::coord_cartesian(ylim=c(0,10)) + 
  ggplot2::theme_classic() +
  ggplot2::geom_hline(yintercept=0, linetype="dashed") +
  ggplot2::scale_y_continuous(name="Number of Features Recovered", breaks=seq(0,10,2), labels=seq(0,10,2)) + 
  ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                 axis.title.x=ggplot2::element_text(size=20),
                 axis.title.y=ggplot2::element_text(size=20),
                 axis.text.y=ggplot2::element_text(size=15),
                 axis.text.x=ggplot2::element_text(size=15),
                 strip.text=ggplot2::element_text(size=15)) +
  ggplot2::facet_wrap(~factor(K_text, levels = paste("K = ", 2:mxk, sep="")), ncol=2)
max(KMN$FT)
remove(mxk, n, rpk)
