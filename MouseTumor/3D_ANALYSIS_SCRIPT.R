### Load SPACE package
library(SPACE)   # Note that working directories referenced below will not be the same on your machine!

# Figure 4: Analysis of 3d Tumor Volume with CD4 T, CD8 T, and CD31 Objects

# Panel A: Imaris Screenshot of 3d Rendering

# Panel B: SPACE Analysis and Covariation Plot
# Export to PDF at 8.5" x 6.5", edit, convert to png.
CD4 <- read.csv("IBEX9_Tumor_P1_cropped_CD4TPositions_Detailed.csv")
CD8 <- read.csv("IBEX9_Tumor_P1_cropped_CD8TPositions_Detailed.csv")
CD31 <- read.csv("IBEX9_Tumor_P1_cropped_VesselPositions_Detailed.csv")
CD4$Object <- rep(1, nrow(CD4))
CD8$Object <- rep(2, nrow(CD8))
CD31$Object <- rep(3, nrow(CD31))
POS <- rbind(CD4, CD8, CD31)
remove(CD4, CD8, CD31)
RAD <- c(20,30,40,50)
NUM <- round(5 * (808 * 753 * 270) / ((4/3) * pi * RAD^3))
NUM[NUM > nrow(POS)] <- nrow(POS)
CENT <- SPACE::census_table(POS, radii = RAD, sample_size = NUM)
PLST <- CENT[[2]]
CENT <- CENT[[1]]
CMIT <- SPACE::measure_cisMI(CENT, PLST, 3, RAD)
SPACE::learn_pattern(CENT, c("O1.1","O1.2","O1.3"), 50, list(O1 = c("darkcyan","orange","green3")), 
                      patch_list = PLST, toroidal=F, smooth_window = 100)
save(CENT, file = "IBEX9_Tumor_P1_cropped_Census.Rdata")
save(PLST, file = "IBEX9_Tumor_P1_cropped_PatchList.Rdata")
save(CMIT, file = "IBEX9_Tumor_P1_cropped_CisMI.Rdata")

# Make logit function and inverse logit functions
logit <- function(p) {
  return(log(p/(1-p)))
}
invlogit <- function(s) {
  return(exp(s)/(1+exp(s)))
}

# Panel C: Compare Numbers of T Cells Near vs. Far from Vessels
# Export as pdf at 7" x 5", edit, and save as png
TCS <- data.frame(Object = POS[POS$Object < 3, 4], # Type of each T cell & whether it is close to a vessel
                  Close = as.vector(FNN::knnx.dist(POS[POS$Object == 3, 1:3], POS[POS$Object < 3, 1:3], k=1) < 50)) 
PD <- as.data.frame(expand.grid(Object = c(1,2), Close = c(T,F)))
PD$Num <- rep(0, nrow(PD))
for (i in 1:nrow(PD)) {
  PD$Num[i] <- sum(TCS$Object == PD$Object[i] & TCS$Close == PD$Close[i]) / sum(TCS$Close == PD$Close[i])
}
ggplot2::ggplot(PD, ggplot2::aes(x=forcats::fct_inorder(as.factor(Close)), y=Num, fill=as.factor(Object))) +
  ggplot2::geom_bar(position=ggplot2::position_dodge(), stat="identity", colour='black') +
  ggplot2::scale_fill_manual(values=c("darkcyan","orange"), labels=c("CD4+ T", "CD8+ T"), name="Cell Type") + 
  ggplot2::scale_x_discrete(breaks=c(T,F), labels=c("< 50", "> 50"), name="Distance to Nearest Vessel (\u03bcm)") +
  ggplot2::ylab("Fraction of Cells") + ggplot2::theme_bw() + ggplot2::ylim(c(0,1)) +
  ggplot2::theme(text = ggplot2::element_text(size=20))
mod <- glm((Object-1) ~ 1 + Close, data=TCS, family=binomial)
summary(mod)
emmeans::emmeans(mod, c("Close"), infer=T)
2 * (1 - pnorm(18.215, 0, 1))                     # P value for whether CD4 vs. CD8 counts differ far from vessel
2 * (1 - pnorm(6.774))                            # P value for whether CD4 vs. CD8 counts differ close to vessel

# Panel D: Compare Distances of T Cells to Vessels, Near vs. Far from Vessels
# Export as pdf at 7" x 5", edit, and save as png
TCS$Dist <- as.vector(FNN::knnx.dist(POS[POS$Object == 3, 1:3], POS[POS$Object < 3, 1:3], 1))
PD$Dist <- rep(0, nrow(PD))
PD$Dist_SE <- rep(0, nrow(PD))
for (i in 1:nrow(PD)) {
  PD$Dist[i] <- mean(TCS$Dist[TCS$Object == PD$Object[i] & TCS$Close == PD$Close[i]])
  PD$Dist_SE[i] <- sd(TCS$Dist[TCS$Object == PD$Object[i] & TCS$Close == PD$Close[i]]) / 
    sqrt(sum(TCS$Object == PD$Object[i] & TCS$Close == PD$Close[i]))
}
ggplot2::ggplot(PD, ggplot2::aes(x=forcats::fct_inorder(as.factor(Close)), y=Dist, fill=as.factor(Object))) +
  ggplot2::geom_bar(position=ggplot2::position_dodge(), stat="identity", colour='black') +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=Dist-Dist_SE, ymax=Dist+Dist_SE), width=.2,
                         position=ggplot2::position_dodge(0.9), linewidth = 1) +
  ggplot2::scale_fill_manual(values=c("darkcyan","orange"), labels=c("CD4+ T", "CD8+ T"), name="Cell Type") + 
  ggplot2::scale_x_discrete(breaks=c(T,F), labels=c("< 50", "> 50"), name="Distance to Nearest Vessel (\u03bcm)") +
  ggplot2::ylab("Distance to Nearest Vessel (\u03bcm)") + ggplot2::theme_bw() + ggplot2::ylim(c(0,100)) +
  ggplot2::theme(text = ggplot2::element_text(size=20))
mod <- glm(Dist ~ 1 + Object, data=TCS[TCS$Close, ], family=gaussian)
summary(mod)
mod <- glm(Dist ~ 1 + Object, data=TCS[!TCS$Close, ], family=gaussian)
summary(mod)
        
