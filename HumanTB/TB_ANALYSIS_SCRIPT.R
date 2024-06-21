##### Make changes to SPACE code ##############################################
setwd("C:/Users/schromec/Desktop/SPACE/SPACE/R")
devtools::document()
devtools::install()
#restart R session
setwd("C:/Users/schromec/Desktop/SPACE/TB")
library(SPACE)

# # SINGLE USE CODE: trim TB cell phenotype maps to the correct size and reorganize
# load("TB_CPM.Rdata")
# for (i in 1:length(CPMs)) {
#   for (j in 1:length(CPMs[[i]])) {                                              # For each CPM
#     if (dim(CPMs[[i]][[j]])[1] > 964) {                                         # If larger than 964 x 964
#       border <- (dim(CPMs[[i]][[j]])[1] - 964) / 2                              # Find width of border
#       CPMs[[i]][[j]] <- CPMs[[i]][[j]][(border+1):(dim(CPMs[[i]][[j]])[1] - border), # Trim off border
#                                        (border+1):(dim(CPMs[[i]][[j]])[2] - border), ]
#     }                                                                           # Save as 4d array
#     CPMs[[i]][[j]] <- array(CPMs[[i]][[j]], dim=c(dim(CPMs[[i]][[j]])[1:2], 1, 1))
#   }
# }
# TB_CPM <- vector("list",30)
# for (i in 1:length(CPMs)) {                                                     # Make nested list into list
#   for (j in 1:length(CPMs[[i]])) {
#     TB_CPM[[((i-1)*2)+j]] <- CPMs[[i]][[j]]
#   }
# }
# remove(CPMs)
# save(TB_CPM, file="TB_CPM.Rdata")
load("TB_CPM.Rdata")
load("TB_CPM_NMS.Rdata")
load("TB_CPM_PAL.Rdata")

# # SINGLE USE CODE: trim TB max probability maps to the correct size and reorganize
# load("TB_MPM.Rdata")
# for (i in 1:length(MPM)) {
#   for (j in 1:length(MPM[[i]])) {                                               # For each CPM
#     if (dim(MPM[[i]][[j]])[1] > 964) {                                          # If larger than 964 x 964
#       border <- (dim(MPM[[i]][[j]])[1] - 964) / 2                               # Find width of border
#       MPM[[i]][[j]] <- MPM[[i]][[j]][(border+1):(dim(MPM[[i]][[j]])[1] - border), # Trim off border
#                                      (border+1):(dim(MPM[[i]][[j]])[2] - border), ]
#     }                                                                           # Save as 4d array
#     MPM[[i]][[j]] <- array(MPM[[i]][[j]], dim=c(dim(MPM[[i]][[j]])[1:2], 1, 1))
#   }
# }
# TB_MPM <- vector("list",30)
# for (i in 1:length(MPM)) {                                                    # Make nested list into list
#   for (j in 1:length(MPM[[i]])) {
#     TB_MPM[[((i-1)*2)+j]] <- MPM[[i]][[j]]
#   }
# }
# remove(MPM)
# save(TB_MPM, file="TB_MPM.Rdata")
load("TB_MPM.Rdata")
load("TB_MPM_NMS.Rdata")
load("TB_MPM_PAL.Rdata")

# # SINGLE USE CODE: trim raw marker images to the correct size and reorganize
# load("TB_TIF.Rdata")
# for (i in 1:length(TIF)) {
#   for (j in 1:length(TIF[[i]])) {                                               # For each CPM
#     if (dim(TIF[[i]][[j]])[1] > 964) {                                          # If larger than 964 x 964
#       border <- (dim(TIF[[i]][[j]])[1] - 964) / 2                               # Find width of border
#       TIF[[i]][[j]] <- TIF[[i]][[j]][(border+1):(dim(TIF[[i]][[j]])[1] - border), # Trim off border
#                                      (border+1):(dim(TIF[[i]][[j]])[2] - border), ]
#     }                                                                           # Save as 4d array
#     TIF[[i]][[j]] <- array(TIF[[i]][[j]], dim=c(dim(TIF[[i]][[j]])[1:2], 1, dim(TIF[[i]][[j]])[3]))
#   }
# }
# TB_TIF <- vector("list",30)
# for (i in 1:length(TIF)) {                                                      # Make nested list into list
#   for (j in 1:length(TIF[[i]])) {
#     TB_TIF[[((i-1)*2)+j]] <- TIF[[i]][[j]]
#   }
# }
# remove(TIF)
# save(TB_TIF, file="TB_TIF.Rdata")
load("TB_TIF.Rdata")
load("TB_TIF_NMS.Rdata")
load("TB_TIF_PAL.Rdata")

# # SINGLE USE CODE: Establish meta data file
# TB_GRP <- read.csv("Patient_metakey.csv", header=T)
# TB_GRP <- TB_GRP[rep(1:15, each=2), ]
# TB_GRP$FoV <- rep(1:2, times=15)
# save(TB_GRP, file="TB_GRP.Rdata")
load("TB_GRP.Rdata")

# # SINGLE USE CODE: Establish link table
# TB_LNK <- SPACE::load_table("Link_Table.csv", "L")
# save(TB_LNK, file="TB_LNK.Rdata")
load("TB_LNK.Rdata")

# Determine parameters for censusing. Use 10um radius, because 20 and 50 are too large & pseudoreplicate
TB_PAR_RAD <- SPACE::suggest_radii(10, c(X=0.5, Y=0.5, Z=1))
tb_par_num <- vector("list", length(TB_CPM))
for (i in 1:length(tb_par_num)) {
  tb_par_num[[i]] <- SPACE::suggest_number(5, TB_PAR_RAD, list(O1=TB_CPM[[i]], O2=TB_MPM[[i]], S1=TB_TIF[[i]]))
  print(i)
}
TB_PAR_NUM <- round(mean(unlist(tb_par_num)))
#save(TB_PAR_RAD, file="TB_PAR_RAD.Rdata")
#save(TB_PAR_NUM, file="TB_PAR_NUM.Rdata")

# Census the cell phenotype maps, max probability maps, and raw functional markers together at 20um
TB_CEN <- vector("list", length(TB_CPM))
TB_PLS <- vector("list", length(TB_CPM))
for (i in 1:length(TB_CPM)) {
  cen <- SPACE::census_image(list(O1 = TB_CPM[[i]], O2 = TB_MPM[[i]], S1 = TB_TIF[[i]]), 
                              OS_pairs = list(O1 = TB_LNK), radii = TB_PAR_RAD, sample_size = TB_PAR_NUM)
  TB_CEN[[i]] <- cen[[1]]
  TB_PLS[[i]] <- cen[[2]]
  remove(cen)
  print(paste("FINISHED #", i, sep=""))
  print(paste("Number of NAs: ", sum(is.na(TB_CEN[[i]])), sep=""))
}
TB_CEN <- SPACE::standardize_censuses(TB_CEN)
#save(TB_CEN, file="TB_CEN.Rdata")
#save(TB_PLS, file="TB_PLS.Rdata")

# Panel A and B: Investigate alpha diversity of cell types across FoVs
# Export bar graphs as .pdf at 5" x 3" and convert to png.
# Export box plots as .pdf at 4" x 4" and convert to png.
SPACE::alpha_diversity(TB_CPM[[30]], "O", list(O1 = TB_CPM_PAL))
SPACE::alpha_diversity(TB_CPM[[27]], "O", list(O1 = TB_CPM_PAL))
ADIV <- data.frame(Status = TB_GRP$Status, Div = rep(0,30))
for (i in 1:30) {
  ADIV$Div[i] <- SPACE::alpha_diversity(TB_CPM[[i]], "O", list(O1 = TB_CPM_PAL))[[2]]
}
ggplot2::ggplot(data=ADIV) +
  ggplot2::geom_boxplot(ggplot2::aes(x=Status, y=Div), outlier.shape=NA) + 
  ggplot2::geom_jitter(ggplot2::aes(x=Status, y=Div), size=5, alpha=0.5, width = 0.25) + 
  ggplot2::ylab("Alpha Diversity (bits)") + ggplot2::xlab("Clinical Status") + 
  ggplot2::theme_classic() + ggplot2::ylim(c(1.5,3.75)) + 
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), 
                 panel.grid.major=ggplot2::element_blank(),
                 axis.text=ggplot2::element_text(size=12),
                 axis.title=ggplot2::element_text(size=15),
                 plot.title=ggplot2::element_text(size=20))
mod <- lm(Div ~ Status, data=ADIV)
summary(mod)
emmeans::emmeans(mod, pairwise~Status, infer=T, adjust="none")                  

# Panel C and D: Investigate beta diversity of cell types in MEs across FoVs
# Export bar graphs as .pdf at 5" x 3" and convert to png.
# Export box plots as .pdf at 4" x 4" and convert to png.
SPACE::beta_diversity(list(O2 = TB_MPM[[3]], O1 = TB_CPM[[3]]), "O", list(O2 = TB_MPM_PAL, O2 = TB_CPM_PAL))
SPACE::beta_diversity(list(O2 = TB_MPM[[23]], O1 = TB_CPM[[23]]), "O", list(O2 = TB_MPM_PAL, O2 = TB_CPM_PAL))
BDIV <- data.frame(Status = TB_GRP$Status, Div = rep(0,30))
for (i in 1:30) {
  BDIV$Div[i] <- SPACE::beta_diversity(list(O2 = TB_MPM[[i]], O1 = TB_CPM[[i]]), "O", 
                                        list(O2 = TB_MPM_PAL, O2 = TB_CPM_PAL))[[2]]
}
ggplot2::ggplot(data=BDIV) +
  ggplot2::geom_boxplot(ggplot2::aes(x=Status, y=Div), outlier.shape=NA) + 
  ggplot2::geom_jitter(ggplot2::aes(x=Status, y=Div), size=5, alpha=0.5, width = 0.25) + 
  ggplot2::ylab("Beta Diversity (bits)") + ggplot2::xlab("Clinical Status") + 
  ggplot2::theme_classic() + 
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), 
                 panel.grid.major=ggplot2::element_blank(),
                 axis.text=ggplot2::element_text(size=12),
                 axis.title=ggplot2::element_text(size=15),
                 plot.title=ggplot2::element_text(size=20))
mod <- lm(Div ~ Status, data=BDIV)
summary(mod)
emmeans::emmeans(mod, pairwise~Status, infer=T, adjust="none")

# Find out which variables are present in most images: 3+/6 resection, 3+/6 postmortem, and 9+/18 biopsy
MV <- vector("list", length(TB_CEN))
for (i in 1:length(TB_CEN)) {
  MV[[i]] <- colSums(TB_CEN[[i]][TB_CEN[[i]]$Radius == 10, 1:(ncol(TB_CEN[[i]])-4)]) > 0
}
MV <- simplify2array(MV)
MV <- data.frame(Var = rownames(MV), Resection = rowSums(MV[,TB_GRP$Status=="Resection"]),
                Postmortem = rowSums(MV[,TB_GRP$Status=="Postmortem"]),
                Biopsy = rowSums(MV[,TB_GRP$Status=="Biopsy"]))
MV$Use <- (MV$Resection >= 3) & (MV$Postmortem >= 3) & (MV$Biopsy >= 9)

# Measure transMI for ensembles up to depth 3 
# Export as .pdf at 5" x 5" and convert to .png. Also export color palettes.
TB_TMI <- SPACE::measure_transMI(censuses = TB_CEN, groups = TB_GRP[,3,drop=F], depth = 3, radii = 10,
                                  not = MV$Var[!MV$Use])
#save(TB_TMI, file="TB_TMI.Rdata")
TB_TMIS <- SPACE::plot_MI_rank(TB_TMI, 10, 1:3, list(O1 = TB_CPM_PAL, S1 = TB_TIF_PAL), p_thr = 1e-7, 
                                group="Status")
TB_TMI[[1]][order(TB_TMI[[1]]$Pvalue_Status), ]
SPACE::plot_palette(TB_CPM_PAL, horizontal=T)
SPACE::plot_palette(TB_TIF_PAL, horizontal=T)

# Investigate the pattern formed by O1.10
# Export covariation plot as .pdf at 5" x 5" and convert to .png
# Export box plot as .pdf at 4" x 5" and convert to .png
PTRN1 <- SPACE::learn_pattern(TB_CEN, c("O1.10"), radius = 10, 
                      col_pal = list(O1 = TB_CPM_PAL, S1 = TB_TIF_PAL), group = TB_GRP[,3,drop=F], 
                      focal = c("Resection", "Postmortem"), reference = "Biopsy")
VRFY <- data.frame(Status = TB_GRP$Status, True_Frac = rep(0,30), Null_Frac = rep(0,30))
for (i in 1:nrow(VRFY)) {
  if (sum(TB_CEN[[i]]$O1.10) > 0) {
    avg <- sum(TB_CPM[[i]] == 10) / sum(TB_CPM[[i]] > 0)
    VRFY$True_Frac[i] <- sum(TB_CEN[[i]]$O1.10 > 90) / nrow(TB_CEN[[i]])
    rnd_cen <- SPACE::plot_dist(TB_CEN[[i]], "O1.10", 10, patch_list = TB_PLS[[i]])
    VRFY$Null_Frac[i] <- sum(rnd_cen > 90) / length(rnd_cen)
    remove(rnd_cen, avg)
    print(i)
  }
}
VRFY$LogRatio <- log(VRFY$True_Frac / VRFY$Null_Frac)
VRFY$LogRatio[VRFY$True_Frac == 0 | VRFY$Null_Frac == 0] <- 0
ggplot2::ggplot(data=VRFY) +
  ggplot2::geom_boxplot(ggplot2::aes(x=Status, y=LogRatio), outlier.shape=NA) + 
  ggplot2::geom_jitter(ggplot2::aes(x=Status, y=LogRatio), size=5, alpha=0.5, width = 0.25) + 
  ggplot2::ylab("CD68 Mac Aggregates (Log Ratio vs. Null)") + ggplot2::xlab("Clinical Status") + 
  ggplot2::theme_classic() + ggplot2::geom_hline(yintercept = 0, linetype="dashed") + ggplot2::ylim(c(-0.5, 2.5)) +
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), 
                 panel.grid.major=ggplot2::element_blank(),
                 axis.text=ggplot2::element_text(size=12),
                 axis.title=ggplot2::element_text(size=15),
                 plot.title=ggplot2::element_text(size=20))
mod <- lm(LogRatio ~ Status, data=VRFY)
summary(mod)
emmeans::emmeans(mod, pairwise~Status, infer=T, adjust="none")
save(VRFY, file = "TB_PTRN1_VRFY.Rdata")

# Investigate the top pattern involving O1.10 and linked scalar expression. O1.10 and O1.13_S1.9\
# Export covariation plot as .pdf at 6" x 5" and convert to .png
# Export scatter plot as .pdf at 6" x 5" and convert to .png
PTRN2 <- SPACE::learn_pattern(TB_CEN, c("O1.7_S1.3", "O1.10", "O1.13_S1.9"), radius = 10, 
                               col_pal = list(O1 = TB_CPM_PAL, S1 = TB_TIF_PAL), group = TB_GRP[,3,drop=F], 
                               focal = c("Resection", "Postmortem"), reference = "Biopsy",
                               smooth_window = 500)
AN <- dplyr::bind_rows(TB_CEN, .id = "Sample")                                  # Merge all nbhds across images
AN <- AN[,c("O1.10","O1.7_S1.3","O1.13_S1.9","Sample")]                         # Only relevant variables
AN$Sample <- as.numeric(AN$Sample)                                              # Make sample ID numeric
AN$Status <- TB_GRP$Status[AN$Sample]                                           # Get clinical status
AN$ID <- TB_GRP$Patient.ID[AN$Sample]                                           # Get patient ID
AN <- AN[AN$O1.13_S1.9 > 0, ]                                                   # Only nbhds with some O1.13_S1.9
table(AN$Status)                                                                # Very few resection nbhds
AN$Status[AN$Status %in% c("Resection","Postmortem")] <- "RescPost"             # So lump with postmortem
plot(density(AN$O1.13_S1.9))                                                    # O1.13_S1.9 is skewed
AN$Log_O1.13_S1.9 <- log(AN$O1.13_S1.9)                                         # So log transform it
mod <- lmerTest::lmer(Log_O1.13_S1.9 ~ O1.10*Status + (1|ID), data = AN)        # Model O1.13_S1.9 expression
summary(mod)
dummy <- expand.grid(O1.10 = seq(0,100,by=1), Status = c("RescPost","Biopsy"))
pred <- predict(mod, newdata=dummy, re.form=NA, se.fit=TRUE)
dummy$Mean <- pred[[1]]
dummy$Hi <- pred[[1]] + 1.96*pred[[2]]
dummy$Lo <- pred[[1]] - 1.96*pred[[2]]
AN$O1.10[AN$O1.10 == 0 & AN$Status == "Biopsy"] <- runif(sum(AN$O1.10 == 0 & AN$Status == "Biopsy"), min=-1, max=0)
AN$O1.10[AN$O1.10 == 0 & AN$Status == "RescPost"] <- runif(sum(AN$O1.10 == 0 & AN$Status == "RescPost"), min=-2, max=-1)
ggplot2::ggplot() +
  ggplot2::geom_point(ggplot2::aes(x=O1.10, y=Log_O1.13_S1.9, color=Status), data=AN, size=2, alpha = 0.25) + 
  ggplot2::geom_line(ggplot2::aes(x=O1.10, y=Mean, color=Status), data=dummy, linewidth = 2) +
  ggplot2::geom_ribbon(ggplot2::aes(x=O1.10, ymax=Hi, ymin=Lo, fill=Status), data=dummy, alpha=0.2, show.legend=F) +
  ggplot2::ylab("Log PDL1 on CD11b/c/CD206 Macs") + ggplot2::xlab("Local % of CD68 Macs") + 
  ggplot2::scale_color_manual(values = c("dodgerblue2", "coral2"), name="Clinical Status", 
                              labels = c("Biopsy","Resection or\nPostmortem")) +
  ggplot2::scale_fill_manual(values=c("coral2", "dodgerblue2")) +
  ggplot2::theme_classic() + 
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
                 axis.text=ggplot2::element_text(size=12), axis.title=ggplot2::element_text(size=15),
                 legend.title=ggplot2::element_text(size=15), legend.text=ggplot2::element_text(size=12))

# Investigate the top pattern involving O1.10 and linked scalar expression. O1.10 and O1.7_S1.3
# AN <- dplyr::bind_rows(TB_CEN, .id = "Sample")                                  # Merge all nbhds across images
# AN <- AN[,c("O1.10","O1.7_S1.3","O1.13_S1.9","Sample")]                         # Only relevant variables
# AN$Sample <- as.numeric(AN$Sample)                                              # Make sample ID numeric
# AN$Status <- TB_GRP$Status[AN$Sample]                                           # Get clinical status
# AN$ID <- TB_GRP$Patient.ID[AN$Sample]                                           # Get patient ID
# AN <- AN[AN$O1.7_S1.3 > 0, ]                                                    # Only nbhds with some O1.7_S1.3
# table(AN$Status)                                                                # Very few resection nbhds
# AN$Status[AN$Status %in% c("Resection","Postmortem")] <- "RescPost"             # So lump with postmortem
# plot(density(AN$O1.7_S1.3))                                                     # O1.7_S1.3 is skewed
# AN$Log_O1.7_S1.3 <- log(AN$O1.7_S1.3)                                           # So log transform it
# mod <- lmerTest::lmer(Log_O1.7_S1.3 ~ O1.10*Status + (1|ID), data = AN)         # Model O1.7_S1.3 expression
# dummy <- expand.grid(O1.10 = seq(0,100,by=1), Status = c("RescPost","Biopsy"))
# pred <- predict(mod, newdata=dummy, re.form=NA, se.fit=TRUE)
# dummy$Mean <- pred[[1]]
# dummy$Hi <- pred[[1]] + 1.96*pred[[2]]
# dummy$Lo <- pred[[1]] - 1.96*pred[[2]]
# AN$O1.10[AN$O1.10 == 0 & AN$Status == "Biopsy"] <- runif(sum(AN$O1.10 == 0 & AN$Status == "Biopsy"), min=-1, max=0)
# AN$O1.10[AN$O1.10 == 0 & AN$Status == "RescPost"] <- runif(sum(AN$O1.10 == 0 & AN$Status == "RescPost"), min=-2, max=-1)
# ggplot2::ggplot() +
#   ggplot2::geom_point(ggplot2::aes(x=O1.10, y=Log_O1.7_S1.3, color=Status), data=AN, size=2, alpha = 0.25) +
#   ggplot2::geom_line(ggplot2::aes(x=O1.10, y=Mean, color=Status), data=dummy, linewidth = 2) +
#   ggplot2::geom_ribbon(ggplot2::aes(x=O1.10, ymax=Hi, ymin=Lo, fill=Status), data=dummy, alpha=0.2, show.legend=F) +
#   ggplot2::ylab("Log IDO1 Level on Neutrophils") +
#   ggplot2::xlab("Local % of CD68+ Macrophages") +
#   ggplot2::scale_color_manual(values = c("dodgerblue2", "coral2"), name="Clinical Status",
#                               labels = c("Biopsy","Resection or\nPostmortem")) +
#   ggplot2::scale_fill_manual(values=c("coral2", "dodgerblue2")) +
#   ggplot2::theme_classic() +
#   ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(),
#                  axis.text=ggplot2::element_text(size=12), axis.title=ggplot2::element_text(size=15),
#                  legend.title=ggplot2::element_text(size=15), legend.text=ggplot2::element_text(size=12))

# Plot picture of top pattern involving O1.10
table(AN$Sample)                                                                # Choose samples 7, 9, 20, 29
# Sample 7 (Patient 4 FoV 1)
SPACE::plot_image(TB_CPM[[7]], "O", TB_CPM_PAL, objects = c(10,13))
IMG <- TB_TIF[[7]][,,,9,drop=F]
IMG[TB_CPM[[7]] != 13] <- 0
SPACE::plot_image(IMG, "S", "white", enh_cnt = 0)
# Sample 9 (Patient 5 FoV 1)
SPACE::plot_image(TB_CPM[[9]], "O", TB_CPM_PAL, objects = c(10,13))
IMG <- TB_TIF[[9]][,,,9,drop=F]
IMG[TB_CPM[[9]] != 13] <- 0
SPACE::plot_image(IMG, "S", "white", enh_cnt = 0)
# Sample 28 (Patient 14 FoV 2)
SPACE::plot_image(TB_CPM[[28]], "O", TB_CPM_PAL, objects = c(10,13))
IMG <- TB_TIF[[28]][,,,9,drop=F]
IMG[TB_CPM[[28]] != 13] <- 0
SPACE::plot_image(IMG, "S", "white", enh_cnt = 0)
# Sample 29 (Patient 15 FoV 1)
SPACE::plot_image(TB_CPM[[29]], "O", TB_CPM_PAL, objects = c(10,13))
IMG <- TB_TIF[[29]][,,,9,drop=F]
IMG[TB_CPM[[29]] != 13] <- 0
SPACE::plot_image(IMG, "S", "white", enh_cnt = 0)
