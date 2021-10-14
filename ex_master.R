# Set working directory to "Code_and_data" folder
# Outputs are saved in the "results" folder

# clear memory !
rm(list = ls())


# load required libraries
library("EnvStats")        # For normal mixture distribution
library("SpatialExtremes") # For GEV distribution 
library("rmutil")          # For Lap distribution 
library("doParallel")      # For parallel computing
library("survPresmooth")   # To get presmoothed kernel density estimation
library("mvtnorm")         # To get multivariate critical point
library("multcomp")        # To get contrast matrix 
library("mratios")         # To get numerator and denominator contrast matrix
library("doRNG")           # To control random number generation
library("gridExtra")       # To create tables


# Load all required function
source("code/functions.R")      


R = 10000 # Number of Replications
B = 5000  # Number of resampling for bootstrap SCI
# You can change both values of 'R' and 'B' to 1000 for faster running of
# codes in approximately 2 hrs.


###########################################
###   Table 1: Differences of medians   ###
###########################################
t1 <- Sys.time()
N <- c(10, 30, 100, 500)

# Asymptotic SCI: True density 
NtdMVT1  <- matrix(NA, nrow = 4, ncol = 2)
EtdMVT1  <- matrix(NA, nrow = 4, ncol = 2)
CtdMVT1  <- matrix(NA, nrow = 4, ncol = 2)
LtdMVT1  <- matrix(NA, nrow = 4, ncol = 2)
GtdMVT1  <- matrix(NA, nrow = 4, ncol = 2)
NMtdMVT1 <- matrix(NA, nrow = 4, ncol = 2)
for(i in 1:4){NtdMVT1 [i,]<-SimD(N[i], R = R, distn = "Norm", den="td")}
for(i in 1:4){EtdMVT1 [i,]<-SimD(N[i], R = R, distn = "Exp" , den="td")}
for(i in 1:4){CtdMVT1 [i,]<-SimD(N[i], R = R, distn = "Cau" , den="td")}
for(i in 1:4){LtdMVT1 [i,]<-SimD(N[i], R = R, distn = "Lap" , den="td")}
for(i in 1:4){GtdMVT1 [i,]<-SimD(N[i], R = R, distn = "GEV" , den="td")}
for(i in 1:4){NMtdMVT1[i,]<-SimD(N[i], R = R, distn = "NM"  , den="td")}

# Asymptotic SCI: Kernel density 
NormMVT1 <- matrix(NA, nrow = 4, ncol = 2)
ExpMVT1  <- matrix(NA, nrow = 4, ncol = 2)
CauMVT1  <- matrix(NA, nrow = 4, ncol = 2)
LapMVT1  <- matrix(NA, nrow = 4, ncol = 2)
GEVMVT1  <- matrix(NA, nrow = 4, ncol = 2)
NMMVT1   <- matrix(NA, nrow = 4, ncol = 2)
for(i in 1:4){NormMVT1[i,]<-SimD(N[i], R = R, distn="Norm")}
for(i in 1:4){ExpMVT1 [i,]<-SimD(N[i], R = R, distn="Exp" )}
for(i in 1:4){CauMVT1 [i,]<-SimD(N[i], R = R, distn="Cau" )}
for(i in 1:4){LapMVT1 [i,]<-SimD(N[i], R = R, distn="Lap" )}
for(i in 1:4){GEVMVT1 [i,]<-SimD(N[i], R = R, distn="GEV" )}
for(i in 1:4){NMMVT1  [i,]<-SimD(N[i], R = R, distn="NM"  )}

# Bootstrap SCI
NormB1 <- matrix(NA, nrow = 4, ncol = 2)
ExpB1  <- matrix(NA, nrow = 4, ncol = 2)
CauB1  <- matrix(NA, nrow = 4, ncol = 2)
LapB1  <- matrix(NA, nrow = 4, ncol = 2)
GEVB1  <- matrix(NA, nrow = 4, ncol = 2)
NMB1   <- matrix(NA, nrow = 4, ncol = 2)
for(i in 1:4){NormB1[i,] <- BootSim(N[i], distn = "Norm", R = R, B = B)}
for(i in 1:4){ExpB1 [i,] <- BootSim(N[i], distn = "Exp" , R = R, B = B)}
for(i in 1:4){CauB1 [i,] <- BootSim(N[i], distn = "Cau" , R = R, B = B)}
for(i in 1:4){LapB1 [i,] <- BootSim(N[i], distn = "Lap" , R = R, B = B)}
for(i in 1:4){GEVB1 [i,] <- BootSim(N[i], distn = "GEV" , R = R, B = B)}
for(i in 1:4){NMB1  [i,] <- BootSim(N[i], distn = "NM"  , R = R, B = B)}

Sys.time() - t1

fMVT1    <- rbind(NtdMVT1,  EtdMVT1, CtdMVT1, LtdMVT1, GtdMVT1, NMtdMVT1)
fhatMVT1 <- rbind(NormMVT1, ExpMVT1, CauMVT1, LapMVT1, GEVMVT1, NMMVT1  )
BtSim1   <- rbind(NormB1,   ExpB1,   CauB1,   LapB1,   GEVB1,   NMB1    )
n1       <- c(N, N, N, N, N, N)
dat1     <- cbind(n1, fMVT1, fhatMVT1, BtSim1)
pdf("results/Table1.pdf", width = 6, height = 9, 
    title = "Coverage Probabilities (CP) for differences of medians")
grid.table(GetTable(dat1))
dev.off()


#####################################
###   Table 2: Ratios of medians  ###
#####################################
t1 <- Sys.time()
N <- c(10, 30, 100, 500)

# Asymptotic: True density with MVT critical point
NtdMVT2  <- matrix(NA, nrow = 4, ncol = 2)
EtdMVT2  <- matrix(NA, nrow = 4, ncol = 2)
CtdMVT2  <- matrix(NA, nrow = 4, ncol = 2)
LtdMVT2  <- matrix(NA, nrow = 4, ncol = 2)
GtdMVT2  <- matrix(NA, nrow = 4, ncol = 2)
NMtdMVT2 <- matrix(NA, nrow = 4, ncol = 2)
for(i in 1:4){NtdMVT2 [i,]<-SimR(N[i],distn="Norm", den="td", R = R)}
for(i in 1:4){EtdMVT2 [i,]<-SimR(N[i],distn="Exp" , den="td", R = R)}
for(i in 1:4){CtdMVT2 [i,]<-SimR(N[i],distn="Cau" , den="td", R = R)}
for(i in 1:4){LtdMVT2 [i,]<-SimR(N[i],distn="Lap" , den="td", R = R)}
for(i in 1:4){GtdMVT2 [i,]<-SimR(N[i],distn="GEV" , den="td", R = R)}
for(i in 1:4){NMtdMVT2[i,]<-SimR(N[i],distn="NM"  , den="td", R = R)}

# Asymptotic: Kernel density with MVT critical point
NormMVT2 <- matrix(NA, nrow = 4, ncol = 2)
ExpMVT2  <- matrix(NA, nrow = 4, ncol = 2)
CauMVT2  <- matrix(NA, nrow = 4, ncol = 2)
LapMVT2  <- matrix(NA, nrow = 4, ncol = 2)
GEVMVT2  <- matrix(NA, nrow = 4, ncol = 2)
NMMVT2   <- matrix(NA, nrow = 4, ncol = 2)
for(i in 1:4){NormMVT2[i,]<-SimR(N[i],distn="Norm", R = R)}
for(i in 1:4){ExpMVT2 [i,]<-SimR(N[i],distn="Exp" , R = R)}
for(i in 1:4){CauMVT2 [i,]<-SimR(N[i],distn="Cau" , R = R)}
for(i in 1:4){LapMVT2 [i,]<-SimR(N[i],distn="Lap" , R = R)}
for(i in 1:4){GEVMVT2 [i,]<-SimR(N[i],distn="GEV" , R = R)}
for(i in 1:4){NMMVT2  [i,]<-SimR(N[i],distn="NM"  , R = R)}

# Bootstrap
NormB2 <- matrix(NA, nrow = 4, ncol = 2)
ExpB2  <- matrix(NA, nrow = 4, ncol = 2)
CauB2  <- matrix(NA, nrow = 4, ncol = 2)
LapB2  <- matrix(NA, nrow = 4, ncol = 2)
GEVB2  <- matrix(NA, nrow = 4, ncol = 2)
NMB2   <- matrix(NA, nrow = 4, ncol = 2)
for(i in 1:4){NormB2[i,] <- BootSim(N[i], distn="Norm", ratio=T, R=R, B=B)}
for(i in 1:4){ExpB2 [i,] <- BootSim(N[i], distn="Exp" , ratio=T, R=R, B=B)}
for(i in 1:4){CauB2 [i,] <- BootSim(N[i], distn="Cau" , ratio=T, R=R, B=B)}
for(i in 1:4){LapB2 [i,] <- BootSim(N[i], distn="Lap" , ratio=T, R=R, B=B)}
for(i in 1:4){GEVB2 [i,] <- BootSim(N[i], distn="GEV" , ratio=T, R=R, B=B)}
for(i in 1:4){NMB2  [i,] <- BootSim(N[i], distn="NM"  , ratio=T, R=R, B=B)}

Sys.time() - t1

fMVT2    <- rbind(NtdMVT2,  EtdMVT2, CtdMVT2, LtdMVT2, GtdMVT2, NMtdMVT2)
fhatMVT2 <- rbind(NormMVT2, ExpMVT2, CauMVT2, LapMVT2, GEVMVT2, NMMVT2  )
BSim2    <- rbind(NormB2,   ExpB2,   CauB2,   LapB2,   GEVB2,   NMB2    )
n2       <- c(N, N, N, N, N, N)
dat2     <- cbind(n2, fMVT2, fhatMVT2, BSim2)
pdf("results/Table2.pdf", width = 6, height = 9, 
    title = "Coverage Probabilities (CP) for ratios of medians")
grid.table(GetTable(dat2))
dev.off()


###################
###   Figure 1  ###
###################
# open a pdf file
pdf("results/Figure1.pdf", width = 11, height = 5.5, 
    title = "Bias in density estimation at the median for Cauchy(0,1)")
# create plots
par(mfrow = c(2, 4))
FigTopPanel( 10, "Cau")
FigTopPanel( 30, "Cau")
FigTopPanel(100, "Cau")
FigTopPanel(500, "Cau")
FigLowerPanel( 10, "Cau")
FigLowerPanel( 30, "Cau")
FigLowerPanel(100, "Cau")
FigLowerPanel(500, "Cau")
par(mfrow = c(1, 1))
# close the pdf file
dev.off()


###################
###   Figure 2  ###
###################
# open a pdf file
pdf("results/Figure2.pdf", width = 11, height = 5.5, 
    title = "Bias in density estimation at the median for Laplace(0,1)")
# create plots
par(mfrow = c(2, 4))
FigTopPanel( 10, "Lap")
FigTopPanel( 30, "Lap")
FigTopPanel(100, "Lap")
FigTopPanel(500, "Lap")
FigLowerPanel( 10, "Lap")
FigLowerPanel( 30, "Lap")
FigLowerPanel(100, "Lap")
FigLowerPanel(500, "Lap")
par(mfrow = c(1,  1))
# close the pdf file
dev.off()


###################
###   Figure 3  ###
###################
# open a pdf file
pdf("results/Figure3.pdf", width = 11, height = 5.5, 
    title = "Bias in density estimation at the median for Normal Mixture")
# create plots
par(mfrow = c(2, 4))
FigTopPanel( 10, "NM")
FigTopPanel( 30, "NM")
FigTopPanel(100, "NM")
FigTopPanel(500, "NM")
FigLowerPanel( 10, "NM")
FigLowerPanel( 30, "NM")
FigLowerPanel(100, "NM")
FigLowerPanel(500, "NM")
par(mfrow = c(1, 1))
# close the pdf file
dev.off()


###################
###   Figure 4  ###
###################
load("data/phonemic.RData")
pdf("results/Figure4.pdf", width = 5, height = 4.8, 
    title = "Boxplots for improvement scores")
boxplot(y ~ f, data = phonemic, xlab = "", ylab = "Improvement scores")
dev.off()


###################
###   Table  3  ###
###################
load("data/phonemic.RData")
y     <- phonemic$y
f     <- phonemic$f
ni    <- tapply(f, f, length)
cmatt <- contrMat(n = ni, type = "Tukey", base = 3)[c(3, 1, 2),]
cmatd <- cmatt[-1,]
phonem   <- mcpqDCI(y, f, cmat = cmatt)
Estimate <- phonem$estimate
Column67 <- phonem$conf.int
Column23 <- mcpqDCI(y, f, cmat = cmatd)$conf.int
Column45 <- Bootci(y, f, p = 0.5, B = 5000, cmat = cmatd, seed = TRUE)
rownames(Column45) <- rownames(cmatd)
colnames(Column45) <- c("Lower", "Upper")
Column89 <- Bootci(y, f, p = 0.5, B = 5000, cmat = cmatt, seed = TRUE)
rownames(Column89) <- rownames(cmatt)
colnames(Column89) <- c("Lower", "Upper")
Column23 <- round(Column23, 2)
Column67 <- round(Column67, 2)
tab3  <- data.frame(Estimate,rbind(c("",""),cbind(Column23, Column45)),
                    cbind(Column67,Column89))
pdf("results/Table3.pdf", width = 9, height = 2, 
    title = "The 95% SCI for differences of medians")
grid.table(tab3)
dev.off()


###################
###   Table  4  ###
###################
load("data/bnct.RData")
tab4 <- summary(survfit(Surv(time, death) ~ trt, data = bnct))
group <- rep(c("Untreated", "Radiated", "RadiatedBPA"), c(8, 6, 7))
tab4 <- data.frame(group, time = tab4$time, R = tab4$n.risk, E=tab4$n.event, 
                   Shat = tab4$surv, std.err = round(tab4$std.err, 4))
pdf("results/Table4.pdf", width = 5, height = 7, 
    title = "Boxplots for improvement scores")
grid.table(tab4)
dev.off()


###################
###   Table  5  ###
###################
load("data/bnct.RData")
f     <- c(rep("Untreated",10), rep("Radiated",10), rep("RadiatedBPA",10))
y     <- bnct$time
event <- bnct$death 
cmatt <- matrix(c(-1, 1, 0, 1, 0, -1, 0, 1, -1), 3, byrow = TRUE)
rownames(cmatt) <- c("RadiatedBPA - Radiated", "Radiated - Untreated",
                     "RadiatedBPA - Untreated")
cmatd    <- cmatt[-1, ]
cmattB   <- matrix(c(0, -1, 1, -1, 1, 0, -1, 0, 1), 3, byrow = TRUE)
Col45 <- BootSURV(y, f, dat = bnct, event, cmat = cmattB[-1,])
Col89 <- BootSURV(y, f, dat = bnct, event, cmat = cmattB, seed = 4)
Estim <- round(mcpqDCI(y, f, event, Right.Censored = TRUE, 
                       p = 0.5, cmat = cmatt)$estimate, 2)
Col23 <- round(mcpqDCI(y, f, event, Right.Censored = TRUE, 
                       cmat = cmatd)$conf.int, 2)
Col67 <- round(mcpqDCI(y, f, event, Right.Censored = TRUE, 
                       cmat = cmatt)$conf.int, 2)
tab5  <- data.frame(Estim, rbind(c("", ""), cbind(Col23, Col45)),
                    cbind(Col67, Col89))
pdf("results/Table5.pdf", width = 9, height = 2, 
    title = "The 95% SCI for differences of median survival times")
grid.table(tab5)
dev.off()


###################
###   Table  6  ###
###################
load("data/bnct.RData")
f     <- c(rep("Untreated",10), rep("Radiated",10), rep("RadiatedBPA",10))
y     <- bnct$time
event <- bnct$death 
NumCt <- matrix(c(0, 1, 0, 1, 0, 0, 0, 1, 0), 3, byrow = TRUE)
DenCt <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 1), 3, byrow = TRUE)
rownames(NumCt) <- c("RadiatedBPA / Radiated", "Radiated / Untreated",
                     "RadiatedBPA / Untreated")
rownames(DenCt) <- rownames(NumCt)
NumCtB <- matrix(c(0, 0, 1, 0, 1, 0, 0, 0, 1), 3, byrow = TRUE)
DenCtB <- matrix(c(0, 1, 0, 1, 0, 0, 1, 0, 0), 3, byrow = TRUE)
Co45 <- round(BootSURV(y,f, dat = bnct, event, NumC=NumCtB[-1,],
                       DenC=DenCtB[-1,],ratio=T), 2)
Co89 <- round(BootSURV(y,f, dat=bnct, event, NumC=NumCtB,     
                       DenC=DenCtB,     ratio=T), 2)
Esti <- round(mcpqRCI(y, f, event, Right.Censored=T, Num.cmat=NumCt,      
                      Den.cmat = DenCt     )$estimate, 2)
Co23 <- round(mcpqRCI(y, f, event, Right.Censored=T, Num.cmat=NumCt[-1,], 
                      Den.cmat = DenCt[-1,])$conf.int, 2)
Co67 <- round(mcpqRCI(y, f, event, Right.Censored=T, Num.cmat=NumCt,      
                      Den.cmat = DenCt     )$conf.int, 2)
rownames(Esti) <- rownames(NumCt)
tab6  <- data.frame(Esti, rbind(c("",""), cbind(Co23, Co45)),
                    cbind(Co67,Co89))
pdf("results/Table6.pdf", width = 9, height = 2, 
    title = "The 95% SCI for ratios of median survival times")
grid.table(tab6)
dev.off()

