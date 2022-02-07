library(readr)

GLS <- read_delim("/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/logDEG", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)

#GSC <- read_delim("/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/iTN656_mets.txt", 
 #                 "\t", escape_double = FALSE, trim_ws = TRUE)

#GSC <- read_delim("/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/iTN656_subSystems.txt", 
#                  "\t", escape_double = FALSE, trim_ws = TRUE)

#GSC <- read_delim("/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/iTN656_rxnNames.txt", 
#                  "\t", escape_double = FALSE, trim_ws = TRUE)

GSC <- read_delim("/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/iTN656_equations.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)

gsc_nachon <- GSC
names(gsc_nachon) <- c("s", "g")
gsc_nachon <- gsc_nachon[, c(2,1)]

gls_nachon <- GLS
gsc_nachon$g <- as.character(gsc_nachon$g)
gsc_nachon$s <- as.character(gsc_nachon$s)
library(piano)
gsc_nachon_obj <- loadGSC(gsc_nachon)
pvals <- gls_nachon$p
names(pvals) <- gls_nachon$g
directions <- gls_nachon$FC
names(directions) <- gls_nachon$g
nachon_data <- list(gsc = gsc_nachon_obj, pvals = pvals, directions = directions)
gsares <- runGSA(geneLevelStats= pvals , directions=directions,
                 gsc=gsc_nachon_obj, nPerm=500,
                 geneSetStat = "reporter")
dev.off()

hm <- GSAheatmap(gsares)

dev.off()
#View(GSAsummaryTable(gsares))
#GSAsummaryTable(gsares, save=T, 
 #               file="/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/gsares_mets.txt")
#GSAsummaryTable(gsares, save=T, 
#                file="/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/gsares_subSystems.txt")


#GSAsummaryTable(gsares, save=T, 
#                file="/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/gsares_rxnNames.txt")

GSAsummaryTable(gsares, save=T, 
                file="/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM/GSC/gsares_equations.txt")
