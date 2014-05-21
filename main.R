library(ROI.plugin.glpk)
library(ROI)
solvers <- ROI_installed_solvers()
if (!is.na(solvers['symphony'])) {
  .solver <<- 'symphony'
} else if (!is.na(solvers['glpk'])) {
  .solver <<- 'glpk'
} else {
  stop("No ROI Symphony or GLPK plugin installed")
}
source('etrib.R')

#warianty
alts <- read.table(file="alts.csv", sep=",", header=TRUE)
rownames(alts) = alts[,1]
alts <- alts[,2:ncol(alts)]

#granice klas
profs <- read.table(file="profs.csv", sep=",", header=FALSE)
rownames(profs) = profs[,1]
profs <- profs[,2:ncol(profs)]
colnames(profs) <- colnames(alts)

#przykładowe progi
thresholds <- matrix(c(
  0, 0.01, 0, 0.02, FALSE,
  0, 0, 1.9, 0, FALSE,
  0, 0, 1.9, 0, FALSE,
  0, 0, 1.9, 0, FALSE,
  0, 0, 2, 0, FALSE),ncol=5, byrow=TRUE)

## przykładowe przdziały do klas
assigs1 <- matrix(c(
  1, 1,
  5, 2),  , ncol=2, byrow=TRUE)
  
message("--- starting tests, iteration 1")

etrib.init(alts, profs, assigs1, th=thresholds)
