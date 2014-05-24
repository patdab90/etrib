M <- 2
MINEPS <- 1E-10

source('etriutils.R')
source('etributils.R')
source('etribbase.R')


etrib.init <- function(performances, profiles, assignments, th) {
  stopifnot(ncol(performances) == ncol(profiles))
  stopifnot(ncol(assignments) == 2)
  stopifnot(nrow(th) == ncol(performances))
  profiles <- createBorderProfiles(profiles, assignments)
  
  message("--- constructing base model")
  
  etrib <- list()
  etrib <- buildBaseModel(etrib, performances, profiles, assignments, th)
  
  return(etrib)
}


buildBaseModel <- function(etrib, performances, profiles, assignments, th) {
  nAlts <- nrow(performances)
  nCrit <- ncol(performances)
  nAssignments <- nrow(assignments)
  nCats <- nrow(profiles)
  
  
  etrib <- createEpsilonConstraint(etrib)
  etrib <- createB1Constraint(etrib, nCrit)
  etrib <- createB2Constraint(etrib, nCrit, nCats)
  etrib <- createB4Constraint(etrib)
  etrib <- createB5Constraint(etrib, nCrit)
  etrib <- createB6Constraint(etrib, performances, profiles, th)
  
  ##allConst <- combineConstraintsMatrix(ec)
  ##colnames(allConst$lhs) <- getColNames(nAlts, nCrit, nAssignments, nCats)
  return(etrib)
}

createEpsilonConstraint <- function(eb){
  lhs <- matrix(c(1))
  rownames(lhs) <- c("EC")
  dir <- matrix(c(">="))
  rownames(dir) <- c("EC")
  rhs <- matrix(MINEPS)
  rownames(rhs) <- c("EC")
  eb$constr <- list(lhs=lhs, dir=dir, rhs=rhs)
  eb$epsylonIndex <- 1
  return(eb)
}
