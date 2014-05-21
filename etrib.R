M <- 2
MINEPS <- 1E-10

source('etributils.R')

etrib.init <- function(performances, profiles, assignments, th) {
  stopifnot(ncol(performances) == ncol(profiles))
  stopifnot(ncol(assignments) == 2)
  stopifnot(nrow(th) == ncol(performances))
  
  message("--- constructing base model")
  
  etrib <- NA
  etrib <- buildBaseModel(etrib, performances, profiles, assignments, th)
  
  return(etrib)
}


buildBaseModel <- function(etrib, performances, profiles, assignments, th) {
  nAlts <- nrow(performances)
  nCrit <- ncol(performances)
  nAssignments <- nrow(assignments)
  nCats <- nrow(profiles)
  
  eb <- createEpsilonConstraint(etrib)

  eb
  eb.epsilonIndex
  
  ##allConst <- combineConstraintsMatrix(ec)
  ##colnames(allConst$lhs) <- getColNames(nAlts, nCrit, nAssignments, nCats)
  return(eb)
}

createEpsilonConstraint <- function(eb){
  lhs <- matrix(c(1))
  rownames(lhs) <- c("EC")
  dir <- data.frame(c(">="))
  rownames(dir) <- c("EC")
  rhs <- c()
  rownames(rhs) <- c("EC")
  eb.constr <- list(lhs=lhs, dir=dir, rhs=rhs)
  eb.epsilonIndex <- 1
  
}
