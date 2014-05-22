M <- 2
MINEPS <- 1E-10

source('etriutils.R')
source('etributils.R')


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



createB1Constraint <- function(etrib, nCrit) {
  varnames <- c()
  for(j in 1:nCrit){
    name <- paste("w",j, sep="")
    varnames = c(varnames, name)
  }
  etrib$constr$lhs <- etriutils.addVariables(etrib$constr$lhs, varnames)
  
  weigthConst <- matrix(0, ncol=ncol(etrib$constr$lhs), nrow=nrow(etrib$constr$lhs), dimnames=list("B1"))
  etrib$constr$lhs <- rbind(etrib$constr$lhs, weigthConst)
  
  for(name in varnames){
    etrib$constr$lhs[nrow(etrib$constr$lhs), name] <- 1
  } 
  
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(1, ncol=1, nrow=1, dimnames=list("B1")))
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix("==", ncol=1, nrow=1, dimnames=list("B1")))
  
  return(etrib)
}

createB2Constraint <- function(etrib, nCrit, nProf) {
  varnames <- c()
  for(i in 1:nCrit){
    name <- paste0("c", i, "(b", nProf-1, ", b0)")
    varnames = c(varnames, name)
  }
  etrib$B2StartIndex <- ncol(etrib$constr$lhs)
  etrib$constr$lhs <- etriutils.addVariables(etrib$constr$lhs, varnames)
  etrib$B2EndIndex <- ncol(etrib$constr$lhs)
  
  for(i in 1:nCrit){
    cName <- paste0("c", i, "(b", nProf-1, ", b0)")
    wName <- paste0("w", i)
    lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list("B2", colnames(etrib$constr$lhs)))
    lhs[,cName] <- (-1)
    lhs[,wName] <- 1
    etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
    
    etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(0, ncol=1, nrow=1, dimnames=list("B2")))
    etrib$constr$dir <- rbind(etrib$constr$dir, matrix("==", ncol=1, nrow=1, dimnames=list("B2")))
  }
  
  return(etrib)
}
