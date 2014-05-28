M <- 1E+10
MINEPS <- 1E-10

source('etriutils.R')
source('etributils.R')
source('etribbase.R')
source('etribcardinalities.R')
source('etribpairwisecomp.R')


etrib.init <- function(performances, profiles, assignments,monotonicity, th, cardinalities, pCk, pCl) {
  stopifnot(ncol(performances) == ncol(profiles))
  stopifnot(ncol(assignments) == 3)
  stopifnot(ncol(pCk) == 3)
  stopifnot(ncol(pCl) == 3)
  stopifnot(nrow(assignments) < nrow(performances))
  stopifnot(nrow(th) == ncol(performances))
  stopifnot(nrow(th) == length(monotonicity))
  profiles <- createBorderProfiles(profiles, assignments)
  
  message("--- constructing base model")
  
  A <- 1:nrow(performances)
  H <- 1:nrow(profiles)
  J <- 1:ncol(performances)
  
  etrib <- list()
  varnames <- createTRIBVariables(A, H, J, pCk, pCl)
  etrib$constr$lhs <- intiMatrix(varnames)
  etrib$constr$dir <- initDIRMatrix()
  etrib$constr$rhs <- intiRHSMatrix()
  
  etrib <- buildBaseModel(etrib, performances, profiles, monotonicity, th, A, H, J)
  etrib <- buildAEModel(etrib, J, assignments)
  etrib <- buildCCModel(etrib, A, H, J, cardinalities)
  etrib <- buildPConstraint(etrib, J, H, pCk, pCl)
  
  return(etrib)
}


buildBaseModel <- function(etrib, performances, profiles, monotonicity, th, A, H, J) {
  etrib <- createB1Constraint(etrib, J)
  etrib <- createB2Constraint(etrib, J, H)
  etrib <- createB4Constraint(etrib)
  etrib <- createB5Constraint(etrib, J)
  etrib <- createB6Constraint(etrib, performances, profiles,monotonicity, th)
  
  return(etrib)
}

buildAEModel <- function(etrib, J, assignments){
  nAs <- nrow(assignments)
  nrows <- nAs * 2
  rnames <- paste0("AE.",1:nrows)
  lhs <- matrix(0, nrow=nrows, ncol=ncol(etrib$constr$lhs), dimnames=list(rnames,colnames(etrib$constr$lhs)))
  
  row <- 0
  for(i in 1:nAs){
    a <- assignments[i,]
    row <- row + 1
    lhs[row,paste0("c",J,"(a",a[1],",b",a[2],")")] <- 1
    lhs[row,"L"] <- -1
  }
  
  for(i in 1:nAs){
    a <- assignments[i,]
    row <- row + 1
    lhs[row,paste0("c",J,"(a",a[1],",b",a[3],")")] <- 1
    lhs[row,"L"] <- -1
    lhs[row,"e"] <- 1
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  
  dir <- as.matrix(rep(c(">=","<="), nAs))
  rownames(dir) <- rnames
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  rhs <- as.matrix(rep(0, row))
  rownames(rhs) <- rnames
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  
  return(etrib)
}

buildCCModel <- function(etrib, A, H, J, cardinalities){
  etrib <- createCC1Constraints(etrib, A, H)
  etrib <- createCC2Constraints(etrib, A, H, J)
  etrib <- createCC3Constraints(etrib, A, H, J)
  etrib <- createCC4Constraints(etrib, A, cardinalities)
  etrib <- createCC5Constraints(etrib, A, cardinalities)
  return(etrib)
}

buildPConstraint <- function(etrib, J, H, pCk, pCl){
  etrib <- createPC1Constrints(etrib, J, H, pCk)
  etrib <- createPC2Constrints(etrib, J, H, pCk)
  etrib <- createPC3Constrints(etrib, J, H, pCk)
  etrib <- createPU1Constrints(etrib, J, H, pCl)
  etrib <- createPU2Constrints(etrib, J, H, pCl)
  etrib <- createPU3Constrints(etrib, J, H, pCl)
  return(etrib)
}

createTRIBVariables <- function(A, H, J, pCk, pCl){
  varnames <- c()
  
  p <- length(H)
  
  varnames <- paste0("w",J)
  varnames <- c(varnames, c("L"))
  varnames <- c(varnames, c("e"))
  varnames = c(varnames, paste0("c", J, "(b", p, ",b0)"))

  for(a in A){
    for(h in H){
      varnames <- c(varnames, paste0("v(a",a,",h",h,")"))
    }
  }
  
  for (a in A) {
    for (b in c(0,H)) {
      varnames <- c(varnames, paste0('c', J, '(a', a, ',b', b, ')'))
    }
  }
  
  for (b in H) {    
    for (a in A) {
      varnames <- c(varnames, paste0('c', J, '(b', b, ',a', a, ')'))
    }
  }
  
  for(row in nrow(pCk)){
    pck <- pCk[row,]
    varnames <- c(varnames, paste0("v(a",pck[1],",b",pck[2],",>=",pck[3],",h",1:(p-pck[3]),")"))
  }
  
  for(row in nrow(pCl)){
    pcl <- pCl[row,]
    varnames <- c(varnames, paste0("v(a",pcl[1],",b",pcl[2],",<=",pcl[3],",h",1:(p-pcl[3]),")"))
  }
  
  return(varnames)
}

intiMatrix <- function(names){
  lhs <- matrix(0, ncol=length(names), nrow=1,dimnames=list("E",names))
  lhs["E","e"] <- 1
  return(lhs)
}

initDIRMatrix <- function(){
  dir <- matrix(c(">="))
  rownames(dir) <- c("E")
  return(dir)
}

intiRHSMatrix <- function(){
  rhs <- matrix(MINEPS)
  rownames(rhs) <- "E"
  return(rhs)
}

