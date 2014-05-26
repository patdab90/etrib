M <- 2
MINEPS <- 1E-10

source('etriutils.R')
source('etributils.R')
source('etribbase.R')


etrib.init <- function(performances, profiles, assignments,monotonicity, th) {
  stopifnot(ncol(performances) == ncol(profiles))
  stopifnot(ncol(assignments) == 3)
  stopifnot(nrow(th) == ncol(performances))
  stopifnot(nrow(th) == length(monotonicity))
  profiles <- createBorderProfiles(profiles, assignments)
  
  message("--- constructing base model")
  
  etrib <- list()
  etrib <- buildBaseModel(etrib, performances, profiles, assignments,monotonicity, th)
  ##etrib <- createAEConstraint(etrib,  performances, assignments)
  
  return(etrib)
}


buildBaseModel <- function(etrib, performances, profiles, assignments,monotonicity, th) {
  nAlts <- nrow(performances)
  nCrit <- ncol(performances)
  nAssignments <- nrow(assignments)
  nCats <- nrow(profiles)
  
  
  etrib <- createEpsilonConstraint(etrib)
  etrib <- createB1Constraint(etrib, nCrit)
  etrib <- createB2Constraint(etrib, nCrit, nCats)
  etrib <- createB4Constraint(etrib)
  etrib <- createB5Constraint(etrib, nCrit)
  etrib <- createB6Constraint(etrib, performances, profiles,monotonicity, th)
  
  ##allConst <- combineConstraintsMatrix(ec)
  ##colnames(allConst$lhs) <- getColNames(nAlts, nCrit, nAssignments, nCats)
  return(etrib)
}

createAEConstraint <- function(etrib, performances, assignments){
  nCrit <- ncol(performances)
  nrows <- nrow(assignments)*2
  rnames <- paste0("AE.",1:nrows)
  lhs <- matrix(0, nrow=nrows, ncol=ncol(etrib$constr$lhs), dimnames=list(rnames,colnames(etrib$constr$lhs)))
  
  row <- 0
  for(i in 1:nrow(assignments)){
    a <- assignments[i,]
    row <- row + 1
    for(j in 1:nCrit){
      nameL <- paste0("c",j,"(a",a[1],",b",a[2],")")
      lhs[row,nameL] <- 1
    }
    lhs[row,"L"] <- -1
  }
  
  for(i in 1:nrow(assignments)){
    a <- assignments[i,]
    row <- row + 1
    for(j in 1:nCrit){
      nameU <- paste0("c",j,"(a",a[1],",b",a[3],")")
      lhs[row,nameU] <- 1
    }
    lhs[row,"L"] <- -1
    lhs[row,"e"] <- 1
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  
  dir <- as.matrix(rep(c(">=","<="), each=row/2))
  rownames(dir) <- rnames
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  rhs <- as.matrix(rep(0, row))
  rownames(rhs) <- rnames
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  
  return(etrib)
}

createCC1Constraints <- function(etrib, A, H){
  varnames <- c()
  for(a in A){
    for(h in H){
      varnames <- c(varnames, paste0("v(a",a,",h",h,")"))
    }
  }
  
  etrib$constr$lhs <- etriutils.addVariables(etrib$constr$lhs, varnames)
  

  lhs <- matrix(0, nrow=length(A), ncol=ncol(etrib$constr$lhs), dimnames=list(paste0("CC1.",1:length(A)),colnames(etrib$constr$lhs)))
  
  rows <- 0
  for(a in A){
    rows <- rows + 1
    for(h in H){
      v <-  paste0("v(a",a,",h",h,")")
      lhs[rows, v] <- 1
    }
  }
  

  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  
  return(etrib)
}

createCC2Constraints <- function(etrib, A, H, J){
  lhs <- matrix(0, nrow=length(A), ncol=ncol(etrib$constr$lhs), dimnames=list(paste0("CC2.",1:length(A)),colnames(etrib$constr$lhs)))
  
  rows <- 0
  for(a in A){
    rows <- rows + 1
    for(h in 2:length(H)){
      v <-  paste0("v(a",a,",h",h,")")
      lhs[rows, v] <- -1 * M
      lhs[rows, "L"] <- -1
      for(j in J){
          cab1 <- paste0("c",j,"(a",a,",b",h-1,")")
          lhs[rows, cab1] <- 1
      }   
    }
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  ##// TODO dir rhs
  return(etrib)
}

createCC3Constraints <- function(etrib, A, H, J){
  lhs <- matrix(0, nrow=length(A), ncol=ncol(etrib$constr$lhs), dimnames=list(paste0("CC3.",1:length(A)),colnames(etrib$constr$lhs)))
  
  rows <- 0
  for(a in A){
    rows <- rows + 1
    for(h in 1:(length(H)-1){
      v <-  paste0("v(a",a,",h",h,")")
      lhs[rows, v] <- M
      lhs[rows, "L"] <- -1
      lhs[rows, "e"] <- 1
      for(j in J){
        cab1 <- paste0("c",j,"(a",a,",b",h,")")
        lhs[rows, cab1] <- 1
      }   
    }
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  ##// TODO dir rhs
  return(etrib)
}