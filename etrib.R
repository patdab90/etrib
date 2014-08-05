M <- 1E+10
MINEPS <- 1E-10

source('etriutils.R')
source('etribbase.R')
source('etribcardinalities.R')
source('etribpairwisecomp.R')


etrib.init <- function(performances, profiles, assignments, monotonicity, th, cardinalities, pCk, pCl) {
  stopifnot(ncol(performances) == ncol(profiles))
  stopifnot(is.null(assignments) || ncol(assignments) == 3)
  stopifnot((is.null(pCk)) || (ncol(pCk) == 3))
  stopifnot(is.null(pCl) || ncol(pCl) == 3)
  stopifnot(nrow(assignments) < nrow(performances))
  stopifnot(nrow(th) == ncol(performances))
  stopifnot(nrow(th) == length(monotonicity))
  
  message("--- constructing base model")
  
  n <- nrow(performances)
  p <- nrow(profiles)-1
  m <- ncol(performances)
  
  A <- 1:n ## an
  B <- 0:p ## profile
  H <- 1:p ## klasy 
  J <- 1:m ## j
  
  etrib <- list()
  etrib$n <- n
  etrib$p <- p
  etrib$m <- m
  names <- createTRIBVariables(A, H, J, pCk, pCl)
  
  varnames <- names$varnames
  etrib$binary <- names$binaryVarNames
  
  etrib$constr$lhs <- intiMatrix(varnames)
  etrib$constr$dir <- initDIRMatrix()
  etrib$constr$rhs <- intiRHSMatrix()
  
  etrib <- buildBaseModel(etrib, performances, profiles, monotonicity, th, A, H, J)
  etrib <- buildAEModel(etrib, J, assignments)
  etrib <- buildCCModel(etrib, A, H, J, cardinalities)
  etrib <- buildPConstraint(etrib, J, H, pCk, pCl)
  
  return(etrib)
}

etrib.extremeCardinalities <- function(etrib, max){
  colname <- "MAX"
  if(!max) {
    colname <- "MIN"
  }
  extremeCardinalities <- matrix(0, nrow=etrib$p, ncol=1, dimnames=list(paste0("C",1:etrib$p),colname))
  for(h in 1:etrib$p){
    varnames <- c()
    for(a in 1:etrib$n){
      binary <- paste0("v(a",a,",h",h,")")
      varnames <- c(varnames, binary) 
    }
    objFunction <- etriutils.buildObjectiveFunction(varnames=colnames(etrib$constr$lhs),objectivesVarNames=varnames)
    varTypes <- etriutils.getConstraintTypes(varnames=colnames(etrib$constr$lhs),binaryVarNames=etrib$binary) 
    ret <- etriutils.solve(objectioveFunction=objFunction,varTypes=varTypes,
                           lhs=etrib$constr$lhs,dir=etrib$constr$dir, rhs=etrib$constr$rhs,max=max)
    extremeCardinalities[h,] <- ret$objval
  }
  return(extremeCardinalities)
}

etrib.isFeasible <- function(etrib){
  objFunction <- etriutils.buildObjectiveFunction(varnames=colnames(etrib$constr$lhs),objectivesVarNames="e")
  varTypes <- etriutils.getConstraintTypes(varnames=colnames(etrib$constr$lhs),binaryVarNames=etrib$binary) 
  ret <- etriutils.solve(objectioveFunction=objFunction,varTypes=varTypes,
                         lhs=etrib$constr$lhs,dir=etrib$constr$dir, rhs=etrib$constr$rhs,max=TRUE)
  return(ret$status$code == 0 && ret$objval > 0)
}

etrib.preferenceRelation <- function(etrib){
  
  preferenceRelation <- matrix(FALSE,ncol=etrib$n, nrow=etrib$n, 
                               dimnames= list(paste0("a",1:etrib$n), paste0("a", 1:etrib$n)))
  for(a in 1:etrib$n){
    for(b in 1:etrib$n){
      plModel <- buildPLModel(colnames(etrib$constr$lhs), a, b, etrib$p, 1:etrib$m)
      lhs <- rbind(etrib$constr$lhs, plModel$lhs)
      dir <- rbind(etrib$constr$dir, plModel$dir)
      rhs <- rbind(etrib$constr$rhs, plModel$rhs)
      
      objFunction <- etriutils.buildObjectiveFunction(varnames=colnames(etrib$constr$lhs),objectivesVarNames="e")
      varTypes <- etriutils.getConstraintTypes(varnames=colnames(etrib$constr$lhs),binaryVarNames=etrib$binary) 
      ret <- etriutils.solve(objectioveFunction=objFunction,varTypes=varTypes,
                             lhs=lhs,dir=dir, rhs=rhs,max=TRUE)
      
      preferenceRelation[a,b] <- ret$status$code != 0 || ret$objval < 0
    }
  }
  return(preferenceRelation)
}

buildPLModel <- function(varnames, a, b, p, J){
  cardenality <- list(lhs = matrix(nrow=0, ncol=length(varnames)),
                      rhs = matrix(nrow=0, ncol=1), 
                      dir = matrix(nrow=0, ncol=1))
  for(h in 1:(p-1)){
    rnames <- c(paste0("PL1.",h), paste0("PL2.",h))
    lhs <- matrix(0, nrow=2, ncol=length(varnames), dimnames=list(rnames, varnames))
    lhs[1,paste0("c",J,"(a",b,",b",h,")")] = 1
    lhs[,"L"] = -1
    lhs[2,paste0("c",J,"(a",a,",b",h,")")] = 1
    lhs[2,"e"] <- 1 
    
    dir <- matrix(c(">=","<="),nrow=2, ncol=1, dimnames=list(rnames))
    
    rhs <- matrix(0,nrow=2, ncol=1, dimnames=list(rnames))
    
    cardenality$lhs <- rbind(cardenality$lhs, lhs)
    cardenality$dir <- rbind(cardenality$dir, dir)
    cardenality$rhs <- rbind(cardenality$rhs, rhs)
  }
  return(cardenality)
}

etrib.possibleAssigment <- function(etrib){
  possibleRanking <- matrix(FALSE,ncol=etrib$p, nrow=etrib$n, 
                            dimnames= list(paste0("a",1:etrib$n), paste0("C", 1:etrib$p)))
  
  for(i in 1:etrib$p){
    for(j in 1:etrib$n){
  
      paModel <- buildPAModel(colnames(etrib$constr$lhs),a=j,h=i,p=etrib$p, J=1:etrib$m)
      lhs <- rbind(etrib$constr$lhs, paModel$lhs)
      dir <- rbind(etrib$constr$dir, paModel$dir)
      rhs <- rbind(etrib$constr$rhs, paModel$rhs)
      
      objFunction <- etriutils.buildObjectiveFunction(varnames=colnames(etrib$constr$lhs),objectivesVarNames="e")
      varTypes <- etriutils.getConstraintTypes(varnames=colnames(etrib$constr$lhs),binaryVarNames=etrib$binary) 
      ret <- etriutils.solve(objectioveFunction=objFunction,varTypes=varTypes,
                           lhs=lhs,dir=dir, rhs=rhs,max=TRUE)
      
      possibleRanking[j,i] <-  ret$status$code == 0 && ret$objval > 0
    }
  }

  return(possibleRanking)
}




constraintsToString <- function(lhs, dir, rhs){
  res <- matrix("", nrow=nrow(lhs), ncol=1, dimnames=list(rownames(lhs)))
  for(j in 1:nrow(lhs)){
    for(i in 1:ncol(lhs)){
      if(lhs[j,i] != 0){
        if(lhs[j,i] > 0)
          if(lhs[j,i] == 1)
            res[j,] <- paste(res[j,],"+",colnames(lhs)[i])
        else
          res[j,] <- paste(res[j,],"+",lhs[j,i],colnames(lhs)[i])
        else
          res[j,] <- paste(res[j,],lhs[j,i],colnames(lhs)[i])
      }
    }
    res[j,] <- paste(res[j,],dir[j,],rhs[j,])
  }
  return(res)
}

buildPAModel <- function(varnames,a, h, p, J){
  cardenality <- list(lhs = matrix(nrow=0, ncol=length(varnames)),
                      rhs = matrix(nrow=0, ncol=1), 
                      dir = matrix(nrow=0, ncol=1))
  
  colnames(cardenality$lhs) <- varnames
 
  if(h > 1){
    lhs <- matrix(0, nrow=1, ncol=length(varnames), dimnames=list("PA.1", varnames))
    lhs[1,paste0("c",J,"(a",a,",b",h-1,")")] = 1
    lhs[1,"L"] = -1
  
    dir <- matrix(nrow=1, ncol=1, dimnames=list("PA.1"))
    dir[1,] <- ">="
    
    rhs <- matrix(nrow=1, ncol=1, dimnames=list("PA.1"))
    rhs[1,] <- 0
    
    cardenality$lhs <- rbind(cardenality$lhs, lhs)
    cardenality$dir <- rbind(cardenality$dir, dir)
    cardenality$rhs <- rbind(cardenality$rhs, rhs)
    
    
  }
  if(h < p){
    lhs <- matrix(0, nrow=1, ncol=length(varnames), dimnames=list("PA.2", varnames))
    lhs[1,paste0("c",J,"(a",a,",b",h,")")] = 1
    lhs[1,"L"] <- -1
    lhs[1,"e"] <- 1
   
    dir <- matrix(nrow=1, ncol=1, dimnames=list("PA.2"))
    dir[1,] <- "<="
    
    rhs <- matrix(nrow=1, ncol=1, dimnames=list("PA.2"))
    rhs[1,] <- 0
    
    cardenality$lhs <- rbind(cardenality$lhs, lhs)
    cardenality$dir <- rbind(cardenality$dir, dir)
    cardenality$rhs <- rbind(cardenality$rhs, rhs)
  }
  
  return(cardenality)
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
  if(is.null(assignments)) return(etrib)
  nAs <- nrow(assignments)
  nrows <- nAs * 2
  rnames <- paste0("AE.",1:nrows)
  lhs <- matrix(0, nrow=nrows, ncol=ncol(etrib$constr$lhs), dimnames=list(rnames,colnames(etrib$constr$lhs)))
  
  row <- 0
  for(i in 1:nAs){
    a <- assignments[i,]
    row <- row + 1
    lhs[row,paste0("c",J,"(a",a[1],",b",a[2]-1,")")] <- 1
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
  
  dir <- as.matrix(rep(c(">=","<="), each=nAs))
  rownames(dir) <- rnames
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  rhs <- as.matrix(rep(0, row))
  rownames(rhs) <- rnames
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  
  return(etrib)
}

buildCCModel <- function(etrib, A, H, J, cardinalities){
  
  etrib <- createCC1Constraints(etrib, A, H)
  etrib <- createCC2Constraints(etrib, J)
  etrib <- createCC3Constraints(etrib, J)
  if(is.null(cardinalities)) return(etrib)
  etrib <- createCC4Constraints(etrib, A, cardinalities)
  etrib <- createCC5Constraints(etrib, A, cardinalities)
  return(etrib)
}

buildPConstraint <- function(etrib, J, H, pCk, pCl){
  if(!is.null(pCk)){
    etrib <- createPC1Constrints(etrib, J, H, pCk)
    etrib <- createPC2Constrints(etrib, J, H, pCk)
    etrib <- createPC3Constrints(etrib, J, H, pCk)
  }
  if(!is.null(pCl)){
    etrib <- createPU1Constrints(etrib, J, H, pCl)
    etrib <- createPU2Constrints(etrib, J, H, pCl)
    etrib <- createPU3Constrints(etrib, J, H, pCl)
  }
  return(etrib)
}

createTRIBVariables <- function(A, H, J, pCk, pCl){
  varnames <- c()
  
  p <- length(H)
  
  varnames <- paste0("w",J)
  varnames <- c(varnames, c("L"))
  varnames <- c(varnames, c("e"))
  ##varnames = c(varnames, paste0("c", J, "(b", p, ",b0)"))

  binaryVarNames <- c()
  #if(FALSE)
  for(a in A){
    for(h in H){
      binary <- paste0("v(a",a,",h",h,")")
      varnames <- c(varnames, binary) 
      binaryVarNames <- c(binaryVarNames, binary)
    }
  }
  
  for (a in A) {
    for (b in 0:p) {
      varnames <- c(varnames, paste0('c', J, '(a', a, ',b', b, ')'))
    }
  }
  
  for (b in 0:p) {    
    for (a in A) {
      varnames <- c(varnames, paste0('c', J, '(b', b, ',a', a, ')'))
    }
  }
  if(!is.null(pCk))
  for(row in 1:nrow(pCk)){
    pck <- pCk[row,]
    binary <- paste0("v(a",pck[1],",b",pck[2],",>=",pck[3],",h",1:(p-pck[3]),")")
    varnames <- c(varnames, binary)
    binaryVarNames <- c(binaryVarNames, binary)
  }
  
  if(!is.null(pCl))
  for(row in 1:nrow(pCl)){
    pcl <- pCl[row,]
    binary <- paste0("v(a",pcl[1],",b",pcl[2],",<=",pcl[3],",h",1:(p-pcl[3]),")")
    varnames <- c(varnames, binary)
    binaryVarNames <- c(binaryVarNames, binary)
  }
  
  return(list(varnames=varnames, binaryVarNames=binaryVarNames))
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

