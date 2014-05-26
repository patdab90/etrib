createEpsilonConstraint <- function(eb){
  lhs <- matrix(c(1))
  rownames(lhs) <- c("EC")
  colnames(lhs) <- c("e")
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

createB4Constraint <- function(etrib){
  etrib$constr$lhs <- etriutils.addVariables(etrib$constr$lhs, c("L"))
  
  lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list("B4", colnames(etrib$constr$lhs)))
  lhs[,"L"] <- 1
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(1, ncol=1, nrow=1, dimnames=list("B4")))
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix("<=", ncol=1, nrow=1, dimnames=list("B4")))
  
  lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list("B4", colnames(etrib$constr$lhs)))
  lhs[,"L"] <- 1
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(0.5, ncol=1, nrow=1, dimnames=list("B4")))
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix(">=", ncol=1, nrow=1, dimnames=list("B4")))
  return(etrib)
}

createB5Constraint <- function(etrib, nCrit){
  for(j in 1:nCrit){
    name <- paste0("w",j)
    lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list("B5", colnames(etrib$constr$lhs)))
    lhs[,name] <- 1
    etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
    etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(0, ncol=1, nrow=1, dimnames=list("B5")))
    etrib$constr$dir <- rbind(etrib$constr$dir, matrix(">=", ncol=1, nrow=1, dimnames=list("B5")))
    
    lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list("B5", colnames(etrib$constr$lhs)))
    lhs[,name] <- 1
    etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
    etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(1, ncol=1, nrow=1, dimnames=list("B5")))
    etrib$constr$dir <- rbind(etrib$constr$dir, matrix("<=", ncol=1, nrow=1, dimnames=list("B5")))
  }
  return(etrib)
}

createB6Constraint <- function(etrib, performances, profiles,monotonicity, th){
  nAlts <- nrow(performances)
  nCats <- nrow(profiles)
  nCrit <- ncol(performances)
  varnames <- c()
  
  for (j in 1:nCrit) {
    for (a in 1:nAlts) {
      for (b in 1:nCats) {
        varnames = c(varnames, paste0('c', j, '(a', a, ',b', b-1, ')'))
      }
    }
  }
  
  for (j in 1:nCrit) {
    for (b in 1:nCats) {    
      for (a in 1:nAlts) {
        varnames = c(varnames, paste0('c', j, '(b', b-1, ',a', a, ')'))
      }
    }
  }
  
  etrib$constr$lhs <- etriutils.addVariables(constr=etrib$constr$lhs, names=varnames)
  rownames <- paste0("B6.", seq(1:length(varnames)))
  lhs <- matrix(0, ncol=ncol(etrib$constr$lhs), nrow=length(varnames), dimnames=list(rownames, colnames(etrib$constr$lhs)))
  
  row <- 0
  for (j in 1 : nCrit) {
    for (aInd in 1 : nAlts) {
      for (bInd in 1 : nCats) {
        row = row + 1        
        indAB <- paste0('c', j, '(a', aInd, ',b', bInd-1, ')')
        lhs[row,indAB] = 1
        lhs[row,paste0("w",j)] = -1 *
          outranking(performances[aInd,j], profiles[bInd,j],
                     th[j,1], th[j,2], th[j, 3], th[j, 4], monotonicity[j])
      }
    }
  }
  
  for (j in 1:nCrit) {
    for (bInd in 1:nCats) {
      for (aInd in 1:nAlts) {
        row = row + 1        
        indBA <- paste0('c', j, '(b', bInd-1, ',a', aInd, ')')  
        lhs[row,indBA] = 1
        val <- -1 * outranking(profiles[bInd,j], performances[aInd,j],
                               th[j,1], th[j,2], th[j, 3], th[j, 4], monotonicity[j])
        lhs[row,paste0("w",j)] = val
      }
    }
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  
  rownames(lhs) <- rownames
  dir <- matrix("==", nrow=row, ncol=1)
  rownames(dir) <- rownames
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  rhs <- matrix(0, nrow=length(varnames), ncol=1)
  rownames(rhs) <- rownames
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  
  return(etrib)
}

## p is preference threshold
## q is indifference threshold
outranking <- function(x, y, q, qMult, p, pMult, ascending) {
  stopifnot(p >= 0 && q >= 0 && p >= q)
  diff <- y - x
  if (ascending == FALSE) {
    diff = 0 - diff
  }
  
  indif <- q + qMult * x
  pref <- p + pMult * x
  
  if (diff <= indif) {
    return (1)
  } else if (diff >= pref) {
    return(0)
  } else {
    return ((pref-diff) / (pref - indif))
  }
}
