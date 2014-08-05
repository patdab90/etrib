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

createB1Constraint <- function(etrib, J) {
  weigthConst <- matrix(0, ncol=ncol(etrib$constr$lhs), nrow=1, dimnames=list("B1",colnames(etrib$constr$lhs)))
    
  weigths <- paste0("w",J)
  weigthConst[, weigths] <- 1
  etrib$constr$lhs <- rbind(etrib$constr$lhs, weigthConst)
  
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(1, ncol=1, nrow=1, dimnames=list("B1")))
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix("==", ncol=1, nrow=1, dimnames=list("B1")))
  
  return(etrib)
}

createB2Constraint <- function(etrib, J, H) {
  p <- length(H)
  for(j in J){
    wName <- paste0("w", j)
    lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list(paste0("B2.",j), colnames(etrib$constr$lhs)))
    lhs[,wName] <- 1
    etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
    
    etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(1/p, ncol=1, nrow=1, dimnames=list(paste0("B2.",j))))
    etrib$constr$dir <- rbind(etrib$constr$dir, matrix("==", ncol=1, nrow=1, dimnames=list(paste0("B2.",j))))
  }
  
  return(etrib)
}

createB4Constraint <- function(etrib){
  lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list("B4.1", colnames(etrib$constr$lhs)))
  lhs[,"L"] <- 1
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(1, ncol=1, nrow=1, dimnames=list("B4.1")))
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix("<=", ncol=1, nrow=1, dimnames=list("B4.1")))
  
  lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list("B4.2", colnames(etrib$constr$lhs)))
  lhs[,"L"] <- 1
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(0.5, ncol=1, nrow=1, dimnames=list("B4.2")))
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix(">=", ncol=1, nrow=1, dimnames=list("B4.2")))
  return(etrib)
}

createB5Constraint <- function(etrib, J){
  names <- paste0("w",J)
  m <- length(names)
  lhs <- matrix(0, nrow = m*2, ncol=ncol(etrib$constr$lhs),
                dimnames = list(paste0("B5.",1:(m*2)), colnames(etrib$constr$lhs)))
                
  dir <- matrix(rep.int(c(">=","<="),times=m), nrow = m*2, ncol = 1, dimnames = list(paste0("B5.",1:(m*2))))
  
  rhs <- matrix(rep.int(c(0, 1),times=m), nrow = m*2, ncol = 1, dimnames = list(paste0("B5.",1:(m*2))))
 
  row <- 0
  for(name in names){
    row <- row + 1
    lhs[row,name] <- 1
  
    row <- row + 1
    lhs[row,name] <- 1
      
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  
  return(etrib)
}

createB6Constraint <- function(etrib, performances, profiles,monotonicity, th){
  a <- nrow(performances)
  b <- nrow(profiles)-1
  m <- ncol(performances)
  
  nrows <- m*(b+1)*a*2

  rownames <- paste0("B6.", seq(1:nrows))
  lhs <- matrix(0, ncol=ncol(etrib$constr$lhs), nrow=nrows, dimnames=list(rownames, colnames(etrib$constr$lhs)))

  row <- 0
  for (j in 1:m) {
    for (aInd in 1 : a) {
      for (bInd in 0 : b) {
        row = row + 1        
        indAB <- paste0('c', j, '(a', aInd, ',b', bInd, ')')
        lhs[row,indAB] = 1
        val = outranking(performances[aInd,j], profiles[bInd+1,j],
                         th[j,2], th[j,1], th[j, 4], th[j, 3], monotonicity[j])
        lhs[row,paste0("w",j)] = -1 * val
      }
    }
  }
  
  for (j in 1:m) {
    for (bInd in 0:b) {
      for (aInd in 1:a) {
        row = row + 1        
        indBA <- paste0('c', j, '(b', bInd, ',a', aInd, ')')  
        lhs[row,indBA] = 1
        val <-  outranking(profiles[bInd+1,j], performances[aInd,j],
                               th[j,2], th[j,1], th[j, 4], th[j, 3], monotonicity[j])
        lhs[row,paste0("w",j)] = -1 * val
      }
    }
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  
  message(nrow(lhs))
  message(length(rownames))
  message(row)
  
  dir <- matrix("==", nrow=row, ncol=1)
  rownames(dir) <- rownames
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  rhs <- matrix(0, nrow=nrows, ncol=1)
  rownames(rhs) <- rownames
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  
  return(etrib)
}

outranking <- function(x, y, qb, qa, pb, pa, ascending){
  px <- pa * x + pb
  py <- pa * y + pb
  qx <- qa * x + qb
  qy <- qa * y + qb
  diff <- y - x
  if(ascending == FALSE){
    diff <- 0 - diff
  }
  if(diff >= px){
    return(0)
  }else if(diff <= qx){
    return(1)
  }else 
    return((px - diff)/(px - qx))
}

