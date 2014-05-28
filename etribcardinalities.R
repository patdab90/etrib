createCC1Constraints <- function(etrib, A, H){
  rnames = paste0("CC1.",A)
  lhs <- matrix(0, nrow = length(A), ncol=ncol(etrib$constr$lhs),
                dimnames = list(rnames, colnames(etrib$constr$lhs)))
  rows <- 0
  for(a in A){
    rows <- rows + 1
    v <-  paste0("v(a",a,",h",H,")")
    lhs[rows, v] <- 1
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix("==", ncol=1, nrow=rows, dimnames=list(rnames)))
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(1, ncol=1, nrow=rows, dimnames=list(rnames)))
  return(etrib)
}

createCC2Constraints <- function(etrib, A, H, J){
  rnames <- paste0("CC2.",A)
  lhs <- matrix(0, nrow=length(A), ncol=ncol(etrib$constr$lhs), dimnames=list(rnames,colnames(etrib$constr$lhs)))
  
  rows <- 0
  for(a in A){
    rows <- rows + 1
    for(h in 2:length(H)){
      v <-  paste0("v(a",a,",h",h,")")
      lhs[rows, v] <- -1 * M
      lhs[rows, "L"] <- -1
      lhs[rows, paste0("c",J,"(a",a,",b",h-1,")")] <- 1   
    }
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix(">=", ncol=1, nrow=rows, dimnames=list(rnames)))
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(-M, ncol=1, nrow=rows, dimnames=list(rnames)))
  return(etrib)
}

createCC3Constraints <- function(etrib, A, H, J){
  rnames <- paste0("CC3.",A)
  lhs <- matrix(0, nrow=length(A), ncol=ncol(etrib$constr$lhs), dimnames=list(rnames,colnames(etrib$constr$lhs)))
  
  rows <- 0
  for(a in A){
    rows <- rows + 1
    for(h in 1:(length(H)-1)){
      v <-  paste0("v(a",a,",h",h,")")
      lhs[rows, v] <- M
      lhs[rows, "L"] <- -1
      lhs[rows, "e"] <- 1
      lhs[rows, paste0("c",J,"(a",a,",b",h,")")] <- 1
    }   
  }
  
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$dir <- rbind(etrib$constr$dir, matrix("<=", ncol=1, nrow=rows, dimnames=list(rnames)))
  etrib$constr$rhs <- rbind(etrib$constr$rhs, matrix(M, ncol=1, nrow=rows, dimnames=list(rnames)))
  return(etrib)
}

createCC4Constraints <- function(etrib, A, cardinalities){
  specifiedCard <- cardinalities[cardinalities[,2] > 0,]
  ncard <- nrow(specifiedCard)
  
  rnames <- paste0("CC4.",1:ncard)
  lhs <- matrix(0, nrow = ncard, ncol=ncol(etrib$constr$lhs), 
                dimnames = list(rnames, colnames(etrib$constr$lhs)))
  
  rhs <- matrix(0, nrow = ncard, ncol = 1, dimnames = list(rnames))
  
  dir <- matrix(">=", nrow = ncard, ncol = 1, dimnames = list(rnames))
  
  for(row in 1:ncard){
    h <- cardinalities[row, ]
    lhs[row,paste0("v(a",A,",h",h[1],")")] <- 1
    rhs[row,] <- h[2]
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  
  return(etrib)
}

createCC5Constraints <- function(etrib, A, cardinalities){
  specifiedCard <- cardinalities[cardinalities[,3] > 0,]
  ncard <- nrow(specifiedCard)
  
  rnames <- paste0("CC5.",1:ncard)
  lhs <- matrix(0, nrow = ncard, ncol=ncol(etrib$constr$lhs), 
                dimnames = list(rnames, colnames(etrib$constr$lhs)))
  
  rhs <- matrix(0, nrow = ncard, ncol = 1, dimnames = list(rnames))
  
  dir <- matrix("<=", nrow = ncard, ncol = 1, dimnames = list(rnames))
  
  for(row in 1:ncard){
    h <- cardinalities[row, ]
    lhs[row,paste0("v(a",A,",h",h[1],")")] <- 1
    rhs[row,] <- h[3]
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  
  return(etrib)
}