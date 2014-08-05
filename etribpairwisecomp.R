createPC1Constrints <- function(etrib, J, H, pCk){
  nrows <- length(H)*nrow(pCk)
  p <- length(H)
  
  row <- 0
  for(row in 1:nrow(pCk)){
    pc <- pCk[row,]
    for(h in 1:(p-pc[3])){
      row <- row + 1
      lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list(paste0("PC1.", row), colnames(etrib$constr$lhs)))
      dir <- matrix(">=", nrow=1, ncol=1, dimnames=list(paste0("PC1.", row)))
      rhs <- matrix(-M, nrow=1, ncol=1, dimnames=list(paste0("PC1.", row)))
      
      canames <- paste0("c",J,"(a",pc[1],",b",(h-1+pc[3]),")")
      lhs[,canames] <- 1
      lhs[,"L"] <- -1
      lhs[, paste0("v(a",pc[1],",b",pc[2],",>=",pc[3],",h",h,")")] <- -M
      
      etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
      etrib$constr$dir <- rbind(etrib$constr$dir, dir)
      etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
    }
  }
  return(etrib)
}

createPC2Constrints <- function(etrib,J, H, pCk){
  nrows <- length(H)*nrow(pCk)
  p <- length(H)
  
  row <- 0
  for(row in 1:nrow(pCk)){
    pc <- pCk[row,]
    for(h in 1:(p-pc[3])){
      row <- row + 1
      lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs),
                    dimnames=list(paste0("PC2.", row), colnames(etrib$constr$lhs)))
      dir <- matrix("<=", nrow=1, ncol=1, dimnames=list(paste0("PC2.", row)))
      rhs <- matrix(M, nrow=1, ncol=1, dimnames=list(paste0("PC2.", row)))
      
      canames <- paste0("c",J,"(a",pc[1],",b",h,")")
      lhs[,canames] <- 1
      lhs[,"L"] <- -1
      lhs[,"e"] <- 1
      lhs[, paste0("v(a",pc[1],",b",pc[2],",>=",pc[3],",h",h,")")] <- M
      
      etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
      etrib$constr$dir <- rbind(etrib$constr$dir, dir)
      etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
    }
  }
  return(etrib)
}

createPC3Constrints <- function(etrib, J, H, pCk){
  p <- length(H)
  lhs <- matrix(0, nrow=nrow(pCk), ncol=ncol(etrib$constr$lhs), 
                dimnames=list(paste0("PC3.", 1:nrow(pCk)), colnames(etrib$constr$lhs)))
  dir <- matrix(">=", nrow=nrow(pCk), ncol=1, 
                dimnames=list(paste0("PC3.", 1:nrow(pCk))))
  rhs <- matrix(1, nrow=nrow(pCk), ncol=1, 
                dimnames=list(paste0("PC3.", 1:nrow(pCk))))
  
  for(row in 1:nrow(pCk)){
    pc <- pCk[row,]    
    lhs[row, paste0("v(a",pc[1],",b",pc[2],",>=",pc[3],",h",1:(p-pc[3]),")")] <- 1        
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  return(etrib)
}








createPU1Constrints <- function(etrib, J, H, pCl){
  nrows <- length(H)*nrow(pCl)
  p <- length(H)
  
  row <- 0
  for(row in 1:nrow(pCl)){
    pc <- pCl[row,]
    for(h in 1:(p-pc[3])){
      row <- row + 1
      lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs), dimnames=list(paste0("PU1.", row), colnames(etrib$constr$lhs)))
      dir <- matrix("<=", nrow=1, ncol=1, dimnames=list(paste0("PU1.", row)))
      rhs <- matrix(M, nrow=1, ncol=1, dimnames=list(paste0("PU1.", row)))
      
      canames <- paste0("c",J,"(a",pc[1],",b",(h+pc[3]),")")
      lhs[,canames] <- 1
      lhs[,"L"] <- -1
      lhs[,"e"] <- 1
      lhs[, paste0("v(a",pc[1],",b",pc[2],",<=",pc[3],",h",h,")")] <- M
      
      etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
      etrib$constr$dir <- rbind(etrib$constr$dir, dir)
      etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
    }
  }
  return(etrib)
}

createPU2Constrints <- function(etrib,J, H, pCl){
  nrows <- length(H)*nrow(pCl)
  
  p <- length(H)
  
  row <- 0
  for(row in 1:nrow(pCl)){
    pc <- pCl[row,]
    for(h in 1:(p-pc[3])){
      row <- row + 1
      lhs <- matrix(0, nrow=1, ncol=ncol(etrib$constr$lhs),
                    dimnames=list(paste0("PU2.", row), colnames(etrib$constr$lhs)))
      dir <- matrix("<=", nrow=1, ncol=1, dimnames=list(paste0("PU2.", row)))
      rhs <- matrix(-M, nrow=1, ncol=1, dimnames=list(paste0("PU2.", row)))
      
      canames <- paste0("c",J,"(a",pc[1],",b",h-1,")")
      lhs[,canames] <- 1
      lhs[,"L"] <- -1
      lhs[, paste0("v(a",pc[1],",b",pc[2],",<=",pc[3],",h",h,")")] <- -M
      
      etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
      etrib$constr$dir <- rbind(etrib$constr$dir, dir)
      etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
    }
  }
  return(etrib)
}

createPU3Constrints <- function(etrib, J, H, pCl){
  p <- length(H)
  lhs <- matrix(0, nrow=nrow(pCl), ncol=ncol(etrib$constr$lhs), 
                dimnames=list(paste0("PU3.", 1:nrow(pCl)), colnames(etrib$constr$lhs)))
  dir <- matrix(">=", nrow=nrow(pCl), ncol=1, 
            dimnames=list(paste0("PU3.", 1:nrow(pCl))))
  rhs <- matrix(1, nrow=nrow(pCl), ncol=1, 
            dimnames=list(paste0("PU3.", 1:nrow(pCl))))
  
  for(row in 1:nrow(pCl)){
    pc <- pCl[row,]    
    lhs[row, paste0("v(a",pc[1],",b",pc[2],",<=",pc[3],",h",1:(p-pc[3]),")")] <- 1        
  }
  
  etrib$constr$lhs <- rbind(etrib$constr$lhs, lhs)
  etrib$constr$dir <- rbind(etrib$constr$dir, dir)
  etrib$constr$rhs <- rbind(etrib$constr$rhs, rhs)
  return(etrib)
}