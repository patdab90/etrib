createBorderProfiles <- function(prof, assig){
  ##TODO
  return(prof)
}

combineConstraintsMatrix <- function(...) {
  allConst = list(...)
  
  lhs <- c()
  dir <- c()
  rhs <- c()
  
  for (const in allConst) {
    lhs <- rbind(lhs, const$lhs)
    dir <- rbind(dir, const$dir)
    rhs <- rbind(rhs, const$rhs)
  }
  
  return(list(lhs=lhs, dir=dir, rhs=rhs))
}

getColNames <- function(nAlts, nCrit, nAssignments, nCats) {
  res <- paste('w', 1:nCrit, sep='')
  
  for (j in 1:nCrit) {
    for (a in 1:nAlts) {
      for (b in 1:nCats) {
        res = c(res, paste('c', j, '(a', a, ',b', b, ')', sep=''))
      }
    }
  }
  
  for (j in 1:nCrit) {
    for (b in 1:nCats) {    
      for (a in 1:nAlts) {
        res = c(res, paste('c', j, '(b', b, ',a', a, ')', sep=''))
      }
    }
  }
  
  for (j in 1:nCrit) {
    for (b in 1:(nCats-1)) {
      res = c(res, paste('c', j, '(b', b, ',b', (b+1), ')', sep=''))
    }
  }
  
  res = c(res, 'lam')
  res = c(res, 'e')
  
  for (a in 1:nAssignments) {
    res = c(res, paste("ass", a, "L", sep=''), paste("ass", a, "R", sep=''))
  }
  
  posInds <- c("1", "2", "3", "41", "42", "43")
  for (i in posInds) {
    res = c(res, paste("pos", i, "L", sep=''), paste("pos", i, "R", sep=''))
  }
  
  res = c(res, "necL", "necR")
  
  return(res)
}

getNrBaseVars <- function(nAlts, nCrit, nAssignments, nCats) {
  return(getEpsilonIndex(nAlts=nAlts, nCrit=nCrit, nCats=nCats)
         + 2 * nAssignments + 12 + 2)
}

getEpsilonIndex <- function(nAlts, nCrit, nCats) {
  return (getLambdaIndex(nAlts, nCrit, nCats) + 1)
}

getLambdaIndex <- function(nAlts, nCrit, nCats) {
  return (nCrit + nCrit * nAlts * nCats * 2 + nCrit * (nCats-1) + 1)
}

getWjIndex <- function(j) {
  return (j)
}

getCjBhBh1Index <- function (nCrit, nAlts, nCats, j, h) {
  stopifnot(h > 0 && h < nCats)
  stopifnot(j > 0 && j <= nCrit)
  
  offset <- nCrit +  (nCrit * nAlts * nCats * 2)
  tMinus1 <- nCats - 1
  return(offset + (j-1) * tMinus1 + h)
}

## order is j1a1b1, j1a1b2, j2a1b1, ...
getCjABIndex <- function(j, aInd, bInd, nAlts, nCats, nCrit) {
  stopifnot(bInd <= nCats && bInd > 0)
  stopifnot(aInd <= nAlts && aInd > 0)
  stopifnot(j <= nCrit && j > 0)
  
  offset <- nCrit
  index <- (j - 1) * nCats * nAlts + (aInd - 1) * nCats + bInd
  return (offset + index)
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

getCjBAIndex <- function(j, aInd, bInd, nAlts, nCats, nCrit) {
  stopifnot(bInd <= nCats && bInd > 0)
  stopifnot(aInd <= nAlts && aInd > 0)
  stopifnot(j <= nCrit && j > 0)
  
  offset <- nCrit + (nCrit * nAlts * nCats)
  return (offset + (j - 1) * nCats * nAlts + (bInd - 1) * nAlts + aInd)
}