etriutils.addVariables <- function(constr, names){
  ## 
  w <- matrix(0, ncol=length(names), nrow=nrow(constr))
  colnames(w) <- names
  constr <- cbind(constr, w)
  return(constr)
}

etriutils.getVariablesCount <- function(constr){
  return(ncol(constr))
}

etriutils.solve <- function(objectioveFunction, varTypes, lhs, dir, rhs, max){
  obj <- L_objective(objectioveFunction)
  roiConst <- L_constraint(lhs, dir, rhs)
  
  lp <- OP(objective=obj, constraints=roiConst, maximum=max, types=varTypes)
  
  ret <- ROI_solve(lp, .solver)
  
  return(ret)
}

etriutils.getConstraintTypes <- function(varnames, binaryVarNames){
  types <- matrix("C", nrow = 1, ncol=length(varnames))
  colnames(types) <- varnames
  types[,binaryVarNames] <- "B"
  return(types[1,])
}

etriutils.buildObjectiveFunction <- function(varnames, objectivesVarNames){
  row <- matrix(0, ncol=length(varnames))
  colnames(row) <- varnames
  row[,objectivesVarNames] <- 1
  return(row[1,])
}