
#' Compute sensitivity equations of a function symbolically
#' 
#' @param f named vector of type character, the functions
#' @param states Character vector. Sensitivities are computed with respect to initial
#' values of these states
#' @param parameters Character vector. Sensitivities are computed with respect to initial
#' values of these parameters
#' @param inputs Character vector. Input functions or forcings. They are excluded from
#' the computation of sensitivities.
#' @param reduce Logical. Attempts to determine vanishing sensitivities, removes their
#' equations and replaces their right-hand side occurences by 0.
#' @param secondSens Logical. Determines whether the second sensitivities are to be computed.
#' @details The sensitivity equations are ODEs that are derived from the original ODE f.
#' They describe the sensitivity of the solution curve with respect to parameters like 
#' initial values and other parameters contained in f. These equtions are also useful
#' for parameter estimation by the maximum-likelihood method. For consistency with the
#' time-continuous setting provided by \link{adjointSymb}, the returned equations contain
#' attributes for the chisquare functional and its gradient.
#' @return Named vector of type character with the sensitivity equations. Furthermore,
#' attributes "chi" (the integrand of the chisquare functional), "grad" (the integrand
#' of the gradient of the chisquare functional), "forcings" (Character vector of the 
#' additional forcings being necessare to compute \code{chi} and \code{grad}) and "yini" (
#' The initial values of the sensitivity equations) are returned.
#' @example inst/examples/example2.R
#' @example inst/examples/example3.R
#' @export
sensitivitiesSymb <- function(f, states = names(f), parameters = NULL, inputs = NULL, reduce = FALSE, secondSens = FALSE) {
  variables <- names(f)
  states <- states[!states%in%inputs]
  
  if(is.null(parameters)) {
    pars <- getSymbols(f, exclude=c(variables, inputs, "time"))
  } else {
    pars <- parameters[!parameters%in%inputs]
  }
  
  Dyf <- jacobianSymb(f, variables)
  Dpf <- jacobianSymb(f, pars)
  
  df <- length(f)
  dv <- length(variables)
  ds <- length(states)
  dp <- length(pars)
  
  # generate sensitivity variable names and names with zero entries
  mygridY0 <- expand.grid.alt(variables, states)
  mygridP <- expand.grid.alt(variables, pars)
  sensParVariablesY0 <- apply(mygridY0, 1, paste, collapse = ".")
  sensParVariablesP <- apply(mygridP, 1, paste, collapse = ".")
  
  # Write sensitivity equations in matrix form
  Dy0y <- matrix(sensParVariablesY0, ncol = ds, nrow = dv)
  Dpy <- matrix(sensParVariablesP, ncol = dp, nrow = dv)
  
  
  gl <- c(as.vector(prodSymb(matrix(Dyf, ncol = dv), Dy0y)), 
          as.vector(sumSymb(prodSymb(matrix(Dyf, ncol = dv), Dpy), matrix(Dpf, nrow = dv))))
  
  newfun <- gl
  newvariables.grid <- expand.grid.alt(variables, c(states, pars))
  newvariables <- apply(newvariables.grid, 1, paste, collapse=".")
  names(newfun) <- newvariables
  
  # Reduce the sensitivities
  vanishing <- c(sensParVariablesY0, sensParVariablesP[Dpf == "0"])
  if(reduce) {
    newfun <- reduceSensitivities(newfun, vanishing)
    is.zero.sens <- names(newfun) %in% attr(newfun,"is.zero")
  } else {
    is.zero.sens <- rep(FALSE, length(newfun))
  }
  newfun <- newfun[!is.zero.sens]
  
  
  # Append initial values
  initials <- rep(0, length(newfun))
  names(initials) <- newvariables[!is.zero.sens]
  ones <- which(apply(newvariables.grid, 1, function(row) row[1] == row[2]))
  initials[newvariables[ones]] <- 1
  
  
  # Second Sensitivities with initials
  # I have added "outputs" to newfun already because their symbols appear in the 2nd sensitivity equations and need then to be set to zero.
  attr(newfun, "outputs") <- structure(rep(0, length(which(is.zero.sens))), names = newvariables[is.zero.sens])
  if(secondSens){  
    SecondSens <- secondSensitivities(sensParVariablesY0, sensParVariablesP, dv, dp, ds, variables, pars, states, Dy0y, Dpy, Dyf, newfun, is.zero.sens)
    Secondyini <- rep(0,length(SecondSens))
    names(Secondyini) <- names(SecondSens)
  }
  
  # Compute wrss
  pars <- c(pars, states)
  
  statesD <- paste0(states, "D")
  weightsD <- paste0("weight", statesD)
  
  res <- paste0(weightsD,"*(", states, "-", statesD, ")")
  sqres <- paste0(res, "^2")
  chi <- c(chi = paste0(sqres, collapse =" + "))
  
  sensitivities <- lapply(pars, function(p) paste0(states, ".", p))
  names(sensitivities) <- pars
  
  grad <- sapply(pars, function(p) paste0(paste("2*", res, "*", sensitivities[[p]]), collapse=" + "))
  names(grad) <- paste("chi", pars, sep=".")
  grad <- replaceSymbols(newvariables[is.zero.sens], "0", grad)
  
  
  
  attr(newfun, "chi") <- chi
  attr(newfun, "grad") <- grad
  attr(newfun, "outputs") <- structure(rep(0, length(which(is.zero.sens))), names = newvariables[is.zero.sens])
  attr(newfun, "forcings") <- c(statesD, weightsD)
  attr(newfun, "yini") <- initials
  
  if (secondSens) {
    attr(newfun, "secondSens") <- SecondSens
    attr(newfun, "secondYini") <- Secondyini
    attr(newfun, "secondOutputs") <- attr(SecondSens, "secondOutputs")
  }
  
  return(newfun)
  
}

#' Compute the second sensitivities of the ODEs
#' 
#' @param sensParVariablesY0
#' @param sensParVariablesP
#' @param dv
#' @param dp
#' @param ds
#' @param variables
#' @param pars
#' @param states
#' @param Dy0y
#' @param Dpy
#' @param Dyf
#' @param newfun
#' @param is.zero.sens
#' @details Compute second sensitivity equations with respect to parameters and initial values
#' @return the second sensitivities as names character vector with attribute "secondOutputs"
#' 
secondSensitivities<-function(sensParVariablesY0, sensParVariablesP, dv, dp, ds, variables, pars, states, Dy0y, Dpy, Dyf, newfun, is.zero.sens){
  
  # Second Sensitivies of the ODE with respect to parameters d^2f/dp^2
  # g <- newfun[(dv*ds+1):length(newfun)] #first sensitivities
  # the validity of the formulas in the case of vanishing sensitivities is still given, check my notes from late nov/middle of dec 2015
  
  Dyg <- jacobianSymb(newfun[(sum(!is.zero.sens[1:(dv*ds)], na.rm=TRUE)+1):length(newfun)], variables)
  zSum1 <- as.vector(prodSymb(matrix(Dyg,ncol=dv),Dpy)) 
  
  mygridPP<-expand.grid.alt(sensParVariablesP,pars)
  sensParVariablesPP <- apply(mygridPP, 1, paste, collapse=".")
  zSum2 <- matrix(prodSymb(matrix(Dyf,ncol=dv),matrix(sensParVariablesPP, nrow = dv)), ncol=dp) # to consider reduced symmetries. 
  # Rearrange matrix to have it ordered according to the other summands and to be able to use is.zero.sens as line index for vanishing sensitivities.
  zSum2 <- as.vector(zSum2[!is.zero.sens[(dv*ds+1):length(is.zero.sens)],]) # to consider reduced symmetries
  
  zSum3 <- jacobianSymb(newfun[(sum(!is.zero.sens[1:(dv*ds)], na.rm=TRUE)+1):length(newfun)], pars) 
  
  SecondSens <- as.vector(sumSymb(sumSymb(zSum1,zSum2),zSum3))
  names(SecondSens) <- names(zSum3)  
  
  # Second Sensitivities of the ODE with respect to initial values of the states d^2f/dx0^2, see notes of 8.1.16
  Dyg1 <- jacobianSymb(newfun[1:sum(!is.zero.sens[1:(dv*ds)], na.rm=TRUE)], variables) #names of Dyg1 would be suitable if states=variables
  z1Sum1 <- as.vector(prodSymb(matrix(Dyg1,ncol=dv), Dy0y))
  
  mygridY0Y0 <- expand.grid.alt(sensParVariablesY0,states)
  sensParVariablesY0Y0 <- apply(mygridY0Y0, 1, paste, collapse=".")
  z1Sum2 <- matrix(prodSymb(matrix(Dyf,ncol=dv),matrix(sensParVariablesY0Y0, nrow = dv)), ncol=ds) 
  z1Sum2 <- as.vector(z1Sum2[!is.zero.sens[1:(dv*ds)],])
  
  SecondSens1 <- as.vector(sumSymb(z1Sum1,z1Sum2))
  # naming
  if(identical(variables,states)) {
    names(SecondSens1) <- names(Dyg1)
  } else {
    mygridSecondSens1Names <- expand.grid.alt(names(newfun[1:sum(!is.zero.sens[1:(dv*ds)], na.rm=TRUE)]), states)
    SecondSens1Names <- apply(mygridSecondSens1Names, 1, paste, collapse=".")
    names(SecondSens1) <- SecondSens1Names
  }
  
  SecondSens <- c(SecondSens,SecondSens1)
  
  
  # consider symmetries of derivative indices ie d^2x_i/(dp_j dp_k) = d^2x_i/(dp_k dp_j)
  # the vector secondSens is of the form: (ijk), with i changing first, then j, then k
  secondsenssplit <-strsplit(names(SecondSens), ".", fixed= TRUE)
  secondsenssplit.1<-unlist(lapply(secondsenssplit, function(v) v[1]))
  secondsenssplit.2<-unlist(lapply(secondsenssplit, function(v) v[2]))
  secondsenssplit.3<-unlist(lapply(secondsenssplit, function(v) v[3]))
  
  permutedjk<-paste(secondsenssplit.1,secondsenssplit.3,secondsenssplit.2, sep=".") #interchange j and k
  pairs<-unlist(lapply(permutedjk, function(v) match(v, names(SecondSens)))) #find j<->k-pairs and possible reduced Sensitivities
  
  unique.inz<-unlist(lapply((1:length(pairs)), function(s) {ifelse(s<=pairs[s], s, NA)})) # of every j-k-permutation it takes only those with j<=k ist and sets all others to NA,
  # all pairs with one reduced sensitivity are also set to NA
  
  # if a combination like "yi.pj.pj" is zero but it still comes up as a symbol in other second sensitivities, set it to zero
  secondZerosSame <- strsplit(sensParVariablesPP, ".", fixed= TRUE)
  secondZerosSame <- secondZerosSame[unlist(lapply(secondZerosSame, function(v) v[2]==v[3]))]
  secondZerosSame <- unlist(lapply(secondZerosSame, function(v) paste(v[1],v[2],v[3],sep=".")))
  secondZerosSame <- secondZerosSame[!(secondZerosSame %in% names(SecondSens))]
  
  # replace symbolic derivatives like "yi.pj.pk" with "yi.pk.pj" if the first index combination is being eliminated due to symmetry
  origSecondNames <- c(names(SecondSens[is.na(unique.inz)]),permutedjk[is.na(pairs)],secondZerosSame) #the first part are the names symmetries in indices (not zero), 
  # the second part are the names of sens with no partner(zero). their "counterpart" is set to zero in secondZerosNames(below), 
  # the last are the names derivs to the same index (zero)
  
  secondZerosNames <- names(SecondSens[pairs[is.na(unique.inz)]]) #if sens already covered by symmetry, find the index of the corresponding permutation and use its name instead
  secondZerosNames[is.na(secondZerosNames)] <- "0" #set the sensitivities where the permutation is already kicked out by "vanishing" to 0
  secondZerosNames <- c(secondZerosNames, rep_len("0", length.out = sum(is.na(pairs), length(secondZerosSame)))) #apply the rules
  
  SecondSens <- replaceSymbols(origSecondNames, secondZerosNames, SecondSens)
  
  unique.inz <- unique.inz[!is.na(unique.inz)]
  SecondSens <- SecondSens[unique.inz]
  
  
  # Second Sensitivities of the ODE with respect to initial values of the states and parameters d^2f/dx0dp, see notes of 8.1.16
  # Their second and third indices are not covered by symmetry, so they are computed after the symmetry considerations.
  z2Sum1 <- as.vector(prodSymb(matrix(Dyg,ncol=dv), Dy0y))
  
  mygridPY0 <- expand.grid.alt(sensParVariablesP,states)
  sensParVariablesPY0 <- apply(mygridPY0, 1, paste, collapse=".")
  z2Sum2 <- matrix(prodSymb(matrix(Dyf,ncol=dv),matrix(sensParVariablesPY0, nrow = dv)), ncol=ds) 
  statereducemat <- matrix(!is.zero.sens[1:(dv*ds)], nrow = dv) #to kick sensitivites that are reduced because of the state, set corresponding elements to NA
  statereducemat <- sapply((1:dv), function(s){
    matrix(rep(statereducemat[,s], dp), nrow = dv)
  })
  statereducemat <- matrix(statereducemat, nrow=dv)
  z2Sum2[!statereducemat] <- NA
  z2Sum2 <- as.vector(z2Sum2[!is.zero.sens[(dv*ds+1):length(is.zero.sens)],]) # kick sensitivites that are reduced because of the parameters, now z2Sum2 and z2Sum1 have the same indices.
  statereduceindices <- !is.na(z2Sum2) # identify leftover indices that are reduced because of states
  
  SecondSens2 <- as.vector(sumSymb(z2Sum1[statereduceindices],z2Sum2[statereduceindices]))
  
  # naming
  if(identical(variables,states)) {
    names(SecondSens2) <- names(Dyg[statereduceindices])
  } else {
    secondsenssplit <- strsplit(names(Dyg[statereduceindices]), ".", fixed= TRUE)
    secondsenssplit.3 <- unlist(lapply(secondsenssplit, function(v) v[3]))
    names(SecondSens2) <- names(Dyg[statereduceindices][secondsenssplit.3 %in% states])
  }
  
  SecondSens <- c(SecondSens,SecondSens2)
  
  
  # Add outputs, currently programmed as if symmetrical non-zeros are vanishing as well.
  mygridY0P <- expand.grid.alt(sensParVariablesY0,pars)
  sensParVariablesY0P <- apply(mygridY0P, 1, paste, collapse=".")
  secondOutputs <- c(sensParVariablesPP, sensParVariablesY0Y0, sensParVariablesPY0, sensParVariablesY0P)
  attr(SecondSens, "secondOutputs") <- secondOutputs <- structure(rep(0, length(which(!secondOutputs%in%names(SecondSens)))), names = secondOutputs[!secondOutputs%in%names(SecondSens)])
  
  # Double check if all vanishing sensitivities are excluded from the terms, if not, replace them with "0".
  firstOutputs <- names(attr(newfun, "outputs"))
  allOutputs <- c(firstOutputs, names(secondOutputs))
  allOutputs <- allOutputs[-which(allOutputs%in%origSecondNames)]
  zeros<-rep("0", length(allOutputs))
  SecondSens <- replaceSymbols(allOutputs, zeros, SecondSens)
  
  return(SecondSens)
}



#' Compute adjoint equations of a function symbolically
#' 
#' @param f Named vector of type character, the functions
#' @param states Character vector of the ODE states for which observations are available
#' @param inputs Character vector of the "variable" input states, i.e. time-dependent parameters
#' (in contrast to the forcings).
#' @param parameters Character vector of the parameters
#' @details The adjoint equations are computed with respect to the functional 
#' \deqn{(x, u)\mapsto \int_0^T \|x(t)-x^D(t)\|^2 + \|u(t) - u^D(t)\|^2 dt,}{(x, u) -> int( ||x(t) - xD(t)||^2 + ||u(t) - uD(t)||^2, dt),} 
#' where x are the states being constrained
#' by the ODE, u are the inputs and xD and uD indicate the trajectories to be best
#' possibly approached. When the ODE is linear with respect to u, the attribute \code{inputs}
#' of the returned equations can be used to replace all occurences of u by the corresponding
#' character in the attribute. This guarantees that the input course is optimal with
#' respect to the above function.
#' @return Named vector of type character with the adjoint equations. The vector has attributes
#' "chi" (integrand of the chisquare functional), "grad" (integrand of the gradient of the chisquare functional),
#' "forcings" (character vector of the forcings necessary for integration of the adjoint equations) and
#' "inputs" (the input expressed as a function of the adjoint variables).
#' @example inst/examples/example5.R
#' @export
adjointSymb <- function(f, states=names(f), parameters = NULL, inputs=NULL) {
  
  n <- length(f)
  
  adjNames <- paste0("adj", names(f))
  
  ## Compute adjoint sensitivities  
  jac <- matrix(jacobianSymb(f), n)
  
  negadj <- as.vector(prodSymb(t(jac), matrix(adjNames, ncol=1)))
  adj <- paste0("-(", negadj, ")")
  names(adj) <- adjNames
  
  
  negres <- paste0("(", states, "D - ", states, ")") 
  wres <- paste(negres, paste0("weight", states, "D"), sep="*") 
  
  
  adj[paste0("adj", states)] <- paste(adj[paste0("adj", states)], wres, sep="+")
  
  ## Input equations
  u <- NULL
  if(!is.null(inputs)) {
    
    jac <- matrix(jacobianSymb(f, inputs), n )
    u <- as.vector(
      sumSymb(paste0(inputs, "D"), 
              matrix(paste0("-(",  
                            as.vector(
                              prodSymb(t(jac), matrix(adjNames, ncol=1))),
                            ")*eps/weight", inputs, "D"),
                     ncol=1)
      )
    )
    u <-  paste0("(", u, ")")
    names(u) <- inputs
    
  }
  
  ## Forcings required by the BVP solver
  forcings <- paste0(c(states, inputs), "D")
  forcings <- c(forcings, paste0("weight", forcings))
  
  
  
  ## Time-continous log-likelihood
  res <- paste0("(", c(states, inputs), "-", c(states, inputs), "D)")
  sres <- paste0(res, "^2")
  wsres <- paste(paste0("weight", c(states, inputs), "D"),sres,  sep="*") 
  
  chi <- paste(wsres, collapse = " + ")
  names(chi) <- "chi"
  attr(adj, "chi") <- chi
  
  ## Adjoint "gradient wrt parameters" equations
  if(is.null(parameters)) {
    symbols <- getSymbols(f)
    parameters <- symbols[!(symbols%in%inputs) & !(symbols%in%names(f))]
  }
  if(length(parameters)>0) {
    jac <- matrix(jacobianSymb(f, parameters), n)
    grad <- as.vector(prodSymb(matrix(adjNames, nrow=1), jac))
    gradP <- paste0("2*(", grad, ")")
    names(gradP) <- paste("chi", parameters, sep=".")
    attr(adj, "grad") <- gradP  
  }
  
  
  
  attr(adj, "forcings") <- forcings
  attr(adj, "inputs") <- u
  
  
  
  
  return(adj)
  
}


#'@title Forcings data.frame
#'@name forcData
#'@docType data
#' 
NULL

#'@title Time-course data of O, O2 and O3
#'@name oxygenData
#'@docType data 
NULL
