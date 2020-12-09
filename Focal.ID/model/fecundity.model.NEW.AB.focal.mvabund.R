##################################################
# code modified from Mayfield, M.M. and D.B. Stouffer (2017) Higher-
#    order interactions capture unexplained complexity in diverse communties.
#    Nature Ecology & Evolution 1,Article number: 0062. 
#    DOI: https://doi.org/10.1038/s41559-016-0062
##################################################

####
# function to fit the negative-binomial fecundity model
####

fecundity.model <- function(
  data,
  fit.alphas=FALSE,
  fit.betas=FALSE,
  ...
){
  library(mvabund)
  
  # add in a dummy column to separate fecundity from focal (and make it easier to identify interaction coefs)
  data$fecundity <- data$focal
  
  # lets figure out who the observed competitors are
  competitors <- colnames(data)[!colnames(data) %in% c("focal","fecundity","seeds","site","quadrat")]
  
  # start with no competition model
  base.model.formula <- "seeds ~ 0 + fecundity"
  
  # add the site fixed effect
  if(nlevels(data$site)>1){
    base.model.formula <- paste0(base.model.formula," + fecundity:site")
  }
  
  # in case it's the null
  model.formula <- base.model.formula
  
  # add pairwise coefficients for all competitor species
  all.alphas <- c()
  if(fit.alphas | fit.betas){
    # all competitors are allowed an alpha
    all.alphas <- competitors
    
    # add the alphas to the model formula
    model.formula <- paste0(model.formula, " + ", paste0(paste0("focal",":",all.alphas), collapse=" + "))
    
  }
  
  # determine potential higher order interactions if so desired and add them to the model specification
  all.betas <- c()
  if(fit.betas){
    # betas between conspecific neighbors
    possible.intrabetas <- all.alphas
    possible.intrabetas <- unlist(lapply(possible.intrabetas, function(x){paste0("I(",x," * (",x," - 1)/2)")}))
    
    # betas between heterospecific neighbors
    possible.interbetas <- combn(all.alphas,2)
    possible.interbetas <- apply(possible.interbetas,2,paste,collapse=":")
    
    # combine all betas together into a single variable
    all.betas <- c(possible.intrabetas, possible.interbetas)
    
    # add the betas to the model formula
    model.formula <- paste0(model.formula, " + ", paste0("focal:", all.betas, collapse=" + "))
    
  }
  
  
  full.model.formula <- model.formula
  m <- manyglm(
    as.formula(model.formula),
    family="negative.binomial",
    data=data,model=T,
    maxiter=100,
    maxiter2=100,
    ...)
  
  
  return(m)		
}
