safeFixCatDir <- function(basis.set, modelList) {
  data <- as.data.frame(modelList$data)
  
  for (i in seq_along(basis.set)) {
    var <- basis.set[[i]][2]
    
    # Skip malformed or composite terms
    if (length(var) != 1 || grepl("\\(", var)) {
      next
    }
    
    # Defensive: skip if var is not a column
    if (!var %in% colnames(data)) {
      next
    }
    
    # Safe check
    if (is.factor(data[[var]]) || is.character(data[[var]])) {
      levels.var <- levels(as.factor(data[[var]]))
      conds <- basis.set[[i]][-(1:2)]
      
      # Remove from conditioning set if present
      new.conds <- conds[!conds %in% levels.var]
      basis.set[[i]] <- c(basis.set[[i]][1:2], new.conds)
    }
  }
  
  return(basis.set)
}



basisSet2 <- function (modelList.with.data, direction = NULL, interactions = FALSE) {
  amat <- getDAG(modelList.with.data)
  b <- lapply(1:nrow(amat), function(i) {
    lapply(i:ncol(amat), function(j) {
      if (amat[i, j] != 0 | i == j)
        NULL
      else {
        cvars <- c(rownames(amat)[amat[, rownames(amat)[i],
                                       drop = FALSE] == 1], rownames(amat)[amat[,
                                                                                rownames(amat)[j], drop = FALSE] == 1])
        cvars <- cvars[!duplicated(cvars)]
        c(rownames(amat)[i], rownames(amat)[j], cvars)
      }
    })
  })
  b <- unlist(b, recursive = FALSE)
  b <- b[!sapply(b, is.null)]
  if (length(b) > 0) {
    b <- piecewiseSEM:::filterExogenous(b, modelList.with.data, amat)
    b <- piecewiseSEM:::filterSmoothed(b, modelList.with.data)
    b <- piecewiseSEM:::filterExisting(b, modelList.with.data)
    b <- piecewiseSEM:::filterInteractions(b, interactions)
    b <- piecewiseSEM:::removeCerror(b, modelList.with.data)
    b <- piecewiseSEM:::reverseAddVars(b, modelList.with.data, amat)
    b <- piecewiseSEM:::reverseNonLin(b, modelList.with.data, amat)
    # b <- piecewiseSEM:::fixCatDir(b, modelList.with.data) # Here is the cause of the problem!
    b <- safeFixCatDir(b, modelList.with.data)
    if (!is.null(direction))
      b <- piecewiseSEM:::specifyDir(b, direction)
  }
  class(b) <- "basisSet"
  return(b)
}



multigroup2 <- function(
    modelList,
    group,
    standardize = "scale",
    standardize.type = "latent.linear",
    test.type = "III") {
  name <- deparse(match.call()$modelList)
  data <- modelList$data
  
  modelList.with.data <- modelList # My solution
  modelList <- piecewiseSEM:::removeData(modelList, formulas = 1) # piecewiseSEM:::fixCatDir(b, modelList) caused the issue
  
  intModelList <- lapply(modelList, function(i) {
    rhs2 <-
      ifelse(grepl("merMod", class(i)) == T, paste(paste(paste(piecewiseSEM:::all.vars_trans(i)[-1], "*", group), collapse = " + "),  "+ (", findbars(formula(i)),")"), 
             paste(paste(piecewiseSEM:::all.vars_trans(i)[-1], "*", group), collapse = " + "))
     
     # ifelse(class(i) == "lmerMod", paste(paste(paste(piecewiseSEM:::all.vars_trans(i)[-1], "*", group), 
      #                                           collapse = " + "),  "+ (", findbars(formula(i)),")"), 
      #        paste(paste(piecewiseSEM:::all.vars_trans(i)[-1], "*", group), collapse = " + ")) 
   
     # paste(paste(piecewiseSEM:::all.vars_trans(i)[-1], "*", group),
    #       collapse = " + "
    # )
    i <- update(i, formula(paste(". ~ ", rhs2)))
    return(i)
  })
  newModelList <- lapply(unique(data[, group]), function(i) {
    update(as.psem(modelList),
           data = data[data[, group] == i, ]
    )
  })
  names(newModelList) <- unique(data[, group])
  coefsList <- lapply(
    newModelList,
    coefs,
    standardize,
    standardize.type,
    test.type
  )
  names(coefsList) <- unique(data[, group])
  coefTable <- coefs(
    modelList, standardize, standardize.type,
    test.type
  )
  anovaTable <- anova(as.psem(intModelList))[[1]]
  anovaInts <- anovaTable[grepl(":", anovaTable$Predictor), ]
  global <- anovaInts[anovaInts$P.Value >= 0.05, c(
    "Response",
    "Predictor"
  )]
  global$Predictor <-
    sub(":", "\1", sub(group, "\1", global$Predictor))
  if (nrow(global) == nrow(anovaInts)) {
    newCoefsList <- list(global = coefTable)
  } else {
    newCoefsList <- lapply(names(coefsList), function(i) {
      ct <- as.matrix(coefsList[[i]])
      idx <- which(
        apply(ct[, 1:2], 1, paste, collapse = "") %in%
          apply(global[, 1:2], 1, paste, collapse = "")
      )
      ct[idx, ] <- as.matrix(coefTable[idx, ])
      ct <- cbind(ct, ifelse(1:nrow(ct) %in% idx, "c",
                             ""
      ))
      for (j in 1:nrow(ct)) {
        if (ct[j, ncol(ct)] == "c") {
          model <- modelList[[which(sapply(
            piecewiseSEM:::listFormula(modelList),
            function(x) {
              piecewiseSEM:::all.vars.merMod(x)[1] == ct[
                j,
                "Response"
              ]
            }
          ))]]
          data. <- data[data[, group] == i, ]
          sd.x <-
            piecewiseSEM:::GetSDx(model, modelList, data., standardize)
          sd.x <- sd.x[which(names(sd.x) == ct[j, "Predictor"])]
          sd.y <- piecewiseSEM:::GetSDy(
            model, data., standardize,
            standardize.type
          )
          new.coef <- as.numeric(ct[j, "Estimate"]) *
            (sd.x / sd.y)
          ct[j, "Std.Estimate"] <- ifelse(length(new.coef) >
                                            0, round(as.numeric(new.coef), 4), "-")
        }
      }
      ct <- as.data.frame(ct)
      ct[is.na(ct)] <- "-"
      names(ct)[(ncol(ct) - 1):ncol(ct)] <- ""
      return(ct)
    })
    names(newCoefsList) <- names(coefsList)
  }
  
  if (nrow(global) == nrow(anovaInts)) {
    gof <- piecewiseSEM::fisherC(modelList)
  } else {
    # b <- piecewiseSEM:::basisSet(modelList) # Here is the cause of the problem!
    b <- basisSet2(modelList.with.data)
    cf <- coefTable[coefTable$Response %in% global$Response & coefTable$Predictor %in% global$Predictor, ]
    b <- lapply(
      X = b,
      FUN = function(i) {
        for (j in 3:length(i)) {
          value <- cf[cf$Response == i[2] & cf$Predictor == i[j], "Estimate"]
          if (length(value) != 0) {
            i[j] <- paste0("offset(", value, "*", i[j], ")")
          }
        }
        return(i)
      }
    )
    if (length(b) == 0) {
      b <- NULL
    }
    gof <- fisherC(modelList, basis.set = b)
  }
  
  ret <- list(
    name = name,
    group = group,
    Cstat = gof,
    global = global,
    anovaInts = anovaInts,
    group.coefs = newCoefsList
  )
  class(ret) <- "multigroup.psem"
  return(ret)
}




# additional lavaan functions ####

############## #
#### SET UP ####
############## #

library(lavaan)

################# #
#### FUNCTIONS ####
################# #

# get a list of the constrained paths (Serainas  doing)
get_constrained_paths <- function(fit) {
  # Extract parameter table
  pt <- parameterTable(fit)
  
  # Keep only rows with non-empty labels
  labeled <- pt[pt$label != "", ]
  
  # Find labels that are used more than once (i.e., constrained parameters)
  constraint_labels <- labeled$label[duplicated(labeled$label) | duplicated(labeled$label, fromLast = TRUE)]
  
  # Filter rows that have those constraint labels
  constrained <- labeled[labeled$label %in% constraint_labels, ]
  
  # Create a human-readable path string
  constrained$path <- with(constrained, paste(lhs, op, rhs))
  
  # Group by label
  path_list <- split(constrained$path, constrained$label)
  
  return(path_list)
}

##### path extraction function by Eric Allan ####
## a function to extract all paths from a fitted mutligroup and multilevel model, these are output in the form of a list with the paths at each level
## directed paths and covariances are returned but not variances or intercepts
extract_paths <- function(fit, intercepts = FALSE) {
  
  pp <- parTable(fit)
  
  ## select directed paths and covariances
  pp2 <- pp[pp$op=="~"|pp$op=="~~",]
  ## remove non free parameters (exogenous covariances and variances)
  pp3 <- pp2[pp2$user==1, ]
  
  lhs <- pp3[,"lhs"]
  rhs <- pp3[,"rhs"]
  ops <- pp3[,"op"]
  
  
  ## check for existing constraints
  if("label"%in%names(pp3)){
    
    ## add the label for constrained paths
    label <- pp3[,"label"]
    cons <- rhs[which(label!="")]
    conlab <- label[which(label!="")]
    cons2 <- paste(conlab, cons, sep ="*")
    
    rhs[which(label!="")] <- cons2
  }
  
  else{}
  ## make a vector of all paths, these are repeated across groups
  all.paths <- paste(lhs, ops, rhs, sep="")
  
  all.groups <- pp3[,"group"] ## list of groups for each path
  
  #unique.paths <- all.paths[all.groups==1] ## just select paths for first group (must be same as in all other groups)
  
  # all.levels <- pp3[,"level"]
  # levels <- all.levels[all.groups==1]
  
  path.list <- split(all.paths,all.groups) ## make a list with paths at different levels
  
  # path.list2 <- lapply(path.list, split, levels)
  
  if(intercepts == TRUE){
    ii <- pp[pp$op=="~1",] ## extract matrix with intercepts
    
    lhsi <- ii[,"lhs"]
    # lvii <- ii[,"level"]
    grii <- ii[,"group"]
    rhsi <- rep("~1", length(lhsi))
    
    ## check for existing constraints
    if("label"%in%names(ii)){
      
      ## add the label for constrained paths
      label <- ii[,"label"]
      conlab <- label[which(label!="")]
      cons2 <- paste("~", conlab, "*1", sep ="")
      
      rhsi[which(label!="")] <- cons2
    }
    
    else{}
    
    ## remake the intercepts as terms
    ints <- paste(lhsi,rhsi, sep ="")
    
    ## now add to the path list
    for(j in unique(all.groups)){
      
      path.list[[j]][[1]] <- c(path.list[[j]][[1]], ints[lvii==1&grii==j])
      path.list[[j]][[2]] <- c(path.list[[j]][[2]], ints[lvii==2&grii==j])  
    }
    
  }
  else{}
  
  # Return the unique paths
  return(path.list)
}

##### single path constraining function by Eric Allan, adapted by Seraina Cappelli ####

## a function to constrain a single path in a multigroup model. You need to supply the fitted model, the path to constrain and the level at which the path is fitted
## you can give the constrained path a new name with coef_name and specify which groups to constrain, e.g. c(1,2) to constrain path between groups 1 and 2. 
## output is a new model string with the constraint, you then need to fit this model with lavaan or sem function

constrain_path <- function(fit, path, coef_name, groups_constrain, intercepts = FALSE) {
  ## first extract the paths
  model_paths <- extract_paths(fit, intercepts = intercepts)
  model_pathsc <- model_paths
  
  ## split the path to constrain
  
  pt2 <- unlist(strsplit(path, "(~{1,2})"))
  if (grepl("~~", path)) {
    new_path <- paste(pt2[1], "~~", coef_name, "*", pt2[2], sep = "")
  } else {
    new_path <- paste(pt2[1], "~", coef_name, "*", pt2[2], sep = "")
  }
  
  # pt2 <- unlist(strsplit(path,"~"))
  # new_path <- paste(pt2[1],"~",coef_name, "*", pt2[2],sep="")
  
  ## list for model string
  
  con_ms <- list()
  
  ## add the constraint for selected groups
  for(i in groups_constrain){
    mp <- model_pathsc[[i]]
    
    con_paths <- gsub(pattern = paste("^",path,"$",sep=""), replacement = new_path, x = mp) 
    
    
    ## back into model string format
    
    ## add in the group text
    msg <- paste("Group: ", i, "\n", paste(con_paths, collapse = "\n"), sep = "")
    
    con_ms[i] <- msg
    
  }
  ## extract the number of groups
  total_groups <- lavInspect(fit,"ngroups")
  gg <- 1:total_groups
  gg2 <- gg[-groups_constrain]
  
  ## make original model string without constraint (for groups that should not be constrained)
  uncon_ms <- list()
  
  for(g in gg2){
    mso <- paste(model_paths[[g]], collapse = "\n") 
    
    ## add in the group text
    msgo <- paste("Group: ", g, "\n", mso, sep = "")
    
    uncon_ms[g] <- msgo
    
  }
  
  ## final string
  fin_ms <- rep(list(1),total_groups)
  
  fin_ms[groups_constrain] <- con_ms[groups_constrain]
  fin_ms[gg2] <- uncon_ms[gg2]
  
  final_model_string <- paste(fin_ms, collapse = "\n")
  
  return(final_model_string)
}

##### constraining test function adapted from Eric Allan, adapted by Seraina Cappelli ####
# note: Eric did a step wise approach to the constraining: first test each path
# whether it could be constrained in a fully unconstrained SEM. The statistics
# of these tests are stored in the table model_compare. I usually constrain my 
# models based on this. Eric decided to take a second step: Sequentially constrain 
# first the path with the lowest, then second lowest, then third lowest, etc. 
# p-value in the model_compare table, each time check again if constraining this 
# path in a partially constrained model would make a worse model fit and only 
# constrain the path if it can also be safely constrained in an already partially 
# constrained model. Both approaches are valid. I will stick with the first.
# I might have to modify this function or I will manually run it down to the 
# model_compare table.

## a function to stepwise constrain all paths in a multigroup (and multilevel) SEM. 
## It starts by constraining each path in turn and calculates p values. Then it starts with the least significant path and constrains it then proceeds to the next and so on
## you need to specify the fitted model and the dataset to which it was fitted
## output is the fitted best model, the model string for the best model, 
## a list with all models fitted in the stepwise reduction and p-values and Chi squared for single path constraint from full model

test_constraints <- function(fit, groups_constrain) {
  # Fit the configural model (no constraints)
  fit_configural <- fit
  group_var <- lavInspect(fit,"group")
  data <- lavInspect(fit, "data")
  for(name in names(data)){
    data[[name]] <- data.frame(data[[name]])
    data[[name]]$group <- name
  }
  data <- bind_rows(data, .id = "group")
  colnames(data)[colnames(data) == "group"] <- group_var
  
  constrained.paths <- get_constrained_paths(fit)
  
  path.list <- extract_paths(fit_configural)
  path.list <- path.list[[1]] ## just take the first group, here we therefore assume that no constraining has yet been done 
  path.list <- setdiff(path.list, constrained.paths)
  
  
  
  # Initialize a list to store the fit of constrained models
  fit_list <- list()
  model_compare <- data.frame(Path = unlist(path.list), ChiSquare = NA, pValue = NA, stringsAsFactors = FALSE)
  
  groups <- groups_constrain
  
  
  # First, constrain each path individually and store the results, do this for levels separately
  for(i in 1:length(path.list)){
    
    path <- path.list[[i]]
    
    coef_name <- paste("coef", i, sep = "_")
    
    constrained_model_string <- constrain_path(fit_configural, path, coef_name, groups_constrain = groups)
    fit_constrained <- sem(constrained_model_string, data = data, group = group_var)
    fit_list[[path]] <- fit_constrained
    
    compare <- anova(fit_configural, fit_constrained)
    chi_square <- compare[2, "Chisq diff"]
    p_value <- compare[2, "Pr(>Chisq)"]
    
    model_compare[which(model_compare$Path == path), "ChiSquare"] <- chi_square
    model_compare[which(model_compare$Path == path), "pValue"] <- p_value
  }
  
  return(model_compare)
  
}


