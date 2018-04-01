tofactor <- function(factorVars, df){
  # Coerce variable vector to factor
  # !! OVERWRITES GLOBAL ENV !! #
  df[factorVars] <<- lapply(df[factorVars], factor)
}

addtime <- function(birth, surgery, DFSdate, lastSeen, df, dateForm= 'ymd'){
  eval(parse(text = paste0('dateFun <- lubridate::', dateForm)))
  suppressMessages(attach(df))
  age <- (dateFun(get(lastSeen)) - dateFun(get(birth)))/365.25
  cat("Minimum age:", min(age), "; Observation #", which(age == min(age)))
  OS <- dateFun(get(lastSeen)) - dateFun(get(surgery))
  cat("\nMinimum OS:", min(OS[!is.na(OS)]), "; Observation #", which(OS == min(OS[!is.na(OS)])))
  DFS <- dateFun(get(DFSdate)) - dateFun(get(surgery))
  DFS[is.na(DFS)] <- dateFun(get(lastSeen)[is.na(DFS)]) - dateFun(get(surgery)[is.na(DFS)]) 
  cat("\nMinimum DFS:", min(DFS[!is.na(DFS)]), "; Observation #", which(DFS == min(DFS[!is.na(DFS)]))) 
  detach(df)
  output <- data.frame(age, OS, DFS, df)
  return(output)
}

cpoints <- function(variables, exit, df, dir = "<", how = "Youden", healthy = 0){
  # dir = "<" -> low healthy, dir = ">" -> high healthy
  require(OptimalCutpoints)
  cutpoints <- sapply(variables, function(x){optimal.cutpoints(X = x, status = exit, tag.healthy = healthy, methods = how, direction = dir, data = df)[[1]]})
  optpoints <- data.frame(optimalCut = sapply(cutpoints, function(x){x$optimal.cutoff$cutoff}))
  rownames(optpoints) <- variables
  aucs <- sapply(cutpoints, function(x){x$measures.acc$AUC})
  colnames(aucs) <- variables
  output <- list(optpoints, aucs)
  return(output)
}

addhl <- function(variables, cutpoints, df){
  # Add binary variables based on optimal cutpoints
  # returns new data frame with vars at the end
  points <- cutpoints[[1]]
  newvars <- sapply(variables, function(x){
    ifelse(df[[x]]>points[x, ], 1, 0)
  })
  colnames(newvars) <- paste0(colnames(newvars), ".h")
  output <- data.frame(df, newvars)
  return(output)
}

vcoxu <- function(variables,  exit, df, duration = 'OS', rounding = 4, diagnostics = TRUE){
  # Vectorized univariate cox regressions
  # Pass vector of variable names, duration variable, status variable, dataframe and rounding digits
  require('survival')
  form <- paste(
    'Surv( time =', duration, ', event =', exit, ')~')
  formulae <- sapply(variables,function(x)as.formula( paste( form, x )))
  models <- lapply(formulae, function(x){
    coxph( x, data = df )})
  results <- lapply(models,function(x){
    summary(x)})
  if(diagnostics){
    print(results)
  }
  coefs <- round(unlist(sapply(results, function(x){
    x$coefficients[, 2]})), digits = rounding)
  lower <- round(unlist(sapply(results, function(x){
    x$conf.int[, 3]})), digits = rounding)
  upper <- round(unlist(sapply(results, function(x){
    x$conf.int[, 4]})), digits = rounding)
  ci <- paste(lower, '-', upper)
  pvals <- round(unlist(sapply(results, function(x){
    x$coefficients[, 5]})), digits = rounding)
  pvals <- ifelse(pvals<0.05, paste0(pvals, '*'), pvals)
  output <- data.frame("p-value univariate" = pvals, "RR" = coefs, "CI 95 perc" = ci, stringsAsFactors = F)
  return(output)
}

vcoxm <- function(variables, controls,  exit, df, duration = OS, rounding = 4, diagnostics = TRUE){
  # same as vcoxu, variables is vector of variables of interest
  # controls are included in every model
  controls <- paste(controls, collapse = " + ")
  formulae <- sapply(variables, function(x){
    as.formula(
    paste('Surv( time =', duration, ', event =', exit, ')~',
          x, '+', controls))})
  models <- lapply(formulae, function(x){
    coxph(x, data = df)})
  results <- lapply(models,function(x){summary(x)})
  if(diagnostics){
    print(results)
  }
  coefs <- round(sapply(results, function(x){
    x$coefficients[, 2]}), digits = rounding)
  lower <- as.data.frame(round(sapply(results, function(x){
    x$conf.int[, 3]}), digits = rounding))
  upper <- as.data.frame(round(sapply(results, function(x){
    x$conf.int[, 4]}), digits = rounding))
  ci <- sapply(variables, function(x){
    paste(lower[, x], '-', upper[, x])})
  pvals <- round(sapply(results, function(x){
    x$coefficients[, 5]}), digits = rounding)
  pvals <- ifelse(pvals<0.05, paste0(pvals, '*'), pvals)
  output <- data.frame("p-value univariate" = pvals, "RR" = coefs, "CI 95 perc" = ci, stringsAsFactors = F)
  rownames(output)[1] <- "Variable"
  order <- c(sapply(variables, function(x){grep(colnames(output), pattern = x)}))
  output <- output[ , order]
  return(output)
}

coxexport <- function(vcoxout, names, univariate = T, overall = T){
  # Write cox regression results to csv file for excel 
  # Pass output of vcox*, univariate T/F, overall surv T/F
  modelName <- ifelse(univariate, "Univariate", "Multivariate")
  cNames <- c("p-value" , "RR", "CI")
  if(univariate == F){
  cNames <- rep(cNames, times = length(names))
  varNames <- seq(1, length(cNames), by = 3)
  cNames[varNames] <- paste0("p-value (", names, ")")
  }
  tableName <- ifelse(overall, "Overall survival", "Disease free survival")
  vcoxout <- rbind(cNames, vcoxout)
  rownames(vcoxout)[1] <- tableName
  colnames(vcoxout) <- NULL
  print(vcoxout)
  fileName <- ifelse(univariate, "table2uni.csv", "table2multi.csv")
  write.csv(vcoxout, file = fileName)
  cat("File written to", paste0(getwd(),"/", fileName))
}
