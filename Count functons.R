
# gestMultiple was modified so that Ybin is always T, so it follows binary outcome procedure that suits to count outcomes as well.

gestMultipleCount = function (data, idvar, timevar, Yn, An, Cn = NA, outcomemodels, 
          propensitymodel, censoringmodel = NULL, type, EfmVar = NA, 
          cutoff = NA, ...) 
{
  if (!is.data.frame(data)) 
    (stop("Either no data set has been given, or it is not in a data frame."))
  if (is.na(EfmVar) && type %in% c(2, 4)) 
    (stop("Type 2 or 4 is specified but argument EfmVar not specified."))
  if (!is.na(EfmVar) && !is.numeric(data[, EfmVar]) && type %in% 
      c(2, 4)) 
    (stop("Effect modification is only supported for a continuous covariate, or binary covariate written as an as.numeric() 0,1 vector"))
  if (!is.na(Cn) == TRUE && !is.numeric(data[, Cn])) 
    (stop("A censoring indicator must be written as an as.numeric() 0,1 vector, with 1 indicating censoring."))
  if (!is.null(censoringmodel)) 
    (warning("Variables included in censoringmodel should ideally be included in propensitymodel else propensity scores may be invalid."))
  if (!is.factor(data[, Yn]) && !is.numeric(data[, Yn])) 
    (stop("Outcome Yn must be an as.numeric() continuous variable, or if binary, an as.numeric() 0 1 variable."))
  if (!is.factor(data[, Yn]) && !is.numeric(data[, Yn])) 
    (stop("Exposure An must be either an as.factor() categorical variable, or an as.numeric() variable. If Binary, it must be set either as a two category as.factor() variable or a numeric 0 1 variable."))
  Ybin <- TRUE
  Abin <- FALSE
  Acat <- FALSE
  if (setequal(unique(data[, An][!is.na(data[, An])]), c(0, 
                                                         1)) && is.numeric(data[, An])) 
    (Abin <- TRUE)
  if (is.factor(data[, An])) 
    (Acat <- TRUE)
  if (!is.numeric(data[, timevar])) 
    (stop("timevar must be as as.numeric() variable starting at 1"))
  if (is.na(min(data[, idvar]))) 
    (stop("idvar must not contain any missing values"))
  if (min(data[, timevar]) != 1) 
    (stop("timevar must be as as.numeric() variable starting at 1. It must also not contain any missing values"))
  if (nrow(data) != (length(unique(data[, idvar])) * max(data[, 
                                                              timevar]))) 
    (stop("There must a a row entry for each individual at each time period. For those with entries missing or censored at a time point, add\n                                                                       rows of missing values except for the time and id variable. Consider using the function FormatData."))
  T <- max(data[, timevar])
  if (is.na(cutoff) == TRUE) {
    cutoff <- T
  }
  data$int <- 1
  lmp <- formula(propensitymodel)
  if (Acat == TRUE) {
    modp <- multinom(lmp, data = data)
  }
  else if (Abin == TRUE) {
    modp <- glm(lmp, family = "binomial", data = data)
  }
  else {
    modp <- glm(lmp, family = "gaussian", data = data)
  }
  if (Acat == TRUE) {
    props <- predict(modp, type = "probs", newdata = data)
    if (nlevels(data[, An]) == 2) {
      data$prs <- props
    }
    else {
      data$prs <- props[, -1]
    }
  }
  else {
    props <- predict(modp, type = "response", newdata = data)
    data$prs <- props
  }
  cps <- NA
  if (is.na(Cn)) {
    data$w <- 1
  }
  else {
    lmc <- formula(censoringmodel)
    modc <- glm(lmc, family = "binomial", data = data)
    cps <- 1 - predict(modc, type = "response", newdata = data)
    data$cps <- cps
    data[, paste(Cn, "0", sep = "")] <- as.integer(!data[, 
                                                         Cn])
    data$cprod <- data$cps
    data[is.na(data$cprod) == TRUE, "cprod"] <- 1
    data$w <- data[, paste(Cn, "0", sep = "")]/data$cps
  }
  data$H <- data[, Yn]
  dc <- data
  dc$cntstep <- 1
  dcom <- data[complete.cases(data), ]
  for (i in 2:cutoff) {
    d2 <- data[data[, timevar] %in% seq(1, T - (i - 1), 
                                        by = 1), ]
    d2$cntstep <- i
    dc <- rbind(dc, d2)
  }
  dc <- dc[order(dc[, idvar], dc[, timevar]), ]
  if (type == 1) {
    z <- c("int")
    timevarying <- FALSE
  }
  else if (type == 2) {
    z <- c("int", EfmVar)
    timevarying <- FALSE
  }
  else if (type == 3) {
    z <- c("int")
    timevarying <- TRUE
  }
  else if (type == 4) {
    z <- c("int", EfmVar)
    timevarying <- TRUE
  }
  par1 <- paste(eval(An), eval(z), sep = ":")
  par1[par1 == paste(eval(An), "int", sep = ":")] <- paste(eval(An))
  par2 <- paste("prs", eval(z), sep = ":")
  par2[par2 == paste("prs", "int", sep = ":")] <- paste("prs")
  if (Ybin == TRUE) {
    family <- Gamma(link = "log")
  }
  else {
    family <- gaussian
  }
  for (i in 1:length(outcomemodels)) {
    outcomemodels[[i]] <- formula(outcomemodels[[i]])
    termlabs <- attr(terms(outcomemodels[[i]]), which = "term.labels")
    if (identical(as.numeric(length(which(termlabs == An))), 
                  0)) {
      stop("Every formula in outcomemodels must have an An term")
    }
    if (type %in% c(2, 4)) {
      if (identical(as.numeric(length(which(termlabs == 
                                            paste(An, EfmVar, sep = ":")))), 0)) {
        stop("For types 2 and 4. Every formula in outcomemodels must have an An term, Efm term, and\n                                 an An:Efm term. The An term must appear before any EfmVar term in each formula in outcomemodels.\n                                 Or there must be an An*EfmVar term")
      }
      if (!identical(as.numeric(length(which(termlabs == 
                                             EfmVar))), 0) && which(termlabs == EfmVar) < 
          which(termlabs == An)) 
        (stop("For types 2 and 4. Every formula in outcomemodels must have an An term, Efm term, and\n                                  an An:Efm term. The An term must appear before any EfmVar term in each formula in outcomemodels. Or there must be an An*EfmVar term"))
    }
    outcomemodels[[i]] <- reformulate(c(termlabs, par2), 
                                      response = "H")
  }
  if (timevarying == FALSE) {
    lmy <- outcomemodels[[1]]
    out <- summary(geem(terms(lmy), family = family, id = dcom[, 
                                                               idvar], data = dcom, weights = dcom$w))
    if (Acat == TRUE) {
      nam1 <- paste(An, levels(data[, An])[-1], sep = "")
      nam2 <- apply(expand.grid(nam1, z[-1]), 1, paste, 
                    collapse = ":")
      Acoef <- c(nam1, nam2)
      psi0 <- out$beta[match(Acoef, out$coefnames)]
      names(psi0) <- Acoef
      psicat <- as.list(NULL)
      for (l in 2:nlevels(data[, An])) {
        psicat[[l]] <- psi0[grep(levels(data[, An])[l], 
                                 Acoef)]
      }
      psicat[[1]] <- rep(0, length(psicat[[2]]))
    }
    else {
      psi <- out$beta[match(par1, out$coefnames)]
      names(psi) <- par1
    }
    if (Acat == TRUE) {
      i <- 2
      while (i <= cutoff && i <= T) {
        j <- 1
        dc$psiZA <- 0
        while (j <= (i - 1)) {
          if (length(z) == 1) {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, 
                                                         An])[l] & !is.na(dc[, An]), "psiZA"] <- psicat[[l]]
            }
          }
          else {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, 
                                                         An])[l] & !is.na(dc[, An]), "psiZA"] <- rowSums(sweep(dc[dc$cntstep == 
                                                                                                                    j & dc[, An] == levels(dc[, An])[l] & 
                                                                                                                    !is.na(dc[, An]), z], 2, psicat[[l]], 
                                                                                                               "*"))
            }
          }
          j <- j + 1
        }
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            if (is.na(Cn) == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "cprod"] <- dc[dc$cntstep == (j - 
                                                      1) & dc[, timevar] == k, "cps"] * dc[dc$cntstep == 
                                                                                             (j - 1) & dc[, timevar] == (k + 1), 
                                                                                           "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "w"] <- data[data[, timevar] == (k + 
                                                         j - 1), paste(Cn, "0", sep = "")]/dc[dc$cntstep == 
                                                                                                j & dc[, timevar] == k, "cprod"]
            }
            if (Ybin == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                   dc[, timevar] == (k + 1), "H"] - dc[dc$cntstep == 
                                                                         (j - 1) & dc[, timevar] == (k + 1), 
                                                                       "psiZA"]
            }
            if (Ybin == TRUE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                   dc[, timevar] == (k + 1), "H"] * exp(-dc[dc$cntstep == 
                                                                              (j - 1) & dc[, timevar] == (k + 1), 
                                                                            "psiZA"])
            }
          }
          j <- j + 1
        }
        dt <- dc[dc$cntstep %in% seq(1, i, by = 1), 
        ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH), family = family, 
                            id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi0 <- out$beta[match(Acoef, out$coefnames)]
        names(psi0) <- Acoef
        psicat <- as.list(NULL)
        for (l in 2:nlevels(data[, An])) {
          psicat[[l]] <- psi0[grep(levels(data[, An])[l], 
                                   Acoef)]
        }
        psicat[[1]] <- rep(0, length(psicat[[2]]))
        i <- i + 1
      }
      results <- list(psi = unlist(psicat[-1]), Data = as_tibble(data[, 
                                                                      !names(data) %in% c("H", "psiZA")]), PropensitySummary = summary(data$prs), 
                      CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
    else {
      i <- 2
      while (i <= cutoff && i <= T) {
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            if (is.na(Cn) == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "cprod"] <- dc[dc$cntstep == (j - 
                                                      1) & dc[, timevar] == k, "cps"] * dc[dc$cntstep == 
                                                                                             (j - 1) & dc[, timevar] == (k + 1), 
                                                                                           "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "w"] <- data[data[, timevar] == (k + 
                                                         j - 1), paste(Cn, "0", sep = "")]/dc[dc$cntstep == 
                                                                                                j & dc[, timevar] == k, "cprod"]
            }
            if (Ybin == FALSE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] - dc[dc$cntstep == 
                                                                           (j - 1) & dc[, timevar] == (k + 1), 
                                                                         An] * psi * data[data[, timevar] == 
                                                                                            (k + 1), z]
              }
              else {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] - rowSums(dc[dc$cntstep == 
                                                                                   (j - 1) & dc[, timevar] == (k + 1), 
                                                                                 An] * sweep(dc[dc$cntstep == (j - 
                                                                                                                 1) & dc[, timevar] == (k + 1), z], 
                                                                                             2, psi, "*"))
              }
            }
            if (Ybin == TRUE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] * exp(dc[dc$cntstep == 
                                                                               (j - 1) & dc[, timevar] == (k + 1), 
                                                                             An] * -psi * data[data[, timevar] == 
                                                                                                 (k + 1), z])
              }
              else {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] * exp(-rowSums(dc[dc$cntstep == 
                                                                                        (j - 1) & dc[, timevar] == (k + 1), 
                                                                                      An] * sweep(dc[dc$cntstep == (j - 
                                                                                                                      1) & dc[, timevar] == (k + 1), z], 
                                                                                                  2, psi, "*")))
              }
            }
          }
          j <- j + 1
        }
        dt <- dc[dc$cntstep %in% seq(1, i, by = 1), 
        ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH), family = family, 
                            id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi <- out$beta[match(par1, out$coefnames)]
        names(psi) <- par1
        i <- i + 1
      }
      results <- list(psi = psi, Data = as_tibble(data[, 
                                                       !names(data) %in% c("H")]), PropensitySummary = summary(data$prs), 
                      CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
  }
  else if (timevarying == TRUE) {
    if (Acat == TRUE) {
      lmy <- outcomemodels[[1]]
      out <- summary(geem(terms(lmy), family = family, 
                          id = dcom[, idvar], data = dcom, weights = dcom$w))
      nam1 <- paste(An, levels(data[, An])[-1], sep = "")
      nam2 <- apply(expand.grid(nam1, z[-1]), 1, paste, 
                    collapse = ":")
      Acoef <- c(nam1, nam2)
      psi0 <- out$beta[match(Acoef, out$coefnames)]
      names(psi0) <- Acoef
      psicatlist <- as.list(NULL)
      psicatresult <- as.list(NULL)
      psicat <- as.list(NULL)
      for (l in 2:nlevels(data[, An])) {
        psicat[[l]] <- psi0[grep(levels(data[, An])[l], 
                                 Acoef)]
      }
      psicat[[1]] <- rep(0, length(psicat[[2]]))
      psicatlist[[1]] <- psicat
      psicatresult[[1]] <- psicat[-1]
      i <- 2
      while (i <= cutoff && i <= T) {
        j <- 1
        dc$psiZA <- 0
        while (j <= (i - 1)) {
          if (length(z) == 1) {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, 
                                                         An])[l] & !is.na(dc[, An]), "psiZA"] <- psicatlist[[j]][[l]]
            }
          }
          else {
            for (l in 1:nlevels(dc[, An])) {
              dc[dc$cntstep == j & dc[, An] == levels(dc[, 
                                                         An])[l] & !is.na(dc[, An]), "psiZA"] <- rowSums(sweep(dc[dc$cntstep == 
                                                                                                                    j & dc[, An] == levels(dc[, An])[l] & 
                                                                                                                    !is.na(dc[, An]), z], 2, psicatlist[[j]][[l]], 
                                                                                                               "*"))
            }
          }
          j <- j + 1
        }
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            if (is.na(Cn) == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "cprod"] <- dc[dc$cntstep == (j - 
                                                      1) & dc[, timevar] == k, "cps"] * dc[dc$cntstep == 
                                                                                             (j - 1) & dc[, timevar] == (k + 1), 
                                                                                           "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "w"] <- data[data[, timevar] == (k + 
                                                         j - 1), paste(Cn, "0", sep = "")]/dc[dc$cntstep == 
                                                                                                j & dc[, timevar] == k, "cprod"]
            }
            if (Ybin == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                   dc[, timevar] == (k + 1), "H"] - dc[dc$cntstep == 
                                                                         (j - 1) & dc[, timevar] == (k + 1), 
                                                                       "psiZA"]
            }
            if (Ybin == TRUE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                   dc[, timevar] == (k + 1), "H"] * exp(-dc[dc$cntstep == 
                                                                              (j - 1) & dc[, timevar] == (k + 1), 
                                                                            "psiZA"])
            }
          }
          j <- j + 1
        }
        dt <- dc[dc$cntstep %in% i, ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH), family = family, 
                            id = dtcom[, idvar], data = dtcom, weights = dtcom$w))
        psi0 <- out$beta[match(Acoef, out$coefnames)]
        names(psi0) <- Acoef
        psicat <- as.list(NULL)
        for (l in 2:nlevels(data[, An])) {
          psicat[[l]] <- psi0[grep(levels(data[, An])[l], 
                                   Acoef)]
        }
        psicat[[1]] <- rep(0, length(psicat[[2]]))
        psicatlist[[i]] <- psicat
        psicatresult[[i]] <- psicat[-1]
        i <- i + 1
      }
      nam <- as.vector(NULL)
      for (p in 1:cutoff) {
        nam[p] <- paste("c=", p, sep = "")
      }
      names(psicatresult) <- nam
      results <- list(psi = unlist(psicatresult), Data = as_tibble(data[, 
                                                                        !names(data) %in% c("H", "psiZA")]), PropensitySummary = summary(data$prs), 
                      CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
    else {
      lmy <- outcomemodels[[1]]
      out <- summary(geem(terms(lmy), family = family, 
                          id = dcom[, idvar], data = dcom, weights = dcom$w))
      psi <- out$beta[match(par1, out$coefnames)]
      names(psi) <- par1
      psilist <- as.list(NULL)
      psilist[[1]] <- psi
      i <- 2
      while (i <= cutoff && i <= T) {
        j <- 2
        while (j <= i) {
          for (k in 1:(T - (j - 1))) {
            if (is.na(Cn) == FALSE) {
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "cprod"] <- dc[dc$cntstep == (j - 
                                                      1) & dc[, timevar] == k, "cps"] * dc[dc$cntstep == 
                                                                                             (j - 1) & dc[, timevar] == (k + 1), 
                                                                                           "cprod"]
              dc[dc$cntstep == j & dc[, timevar] == 
                   k, "w"] <- data[data[, timevar] == (k + 
                                                         j - 1), paste(Cn, "0", sep = "")]/dc[dc$cntstep == 
                                                                                                j & dc[, timevar] == k, "cprod"]
            }
            if (Ybin == FALSE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] - dc[dc$cntstep == 
                                                                           (j - 1) & dc[, timevar] == (k + 1), 
                                                                         An] * psilist[[j - 1]] * data[data[, 
                                                                                                            timevar] == (k + 1), z]
              }
              else {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] - rowSums(dc[dc$cntstep == 
                                                                                   (j - 1) & dc[, timevar] == (k + 1), 
                                                                                 An] * sweep(dc[dc$cntstep == (j - 
                                                                                                                 1) & dc[, timevar] == (k + 1), z], 
                                                                                             2, psilist[[j - 1]], "*"))
              }
            }
            if (Ybin == TRUE) {
              if (length(z) == 1) {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] * exp(dc[dc$cntstep == 
                                                                               (j - 1) & dc[, timevar] == (k + 1), 
                                                                             An] * -psilist[[j - 1]] * data[data[, 
                                                                                                                 timevar] == (k + 1), z])
              }
              else {
                dc[dc$cntstep == j & dc[, timevar] == 
                     k, "H"] <- dc[dc$cntstep == (j - 1) & 
                                     dc[, timevar] == (k + 1), "H"] * exp(-rowSums(dc[dc$cntstep == 
                                                                                        (j - 1) & dc[, timevar] == (k + 1), 
                                                                                      An] * sweep(dc[dc$cntstep == (j - 
                                                                                                                      1) & dc[, timevar] == (k + 1), z], 
                                                                                                  2, psilist[[j - 1]], "*")))
              }
            }
          }
          j <- j + 1
        }
        dt <- dc[dc$cntstep %in% i, ]
        dtcom <- dt[complete.cases(dt), ]
        lmH <- outcomemodels[[i]]
        out <- summary(geem(terms(lmH, keep.order = T), 
                            family = family, id = dtcom[, idvar], data = dtcom, 
                            weights = dtcom$w))
        psi <- out$beta[match(par1, out$coefnames)]
        names(psi) <- par1
        psilist[[i]] <- psi
        i <- i + 1
      }
      nam <- as.vector(NULL)
      for (p in 1:cutoff) {
        nam[p] <- paste("c=", p, sep = "")
      }
      names(psilist) <- nam
      results <- list(psi = unlist(psilist), Data = as_tibble(data[, 
                                                                   !names(data) %in% c("H")]), PropensitySummary = summary(data$prs), 
                      CensoringSummary = summary(cps))
      class(results) <- "Results"
      return(results)
    }
  }
}

# bootstrap function suitable for the above function for count data, gestfunc replaced with gestMultipleCount

gestbootCount= function ( data, idvar, timevar, Yn, An, Cn, outcomemodels, 
          propensitymodel, censoringmodel = NULL, type, EfmVar = NA, 
          cutoff, bn, alpha = 0.05, onesided = "twosided", seed = NULL, 
          ...) 
{
  t0 <- gestMultipleCount(data = data, idvar = idvar, timevar = timevar, 
                 Yn = Yn, An = An, Cn = Cn, outcomemodels = outcomemodels, 
                 propensitymodel = propensitymodel, censoringmodel = censoringmodel, 
                 type = type, EfmVar = EfmVar, cutoff = cutoff)$psi
  nams <- names(t0)
  data <- data[order(data[, idvar], data[, timevar]), ]
  Data <- data %>% nest_legacy(-all_of(idvar))
  if (!is.null(seed)) 
    set.seed(seed)
  bs <- bootstraps(Data, times = bn)
  results1 <- as.list(NULL)
  results <- as.list(NULL)
  for (j in 1:bn) {
    tryCatch({
      b <- as.data.frame(as_tibble(bs$splits[[j]]) %>% 
                           unnest_legacy())
      b[, idvar] <- rep(1:length(unique(data[, idvar])), 
                        each = max(b[, timevar]))
      results1[[j]] <- gestMultipleCount(data = b, idvar = idvar, 
                                timevar = timevar, Yn = Yn, An = An, Cn = Cn, 
                                outcomemodels = outcomemodels, propensitymodel = propensitymodel, 
                                censoringmodel = censoringmodel, type = type, 
                                EfmVar = EfmVar, cutoff = cutoff)$psi
      results[[j]] <- unlist(results1[[j]])
    }, error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    })
  }
  if (length(unlist(results)) < bn) 
    (warning("One or more bootstrapped datasets failed to obtain a fitted causal parameter. Consider removing terms from Lny to avoid collinearity, or assess the sparseness of the data."))
  mean <- colMeans(do.call(rbind, results))
  resultsmat <- do.call(rbind, results)
  ci.quant <- function(x = NA) {
    return(quantile(x, probs = c(alpha/2, 1 - alpha/2)))
  }
  ci.quant.bonf <- function(x = NA) {
    return(quantile(x, probs = c(alpha/(2 * length(unlist(t0))), 
                                 1 - alpha/(2 * length(unlist(t0))))))
  }
  ci.quant.upper <- function(x = NA) {
    return(quantile(x, probs = c(1 - alpha)))
  }
  ci.quant.bonf.upper <- function(x = NA) {
    return(quantile(x, probs = c(1 - alpha/(length(unlist(t0))))))
  }
  ci.quant.lower <- function(x = NA) {
    return(quantile(x, probs = c(alpha)))
  }
  ci.quant.bonf.lower <- function(x = NA) {
    return(quantile(x, probs = c(alpha/(length(unlist(t0))))))
  }
  results.sort <- apply(resultsmat, 2, sort)
  if (onesided == "twosided") {
    conf.quant <- t(apply(results.sort, 2, ci.quant))
    conf.quant.bonf <- t(apply(results.sort, 2, ci.quant.bonf))
  }
  else if (onesided == "upper") {
    conf.quant <- t(apply(results.sort, 2, ci.quant.upper))
    conf.quant.bonf <- t(apply(results.sort, 2, ci.quant.bonf.upper))
  }
  else if (onesided == "lower") {
    conf.quant <- t(apply(results.sort, 2, ci.quant.lower))
    conf.quant.bonf <- t(apply(results.sort, 2, ci.quant.bonf.lower))
  }
  results <- list(original = t0, mean.boot = mean, conf = conf.quant, 
                  conf.Bonferroni = conf.quant.bonf, boot.results = as_tibble(resultsmat))
  class(results) <- "Results"
  return(results)
}

# Function producing example data with count outcome. Following changes were made:
# 1) Outcome distribution was changed from normal to Poisson,
# 2) In U, L.1, L.2 and L.3, standard deviation was decreased from 1 to 0.1.
# 3) In L.1, L.2 and L.3, "1 +" was removed from the mean formula.
# The aim of latter two was to produce Y following Poisson distribution without extremely high numbers

dataexamples2 = function (n = 1000, seed = NULL, Censoring = FALSE) 
{
  if (!is.null(seed)) 
    set.seed(seed)
  n <- n
  id <- seq(1, n, by = 1)
  U <- rnorm(n, 0, 0.1)
  L.1 <- rnorm(n,  U, 0.1)
  a.1 <- plogis(1 + 0.1 * L.1)
  A.1 <- rbinom(n, 1, a.1)
  Y.1 <- rpois(n, exp( 1 + A.1 + L.1 + U))
  L.2 <- rnorm(n, (A.1/2) + L.1 + U, 0.1)
  a.2 <- plogis(1 + 0.1 * L.2 + 0.1 * A.1)
  A.2 <- rbinom(n, 1, a.2)
  Y.2 <- rpois(n, exp(1 + (A.1/2) + A.2 + L.1 + L.2 + U))
  L.3 <- rnorm(n,  (A.2/2) + L.2 + U, 0.11)
  a.3 <- plogis(1 + 0.1 * L.3 + 0.1 * A.2)
  A.3 <- rbinom(n, 1, a.3)
  Y.3 <- rpois(n, exp( 1 + (A.2/2) + A.3 + L.1 + L.2 + L.3 + U))
  Y <- rnorm(n, 1 + (A.2/2) + A.3 + L.1 + L.2 + L.3 + U, 1)
  if (Censoring == TRUE) {
    C.1 <- rbinom(n, 1, plogis(-1 + 0.001 * A.1 + 0.001 * 
                                 L.1))
    C.2 <- rbinom(n, 1, plogis(-1 + 0.001 * A.2 + 0.001 * 
                                 L.2))
    C.3 <- rbinom(n, 1, plogis(-1 + 0.001 * A.3 + 0.001 * 
                                 L.3))
    C.2[C.1 == 1] <- 1
    C.3[C.2 == 1] <- 1
    Y[C.3 == 1] <- NA
    Y.3[C.3 == 1] <- NA
    A.3[C.2 == 1] <- NA
    L.3[C.2 == 1] <- NA
    Y.2[C.2 == 1] <- NA
    A.2[C.1 == 1] <- NA
    L.2[C.1 == 1] <- NA
    Y.1[C.1 == 1] <- NA
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, 
                              L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c("A.1", 
                                                      "A.2", "A.3", "L.1", "L.2", "L.3", "C.1", "C.2", 
                                                      "C.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagest <- dl
    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, 
                              A.3, L.1, L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c("Y.1", 
                                                      "Y.2", "Y.3", "A.1", "A.2", "A.3", "L.1", "L.2", 
                                                      "L.3", "C.1", "C.2", "C.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagestmult <- dl
  }
  else {
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, 
                              L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c("A.1", 
                                                      "A.2", "A.3", "L.1", "L.2", "L.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagest <- dl
    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, 
                              A.3, L.1, L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c("Y.1", 
                                                      "Y.2", "Y.3", "A.1", "A.2", "A.3", "L.1", "L.2", 
                                                      "L.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagestmult <- dl
  }
  n <- n
  set.seed(seed)
  id <- seq(1, n, by = 1)
  U <- rnorm(n, 0, 1)
  L.1 <- rnorm(n, 1 + U, 1)
  a.1 <- plogis(1 + 0.1 * L.1)
  A.1 <- as.vector(NULL)
  for (i in 1:n) {
    A.1[i] <- sample(letters[1:3], 1, replace = TRUE, prob = c(1 - 
                                                                 (3 * a.1[i])/5, a.1[i]/5, 2 * (a.1[i]/5)))
  }
  A.1 <- as.factor(A.1)
  A.1par <- as.vector(NULL)
  A.1par[A.1 == letters[1]] <- 1
  A.1par[A.1 == letters[2]] <- 2
  A.1par[A.1 == letters[3]] <- 3
  Y.1 <- rnorm(n, 1 + A.1par + L.1 + U, 1)
  L.2 <- rnorm(n, 1 + A.1par/2 + L.1 + U, 1)
  a.2 <- plogis(1 + 0.1 * L.2 + A.1par)
  A.2 <- as.vector(NULL)
  for (i in 1:n) {
    A.2[i] <- sample(letters[1:3], 1, replace = TRUE, prob = c(1 - 
                                                                 (3 * a.2[i]/5), a.2[i]/5, 2 * (a.2[i]/5)))
  }
  A.2 <- as.factor(A.2)
  A.2par <- as.vector(NULL)
  A.2par[A.2 == letters[1]] <- 1
  A.2par[A.2 == letters[2]] <- 2
  A.2par[A.2 == letters[3]] <- 3
  Y.2 <- rnorm(n, 1 + A.1par/2 + A.2par + L.1 + L.2 + U, 1)
  L.3 <- rnorm(n, 1 + A.2par/2 + L.2 + U, 1)
  a.3 <- plogis(1 + 0.1 * L.3 + A.2par)
  A.3 <- as.vector(NULL)
  for (i in 1:n) {
    A.3[i] <- sample(letters[1:3], 1, replace = TRUE, prob = c(1 - 
                                                                 (3 * a.3[i]/5), a.3[i]/5, 2 * (a.3[i]/5)))
  }
  A.3 <- as.factor(A.3)
  A.3par <- as.vector(NULL)
  A.3par[A.3 == letters[1]] <- 1
  A.3par[A.3 == letters[2]] <- 2
  A.3par[A.3 == letters[3]] <- 3
  Y.3 <- rnorm(n, 1 + A.2par/2 + A.3par + L.1 + L.2 + L.3 + 
                 U, 1)
  Y <- rnorm(n, 1 + A.2par/2 + A.3par + L.1 + L.2 + L.3 + 
               U, 1)
  if (Censoring == TRUE) {
    C.1 <- rbinom(n, 1, plogis(-1 + 0.001 * A.1par + 0.001 * 
                                 L.1))
    C.2 <- rbinom(n, 1, plogis(-1 + 0.001 * A.2par + 0.001 * 
                                 L.2))
    C.3 <- rbinom(n, 1, plogis(-1 + 0.001 * A.3par + 0.001 * 
                                 L.3))
    C.2[C.1 == 1] <- 1
    C.3[C.2 == 1] <- 1
    Y[C.3 == 1] <- NA
    Y.3[C.3 == 1] <- NA
    A.3[C.2 == 1] <- NA
    L.3[C.2 == 1] <- NA
    Y.2[C.2 == 1] <- NA
    A.2[C.1 == 1] <- NA
    L.2[C.1 == 1] <- NA
    Y.1[C.1 == 1] <- NA
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, 
                              L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c("A.1", 
                                                      "A.2", "A.3", "L.1", "L.2", "L.3", "C.1", "C.2", 
                                                      "C.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagestcat <- dl
    datagestcat$A <- as.factor(datagestcat$A)
    levels(datagestcat$A) <- letters[1:3]
    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, 
                              A.3, L.1, L.2, L.3, C.1, C.2, C.3, U))
    dl <- reshape(dw, direction = "long", varying = c("Y.1", 
                                                      "Y.2", "Y.3", "A.1", "A.2", "A.3", "L.1", "L.2", 
                                                      "L.3", "C.1", "C.2", "C.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagestmultcat <- dl
    datagestmultcat$A <- as.factor(datagestmultcat$A)
    levels(datagestmultcat$A) <- letters[1:3]
  }
  else {
    dw <- as.data.frame(cbind(id, Y, A.1, A.2, A.3, L.1, 
                              L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c("A.1", 
                                                      "A.2", "A.3", "L.1", "L.2", "L.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagestcat <- dl
    datagestcat$A <- as.factor(datagestcat$A)
    levels(datagestcat$A) <- letters[1:3]
    dw <- as.data.frame(cbind(id, Y.1, Y.2, Y.3, A.1, A.2, 
                              A.3, L.1, L.2, L.3, U))
    dl <- reshape(dw, direction = "long", varying = c("Y.1", 
                                                      "Y.2", "Y.3", "A.1", "A.2", "A.3", "L.1", "L.2", 
                                                      "L.3"))
    dl <- dl[order(dl$id, dl$time), ]
    datagestmultcat <- dl
    datagestmultcat$A <- as.factor(datagestmultcat$A)
    levels(datagestmultcat$A) <- letters[1:3]
  }
  datas <- list(datagest = datagest, datagestmult = datagestmult, 
                datagestcat = datagestcat, datagestmultcat = datagestmultcat)
  return(datas)
}


     