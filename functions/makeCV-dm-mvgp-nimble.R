##
##
##
##
##
##
makeCV <- function (i, model_name=model_name, y_cv=y_cv, y_cv_prop=y_cv_prop, 
                    X_cv=X_cv, params=params, folds=folds, data_source = "pollen") {
  library(rioja)
  library(analogue)
  library(randomForest)
  library(coda)
  library(gjam)
  library(here)
  
  ## setup cross-validation training and test data
  idx_test <- which(folds == i, arr.ind=TRUE)
  y_train <- as.matrix(y_cv[-idx_test, ])
  y_train_prop <- as.matrix(y_cv_prop[-idx_test, ])
  y_test <- as.matrix(y_cv[idx_test, ])
  y_test_prop <- as.matrix(y_cv_prop[idx_test, ])
  X_train <- c(X_cv[-idx_test])
  X_test <- c(X_cv[idx_test])
  
  if (model_name=="MVGP") {
    ## Fit MVGP model
    Rcpp::sourceCpp(here("mcmc", "mcmc-dirichlet-multinomial-mvgp.cpp"))
    out <- mcmc(mcmcRcpp(y_train, X_train, y_test, params, n_chain=i, 
                         file_name=here("model-fit", "progress", 
                                        "cross-validate", "dm-cv-mvgp.txt")))
    Rcpp::sourceCpp(here("functions", "makeCRPS.cpp"))
    CRPS <- makeCRPS(out$X, X_test, 
                     params$n_mcmc/params$n_thin)
    X_mean <- apply(out$X, 2, mean)
    MSPE <- (X_mean - X_test)^2
    MAE  <- abs(apply(out$X, 2, median) -  X_test)
    X_025 <- apply(out$X, 2, quantile, prob = 0.025)
    X_975 <- apply(out$X, 2, quantile, prob = 0.975)
    coverage <- (X_test >= X_025) & (X_test <= X_975)
    rm(out)
  } else  if (model_name=="GAM") {
    ## Fit GAM model
    Rcpp::sourceCpp(here("mcmc", "mcmc-dm-basis.cpp"))
    out <- mcmc(mcmcRcpp(y_train, X_train, y_test, params, n_chain=i, 
                         file_name=here("model-fit", "progress", 
                                        "cross-validate", "dm-cv-basis.txt")))
    
    Rcpp::sourceCpp(here("functions", "makeCRPS.cpp"))
    CRPS <- makeCRPS(out$X, X_test, 
                     params$n_mcmc/params$n_thin)
    X_mean <- apply(out$X, 2, mean)
    MSPE <- (X_mean - X_test)^2
    MAE  <- abs(apply(out$X, 2, median) -  X_test)
    X_025 <- apply(out$X, 2, quantile, prob = 0.025)
    X_975 <- apply(out$X, 2, quantile, prob = 0.975)
    coverage <- (X_test >= X_025) & (X_test <= X_975)
    rm(out)
  } else if (model_name=="WA") {
    ## WA reconstruction - subset to deal with all zero occurrence species
    zeros_idx <- which(colSums(y_train_prop) == 0)
    if (length(zeros_idx) > 0) {
      modWA <- rioja::WA(y_train_prop[, - zeros_idx], X_train)
      predWA <- predict(modWA, y_test_prop[, - zeros_idx], sse=TRUE, nboot=1000)
    } else {
      ## no data to subset
      modWA <- rioja::WA(y_train_prop, X_train)
      predWA <- predict(modWA, y_test_prop, sse=TRUE, nboot=1000)  
    }
    n_train <- nrow(y_train_prop)
    n_test <- nrow(y_test_prop)
    n_boot <- 1000
    predWABoot <- matrix(0, n_test, n_boot)
    for (i in 1:n_boot) {
      s <- sample(1:n_train, n_train, replace=TRUE)
      y_train_boot <- y_train_prop[s, ]
      X_train_boot <- X_train[s]
      zeros_idx <- which(colSums(y_train_boot) == 0)
      if (length(zeros_idx) > 0) {
        modWABoot <- rioja::WA(y_train_boot[, - zeros_idx], X_train_boot)     
        predWABoot[, i] <- predict(modWABoot, y_test_prop[, - zeros_idx], sse=FALSE, nboot=1)$fit[, 1]
      } else {
        modWABoot <- rioja::WA(y_train_boot, X_train_boot)     
        predWABoot[, i] <- predict(modWABoot, y_test_prop, sse=FALSE, nboot=1)$fit[, 1]
      }
    }
    source(here("functions", "makeCRPSGauss.R"))
    # CRPS <- makeCRPSGauss(predWA$fit[, 1], 
    #                       sqrt(predWA$v1.boot[, 1]^2 + predWA$v2.boot[1]^2), X_test)
    CRPS <- abs(apply(predWABoot, 1, median) - X_test)
    # MAE <- abs(predWA$fit[, 1] - X_test)
    MAE <- abs(apply(predWABoot, 1, median) - X_test)
    MSPE <- (predWA$fit[, 1] - X_test)^2
    coverage <- (X_test >= 
                   (predWA$fit[, 1] - 2*sqrt(predWA$v1.boot[, 1]^2 + predWA$v2.boot[1]^2)) & 
                   (X_test <= (predWA$fit[, 1] + 2*sqrt(predWA$v1.boot[, 1]^2 + predWA$v2.boot[1]^2))))
  } else  if (model_name=="MAT") {
    
    ## Modern analogue technique
    modMAT <- rioja::MAT(as.data.frame(y_train_prop), X_train, k=20, lean=FALSE)
    predMAT <- predict(modMAT, as.data.frame(y_test_prop), k=20, sse=TRUE, n.boot=1000)
    n_train <- nrow(y_train_prop)
    n_test <- nrow(y_test_prop)
    n_boot <- 1000
    predMATBoot <- matrix(0, n_test, n_boot)
    for (i in 1:n_boot) {
      s <- sample(1:n_train, n_train, replace=TRUE)
      y_train_boot <- y_train_prop[s, ]
      X_train_boot <- X_train[s]
      modMATBoot <- MAT(as.data.frame(y_train_boot), X_train_boot, k=20, lean=FALSE)
      predMATBoot[, i] <- predict(modMATBoot, newdata=as.data.frame(y_test_prop), k=20, nboot=1)$fit[, 1]
    }
    
    source(here("functions", "makeCRPSGauss.R"))
    # CRPS <- makeCRPSGauss(
    #   predMAT$fit.boot[, 2], 
    #   sqrt(predMAT$v1.boot[, 2]^2+ predMAT$v2.boot[2]),X_test)
    CRPS <- abs(apply(predMATBoot, 1, median) - X_test)
    MSPE <- ( predMAT$fit.boot[, 2] - X_test)^2
    # MAE <- abs(   predMAT$fit.boot[, 2] - X_test)
    MAE <- abs(apply(predMATBoot, 1, median) - X_test)
    coverage <- 
      ( X_test >= ( predMAT$fit.boot[, 2] -
                      2 * sqrt(predMAT$v1.boot[, 2]^2+ predMAT$v2.boot[2])) & 
          (X_test <= (predMAT$fit.boot[, 2] +
                        2*  sqrt(predMAT$v1.boot[, 2]^2+ predMAT$v2.boot[2]))))
    
  } else if (model_name=="MLRC") {
    ## MLRC reconstruction - subset to deal with all zero occurrence species
    zeros_idx <- which(colSums(y_train_prop) == 0)
    if (length(zeros_idx) > 0) {
      modMLRC <- rioja::MLRC(y_train_prop[, - zeros_idx], X_train)
      predMLRC <- predict(modMLRC, y_test_prop[, - zeros_idx],
                          sse=TRUE, nboot=1000)
    } else {
      modMLRC <- rioja::MLRC(y_train_prop, X_train)
      predMLRC <- predict(modMLRC, y_test_prop, sse=TRUE, nboot=1000)
    }
    n_train <- nrow(y_train_prop)
    n_test <- nrow(y_test_prop)
    n_boot <- 1000
    predMLRCBoot <- matrix(0, n_test, n_boot)
    for (i in 1:n_boot) {
      s <- sample(1:n_train, n_train, replace=TRUE)
      y_train_boot <- y_train_prop[s, ]
      zeros_idx <- which(colSums(y_train_boot) == 0)
      X_train_boot <- X_train[s]
      if (length(zeros_idx) > 0) {
        modMLRCBoot <- MLRC(y_train_boot[, - zeros_idx], X_train_boot, n.cut=0.1)
        predMLRCBoot[, i] <- predict(modMLRCBoot, newdata=y_test_prop[, - zeros_idx], sse=TRUE,
                                     nboot=1, verbose = FALSE)$fit.boot
      } else {
        modMLRCBoot <- MLRC(y_train_boot, X_train_boot, n.cut=0.1)
        predMLRCBoot[, i] <- predict(modMLRCBoot, newdata=y_test_prop, sse=TRUE,
                                     nboot=1, verbose = FALSE)$fit.boot
      }
    }
    source(here("functions", "makeCRPSGauss.R"))
    # CRPS <- makeCRPSGauss(predMLRC$fit[, 1],
    #                       sqrt(predMLRC$v1.boot[, 1]^2 + predMLRC$v2.boot[1]^2),
    #                       X_test)
    CRPS <- abs(apply(predMLRCBoot, 1, median) - X_test) 
    MSPE <- (predMLRC$fit[, 1] - X_test)^2
    # MAE <- abs(predMLRC$fit[, 1] - X_test)
    MAE <- abs(apply(predMLRCBoot, 1, median) - X_test)
    coverage <- ( X_test >= (predMLRC$fit[, 1] - 
                               2*sqrt(predMLRC$v1.boot[, 1]^2 + predMLRC$v2.boot[1]^2))) & 
      (X_test <= (predMLRC$fit[, 1] + 
                    2 * sqrt(predMLRC$v1.boot[, 1]^2 + predMLRC$v2.boot[1]^2)))
  } else if (model_name=="WAPLS") {
    ## WAPLS reconstruction - subset to deal with all zero occurrence species
    zeros_idx <- which(colSums(y_train_prop) == 0)
    if (length(zeros_idx) > 0) {
      modWAPLS <- rioja::WAPLS(y_train_prop[, - zeros_idx], X_train)     
      predWAPLS <- predict(modWAPLS, y_test_prop, sse=TRUE, nboot=1000)
    } else {
      modWAPLS <- rioja::WAPLS(y_train_prop, X_train)     
      predWAPLS <- predict(modWAPLS, y_test_prop, sse=TRUE, nboot=1000)
    }
    source(here("functions", "makeCRPSGauss.R"))
    CRPS <- makeCRPSGauss(predWAPLS$fit[, 1], sqrt(predWAPLS$v1.boot[, 1]),
                          X_test)
    MSPE <- (predWAPLS$fit[, 1] - X_test)^2
    MAE <- abs(predWAPLS$fit[, 1] - X_test)
    coverage <- (
      X_test >=
        (predWAPLS$fit[, 1] - 2*sqrt(predWAPLS$v1.boot[, 1]))) & 
      (X_test <= 
         (predWAPLS$fit[, 1] + 2 * sqrt(predWAPLS$v1.boot[, 1])))
  } else if (model_name=="RF") {
    ## Random Forest
    train <- data.frame(moisture=X_train, y_train)
    test <- data.frame(y_test)
    rf <- randomForest(moisture ~ ., data = train)
    Rcpp::sourceCpp(here("functions", "makeCRPS.cpp"))
    CRPS <- makeCRPS(t(matrix(predict(rf, test, predict.all=TRUE)$individual, 
                              length(idx_test), 500)), X_test, 500)
    # CRPS <- abs(predict(rf, test) - X_test)
    MSPE <- (predict(rf, test) - X_test)^2
    MAE <- abs(predict(rf, test) - X_test)
    rf_CI <- t( apply( predict(rf, test, predict.all=TRUE)$individual, 1,
                       function(x) {
                         quantile(x, c(0.025,0.975))
                       }))
    coverage <- ( (X_test >= rf_CI[, 1]) & (X_test <= rf_CI[, 2]) )
  } else if (model_name=="GJAM") {
    ## GJAM model fit
    idx_hold <- (length(X_cv) - length(idx_test) + 1):length(X_cv)
    Xdf <- data.frame(x=c(X_train, X_test))
    Xdf$x[idx_hold] <- NA
    ydf <- data.frame(as.matrix(rbind(y_train, y_test)))
    colnames(ydf) <- paste("y", 1:dim(y_train)[2], sep="")
    ml <- list(ng = 5000, burnin = 500, typeNames = rep("CC", dim(y_train)[2]), 
               PREDICTX=TRUE)
    # out <- gjam(~ x, Xdf, ydf, ml)
    out <- gjam(~ x + I(x^2), Xdf, ydf, ml)
    # xMu  <- out$prediction$xpredMu        #inverse prediction of x
    # xSd  <- out$prediction$xpredSd
    xMu  <- out$prediction$xpredMu[idx_hold, 2]        #inverse prediction of x
    xSd  <- out$prediction$xpredSd[idx_hold, 2]  

    ##
    source(here("functions", "makeCRPSGauss.R"))
    
    CRPS <- makeCRPSGauss(xMu, xSd, X_test)
    MSPE <- (xMu - X_test)^2
    MAE <- abs(xMu - X_test)
    X_025 <- xMu - 2*xSd
    X_975 <- xMu + 2*xSd
    coverage <- (X_test >= X_025) & (X_test <= X_975)
  } else if (model_name=="BUMMER") {
    n_mcmc <- params$n_mcmc
    n_adapt <- params$n_adapt
    n_thin <- params$n_thin
    
    library(nimble)
    source("~/mvgp/functions/dirichlet-multinomial-nimble.R")
    ## Load the bummer code
    source("~/mvgp/functions/make-bummer.R")
    ## Load the BUMMER model
    source("~/mvgp/nimble/functions/bummer-dm-model.R")
    ## Set up the BUMMER model
    # source("~/mvgp/nimble/functions/bummer-dm-model-initialization-cv.R")
    ## Model Set-up
    ## data
    data.fit <- list(
      Y = rbind(y_train, y_test),
      X = c(X_train, rep(NA, nrow(y_test))))
    
    N <- nrow(data.fit$Y)
    N_obs <- nrow(y_train)
    
    ## constants
    constants.fit = list(
      d=d, 
      N=N, 
      N_obs=N_obs, 
      mu_X=0, 
      s_X=1, 
      count=apply(data.fit$Y, 1, sum))
    
    ## initial conditions
    inits.fit = list(
      alpha=matrix(pmin(rgamma(nrow(data.fit$Y)*d, 1, 1), 5), N, d), 
      a=rgamma(d, 1, 1),
      b=rnorm(d), 
      c=pmax(rgamma(d, 1, 1), 0.5), 
      lambda_a=rgamma(d, 1, 1), 
      s2_a=rgamma(1, 1, 1), 
      mu_b=0, 
      sigma_b=1,
      lambda_c=rgamma(d, 1, 1), 
      s2_c=rgamma(1, 1, 1))
    
    ## dimensions
    dimensions.fit = list(alpha=dim(matrix(NA, N, d)),
                          X=N, 
                          Y=dim(matrix(NA, N, d)))
    
    ## fit model
    model_fit <- nimbleModel(bummer_code, 
                             inits=inits.fit,
                             constants=constants.fit,
                             data=data.fit, 
                             dimensions=dimensions.fit)
    
    
    spec.fit <- configureMCMC(model_fit, thin=n_thin, print=TRUE)
    
    ## Remove non-reflective RW samplers for a and c
    spec.fit$removeSamplers('a')
    spec.fit$removeSamplers('c')
    ## Add in reflective RW samplers for a and c
    for (j in 1:d) {
      spec.fit$addSampler(target = paste('a[', j, ']', sep=''),
                          type = 'RW', 
                          control = list(reflective = TRUE))
      spec.fit$addSampler(target = paste('c[', j, ']', sep=''),
                          type = 'RW', 
                          control = list(reflective = TRUE))
    }
    
    ## Remove RW samplers for X
    spec.fit$removeSamplers('X')
    ## Add ESS sampler for X
    source('~/mvgp/nimble/functions/sampler-ess-univariate.R')
    for (i in (N_obs+1):N) {
      spec.fit$addSampler(target=paste('X[', i, ']', sep=''), type='ess_univariate')
    }
    
    
    ## Add Monitors
    spec.fit$addMonitors(c('a', 'b', 'c', 'X',
                           'alpha')) 
    
    ## Build MCMC
    Rmcmc.fit <- buildMCMC(spec.fit)
    
    ## Compile MCMC
    # Rmcmc.fit$run(n_mcmc)
    cm <- compileNimble(model_fit, showCompilerOutput = FALSE)
    
    Cmcmc.fit <- compileNimble(Rmcmc.fit, project = model_fit,
                               showCompilerOutput = FALSE)
    
    sink(here("model-fit", "progress", 
              "cross-validate", paste0("dm-cv-bummer-", data_source, ".txt")),
         append=TRUE)
    
    cat(paste("Starting Cross-validation", i, "\n"))
    
    samples.fit <- runMCMC(
      mcmc              = Cmcmc.fit, 
      niter             = n_mcmc + n_adapt,
      nchains           = 1, 
      samplesAsCodaMCMC = TRUE)
    
    sink()
    
    ## remove the burn-in
    samples.fit <- samples.fit[-c(1:(n_adapt / n_thin)), ]
    
    n_samples <- nrow(samples.fit)
    X_post <- matrix(0, length(X_test), n_samples)
    for (i in length(X_test):N) {
      X_post[i-N+1, ] <- samples.fit[, paste("X[", i, "]", sep="")]
    }
    
    Rcpp::sourceCpp(here::here("functions", "makeCRPS.cpp"))
    CRPS <- makeCRPS(t(X_post), X_test, n_samples)
    MSPE <- (apply(X_post, 1, mean) - X_test)^2
    MAE <- abs(apply(X_post, 1, mean) - X_test)
    X_025 <- apply(X_post, 1, quantile, prob = 0.025)
    X_975 <- apply(X_post, 1, quantile, prob = 0.975)
    coverage <- (X_test >= X_025) & (X_test <= X_975)
    rm(samples.fit)
  }
  return(list(CRPS=CRPS, MSPE=MSPE, MAE=MAE, coverage=coverage))
}