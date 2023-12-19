library(doParallel)
library(Matrix)
library(glmnet)
library(stringr)
library(dplyr)
library(tictoc)


PeakFuse <- function(X, window = 100000){
  s <- colnames(X)
  
  new.s=str_split_fixed(s,":",2)%>%data.frame()
  
  new.s[,c(3,4)] <- str_split_fixed(new.s$X2,"-",2)
  
  new.s[,5] <- colnames(X)
  
  colnames(new.s) <- c("chr","inerval","start","end","fullname")
  
  X.smooth <- Matrix(rep(0,nrow(X)),nrow=nrow(X),ncol = 1, sparse = T)
  
  colnames(X.smooth)="1"
  
  #Go through all chromosomes
  for (chro in unique(new.s$chr) ) {
    # Get peak position info in current chromosome
    chrdata=subset(new.s, new.s$chr==chro )
    chrdata$start=as.numeric(chrdata$start)
    chrdata$end=as.numeric(chrdata$end)
    
    print(chro)
    
    #If current chromosome only has 1 peak, add it to our data. 
    if (nrow(chrdata)<=1){
      new.peak.name <- chrdata$fullname[1]
      X.smooth <- cbind(X.smooth, X[ , chrdata$fullname[1]])
      colnames(X.smooth)[length(colnames(X.smooth))]=new.peak.name
      next
    }
    
    #Start from first peak, find its neighbors and fuse them. The neighbors will be
    #recorded so that one can avoid duplicate fusion.
    i=1
    processed.record <- c()
    chrdata <- chrdata[order(chrdata$start),]
    while (i <=nrow(chrdata) ){
      if (chrdata$fullname[i]%in%processed.record)
      {i=i+1
      next
      }
      else{
        neighbors.ind <- which( (chrdata$start>chrdata$end[i])&(chrdata$start-chrdata$end[i]<window) ) 
        
        neighbors <- chrdata$fullname[neighbors.ind]
        neighbors <- c(neighbors,chrdata$fullname[i])
        if (length(neighbors)>1){
          X.smooth <- cbind(X.smooth, rowSums(X[ ,neighbors]) )
          processed.record <- c(processed.record,neighbors)
          # new.peak.name <- paste(chro,":",chrdata$start[neighbors.ind[1]],
          #                        "-",chrdata$end[neighbors.ind[length(neighbors.ind)]]
          #                        ,"new")'
          new.peak.name <- "new"
          colnames(X.smooth)[length(colnames(X.smooth))]=new.peak.name
          i=i+1
        }
        else{
          new.peak.name <- chrdata$fullname[i]
          X.smooth <- cbind(X.smooth, X[ , chrdata$fullname[i]])
          colnames(X.smooth)[length(colnames(X.smooth))]=new.peak.name
          i=i+1
        }
      }
      print(i)
    }
    
  }
  
  dim(X.smooth)
  X.smooth.all <- X.smooth[,-1]
  return(X.smooth.all)
}

myic.glmnet <- function (x, y, crit = c("bic", "aic", "aicc", "hqc"), alpha = 1, 
                         ...) 
{
  crit = match.arg(crit)
  n = length(y)
  model = glmnet(x = x, y = y, alpha = alpha, ...)
  coef = coef(model)
  lambda = model$lambda
  df = model$df
  if (alpha == 0) {
    xs = scale(x)
    I = diag(ncol(x))
    xx = t(xs) %*% xs
    for (i in 1:length(lambda)) {
      aux = solve(xx + I * lambda[i])
      df[i] = sum(diag(xs %*% aux %*% t(xs)))
    }
  }
  yhat = cbind(1, x) %*% coef
  residuals = (y - yhat)
  mse = colMeans(residuals^2)
  sse = colSums(residuals^2)
  nvar = df + 1
  bic = n * log(mse) + nvar * log(n)
  aic = n * log(mse) + 2 * nvar
  aicc = aic + (2 * nvar * (nvar + 1))/(n - nvar - 1)
  hqc = n * log(mse) + 2 * nvar * log(log(n))
  sst = (n - 1) * var(y)
  r2 = 1 - (sse/sst)
  adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1))
  crit = switch(crit, bic = bic, aic = aic, aicc = aicc, hqc = hqc)
  selected = best.model = which(crit == min(crit))
  ic = c(bic = bic[selected], aic = aic[selected], aicc = aicc[selected], 
         hqc = hqc[selected])
  result = list(coefficients = coef[, selected], ic = ic, 
                lambda = lambda[selected], nvar = nvar[selected], glmnet = model, 
                residuals = residuals[, selected], fitted.values = yhat[, 
                                                                        selected], ic.range = crit, df = df, call = match.call())
  class(result) = "ic.glmnet"
  return(result)
}

GetPredictionLASSO <- function(X, Y, train.index, f=myic.glmnet){
  X.train <- X[train.index, ]
  Y.train <- Y[train.index, ]
  X.test <- X[-train.index, ]
  Y.test <- Y[-train.index, ]
  
  func.Optimal.lambda <-  function(ii){
    icmodel=f(x=X.train,y=Y.train[,ii],crit = "bic")
    para.num = icmodel[["nvar"]]-1
    lambda = icmodel[["lambda"]]
    if (para.num<5) {
      icmodel=icmodel[["glmnet"]]
      ind = min(which(icmodel[["df"]]>=5))# At least 5 parameters
      lambda=icmodel[["lambda"]][ind]
      para.num=icmodel[["df"]][ind]
    }else{icmodel=icmodel[["glmnet"]]}
    
    y.pred <- predict(icmodel, X.test, s=lambda)
    return(c(y.pred))
  }
  
  cl.cores=8
  cl = makeCluster(cl.cores)
  registerDoParallel(cl)
  tic()
  system.time(result.pre <- foreach(ii=1:ncol(Y.test),
                                    .combine = "rbind",
                                    .packages = c("glmnet")
  ) %dopar% func.Optimal.lambda(ii))
  stopCluster(cl)
  toc()
  result.pre=as.data.frame(t(result.pre))
  colnames(result.pre)=colnames(Y.test)
  row.names(result.pre)=row.names(Y.test)
  result.pre=as.data.frame(lapply(result.pre, as.numeric))
  return(result.pre)
}

#When one predict gene with corresponding peaks, this function should be used.
GetPredictionwithCrspPeaksLASSO <- function(X, Y, train.index, f=myic.glmnet, xuyaoduiyingbiao){
  X.train <- X[train.index, ]
  Y.train <- Y[train.index, ]
  X.test <- X[-train.index, ]
  Y.test <- Y[-train.index, ]
  
  
  func.Optimal.lambda.match <-  function(ii){
    peak.ind <- xuyaoduiyingbiao[[colnames(Y.train)[ii]]]
    
    if(length(peak.ind)<=1){peak.ind <- c(peak.ind, peak.ind+1)}
    
    icmodel=f(x=X.train[, peak.ind], y=Y.train[,ii], crit = "bic")
    para.num = icmodel[["nvar"]]-1
    
    lambda = icmodel[["lambda"]]
    if (para.num<2) {
      icmodel=icmodel[["glmnet"]]
      ind = min(which(icmodel[["df"]]>=2))
      if(is.infinite(ind)){ind = min(which(icmodel[["df"]]>=1))}
      lambda=icmodel[["lambda"]][ind]
      para.num=icmodel[["df"]][ind]
    }else{icmodel=icmodel[["glmnet"]]}
    
    y.pred <- predict(icmodel, X.test[, peak.ind], s=lambda)
    
    return(c(y.pred))
  }
  
  cl.cores=8
  cl = makeCluster(cl.cores)
  registerDoParallel(cl)
  tic()
  result.pre <- foreach(ii=1:ncol(Y.test),
                        .combine = "rbind",
                        .packages = c("glmnet")
  ) %dopar% func.Optimal.lambda.match(ii)
  stopCluster(cl)
  toc()
  result.pre=as.data.frame(t(result.pre))
  colnames(result.pre)=colnames(Y.test)
  row.names(result.pre)=row.names(Y.test)
  result.pre=as.data.frame(lapply(result.pre, as.numeric))
  return(result.pre)
}



#KNN Method
GetPredictionKNN <- function(k = 10, X, Y, train){
  Xlogi.train <- X[train, ] != 0
  Xlogi.test <- X[-train, ] != 0
  
  Y.train <- Y[train, ]
  Y.test <- Y[-train, ]
  
  KnnPre <- function(cell){
    co <- Xlogi.train + Xlogi.test[rep(cell, nrow(Xlogi.train)), ]
    
    a=rowSums(co == 2)
    c.and.d <- rowSums(co == 1)
    simi = a/(a + c.and.d )
    neighbor.index = order(simi,decreasing = T)[1:k]
    
    resForThisTestCell <- colMeans(Y.train[neighbor.index, ])
    return(resForThisTestCell)
  }
  
  
  cl.cores = 8
  cl = makeCluster(6)
  registerDoParallel(cl)
  
  tic()
  
  KNNRes <- foreach(cell = 1:nrow(Y.test),
                    .combine = "rbind",
                    .packages = "Matrix"
  ) %dopar% KnnPre(cell)
  
  stopCluster(cl)
  
  toc()
  
  KNNRes=as.data.frame(KNNRes)
  row.names(KNNRes)=row.names(Y.test)
  
  return(KNNRes) 
}

Metrics <- function(result.pre,Y.test){
  func.corr <- function(ii){
    return( cor(Y.test[,ii],result.pre[,ii]) )
  }
  
  func.spearman <- function(ii){
    return( cor(Y.test[,ii],result.pre[,ii],method="spearman") )
  }
  func.r2 <- function(ii){
    y.bar = mean(Y.test[, ii])
    sst = sum( (Y.test[, ii] - y.bar)^2 )
    ssr = sum( (result.pre[, ii] - y.bar)^2 )
    Rsquare = ssr / sst
  }
  
  cl.cores=8
  cl = makeCluster(cl.cores)
  registerDoParallel(cl)
  Gene.Corr = foreach(ii=1:ncol(result.pre),
                      .combine = "c",
                      .packages = "Matrix"
  )%dopar%func.corr(ii)
  Gene.spearman = foreach(ii=1:ncol(result.pre),
                          .combine = "c",
                          .packages = "Matrix"
  )%dopar%func.spearman(ii)
  Gene.r2 = foreach(ii=1:ncol(result.pre),
                    .combine = "c",
                    .packages = "Matrix"
  )%dopar%func.r2(ii)
  
  stopCluster(cl)
  return(data.frame(Gene=colnames(result.pre),
                    pearson=Gene.Corr,
                    spearman=Gene.spearman,
                    R2 = Gene.r2)
  )
}


OneSplit <- function(seed, prop, X){
  set.seed(seed)
  train <- sample(1:nrow(X), ceiling(prop * nrow(X)) )
  return(train)
}

RepeatSplit <- function(reptimes, X, Y, prop, 
                        metric.chosen = c("pearson", "spearman", "R2"), method = "Lasso",duiyingbiao = NULL, k=10){
  MetricinRepeatation <- data.frame(Gene = colnames(Y))
  for (i in 1:reptimes) {
    print(i)
    train.index = OneSplit(seed = i, X = X, prop = prop)
    if(method == "Lasso"){PredRes = GetPredictionLASSO(X, Y, train.index, f = myic.glmnet)}
    else if(method == "Match"){PredRes = GetPredictionwithCrspPeaksLASSO(X, Y, train.index, 
                                                                         f=myic.glmnet, xuyaoduiyingbiao = duiyingbiao)}
    else if(method == "KNN"){PredRes = GetPredictionKNN(k, X, Y, train.index)}
    
    
    
    MetricOneTime = Metrics(PredRes, Y[-train.index, ])
    MetricinRepeatation <- cbind(MetricinRepeatation, MetricOneTime[, metric.chosen])
  }
  return(MetricinRepeatation)
  
}



#Gene Score
geneScorePara <- function(pos, crsp.list, X, Y){
  #获得所有tile
  tiles_name <- lapply(crsp.list, function(x){
    gsub("chr", "", colnames(X)[x])
  })
  
  
  #所有tile的start
  tiles_start <- lapply(tiles_name, function(x){
    new.s <- str_split_fixed(x, ":", 2)%>%as.data.frame()
    new.s[,c(3,4)] <- str_split_fixed(new.s$V2,"-",2)
    peaks <- new.s[, -2]
    return(as.numeric(as.vector(peaks[, 2])))
  })
  
  
  #计算gene score
  func <- function(i){
    gene <- pos[, 1][i]
    
    
    dis_weight <- exp(-abs(tiles_start[[gene]] - pos[i, "start_position"])/5000) + exp(-1)
    
    dis_weight_scale <- dis_weight * 1/(pos[i, "end_position"] - pos[i, "start_position"])
    
    if(length(dis_weight_scale) == 1){
      dis_weight_scale = 5
    }
    else if(length(table(dis_weight_scale)) == 1){
      dis_weight_scale = seq(1, 5, length.out = length(dis_weight_scale))
    }
    else{
      dis_weight_scale <- 1 + 4 * (dis_weight_scale - min(dis_weight_scale))/(max(dis_weight_scale) - min(dis_weight_scale))
    }
    
    
    if(is.null(dim(X[, crsp.list[[gene]]]))){
      s <- X[, crsp.list[[gene]]] * dis_weight_scale
    }
    else{
      s <- X[, crsp.list[[gene]]] %*% dis_weight_scale
    }
    
    
    return(s)
  }
  
  tic()
  cl <- makeCluster(8)
  registerDoParallel(cl)
  score <- foreach(i = 1:nrow(pos),
                   .packages = "Matrix",
                   .combine = "cbind"
  )%dopar%func(i)
  stopCluster(cl)
  
  toc()
  
  #depth normalization
  M = 10000
  
  colnames(score) <- colnames(Y)
  rownames(score) <- rownames(Y)
  
  score_matrix <- CreateSeuratObject(t(score))
  
  score_matrix <- NormalizeData(score_matrix)
  
  score_matrix <- score_matrix@assays[["RNA"]]@data
  
  return(score_matrix)
}



