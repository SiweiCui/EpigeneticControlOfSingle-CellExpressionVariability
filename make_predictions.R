# Fusing peaks if needed
# It takes about 40 minutes for window = 100000 
Xpf <- PeakFuse(X, window = 100000)
dim(Xpf)


# Repeating train-test split
LassoRepeat <- RepeatSplit(reptimes = 5, X = X, Y = Y, prop = 0.7, method = "Lasso")
LassoMatchRepeat <- RepeatSplit(reptimes = 5, X = X, Y = Y, prop = 0.7, method = "Match", duiyingbiao = crsp.list)
LassoPFRepeat <- RepeatSplit(reptimes = 5, X = Xpf, Y = Y, prop = 0.7, method = "Lasso")
KNNRepeat <- RepeatSplit(reptimes = 5, X = X, Y = Y, prop = 0.7, method = "KNN")


# Gene score
score_matrix = geneScorePara(pos, crsp.list, X, Y)
score_metrics <- Metrics(t(score_matrix), Y)