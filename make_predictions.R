# Fusing peaks if needed
# It takes about 40 minutes for window = 100000 
Xpf <- PeakFuse(X, window = 100000)
dim(Xpf)


# Repeating train-test split
LassoRepeat <- RepeatSplitFixedTrain(reptimes = 5,
                                     X = X,
                                     Y = Y,
                                     train_index = train_index,
                                     method = "Lasso")
LassoMatchRepeat <- RepeatSplitFixedTrain(reptimes = 5,
                                          X = X,
                                          Y = Y,
                                          train_index = train_index,
                                          method = "Match",
                                          duiyingbiao = crsp.list)
LassoPFRepeat <- RepeatSplitFixedTrain(reptimes = 5,
                                       X = Xpf,
                                       Y = Y,
                                       train_index = train_index,
                                       method = "Lasso")
KNNRepeat <- RepeatSplitFixedTrain(reptimes = 5,
                                   X = X,
                                   Y = Y,
                                   train_index = train_index,
                                   method = "KNN")

# Gene score
score_matrix = geneScorePara(pos, crsp.list, X, Y)
score_metrics <- Metrics(t(score_matrix), Y)