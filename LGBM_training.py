from lightgbm import LGBMRegressor
import pyreadr
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import scipy.linalg as la
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler


Data = pyreadr.read_r("Donor31800_Signac_Python.RData")


X_summary = Data["X_smo_sumy"]
Y_summary = Data["Y_smo_sumy"]
X_rownames = Data["X_smo_rownames"]
X_colnames = Data["X_smo_colnames"]
Y_rownames = Data["Y_smo_rownames"]
Y_colnames = Data["Y_smo_colnames"]
X_shape1 = Data["X_shape1"].iloc[0, 0]
X_shape2 = Data["X_shape2"].iloc[0, 0]
Y_shape1 = Data["Y_shape1"].iloc[0, 0]
Y_shape2 = Data["Y_shape2"].iloc[0, 0]
train_index = (Data["train_index"].iloc[1, :] - 1).tolist()


X = sp.coo_matrix( (X_summary["x"], (X_summary["i"]-1 , X_summary["j"]-1)), shape=(X_shape1, X_shape2) ).tocsr()
Y = sp.coo_matrix( (Y_summary["x"], (Y_summary["i"]-1, Y_summary["j"]-1)), shape=(Y_shape1, Y_shape2) ).tocsr()
X.row = X_rownames["X_smo_rownames"]
X.col = X_colnames["X_smo_colnames"]
Y.row = Y_rownames["Y_smo_rownames"]
Y.col = Y_colnames["Y_smo_colnames"]


X_train = X[train_index, :]
X_test = X[[i for i in range(X_shape1) if i not in train_index], :]
Y_train = Y[train_index, :]
Y_test = Y[[i for i in range(X_shape1) if i not in train_index], :]



YtrainF = pd.DataFrame(Y_train.toarray(), 
                       index=train_index,
                        columns=Y_colnames["Y_smo_colnames"])
YtestF = pd.DataFrame(Y_test.toarray(),
                      index= [i for i in range(X_shape1) if i not in train_index],
                        columns=Y_colnames["Y_smo_colnames"])


LGBMSmoothPredict = pd.DataFrame()

for gene in Y_colnames["Y_smo_colnames"]:
    print(np.where(Y_colnames["Y_smo_colnames"] == gene))
    
    
    lgbm_iter = LGBMRegressor(objective="regression", n_estimators=30)
    lgbm_iter.fit( X_train.asfptype(), np.array(YtrainF[[gene]]).reshape(-1) )
    pred = lgbm_iter.predict(X_test.asfptype())
    LGBMSmoothPredict[gene] = pred

  