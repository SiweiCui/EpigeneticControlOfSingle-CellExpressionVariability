import pyreadr
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import scipy.linalg as la
from keras.models import Sequential
from keras.layers import Dense, Dropout
import matplotlib.pyplot as plt
from sklearn.preprocessing import MaxAbsScaler
from keras import backend as K


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
X_train = MaxAbsScaler().fit_transform(X_train)
X_test = MaxAbsScaler().fit_transform(X_test)

def batch_generator(x, y, batch_size):
    number_of_batches = x.shape[0]//batch_size
    counter = 0
    shuffle_index = np.arange(x.shape[0])
    np.random.shuffle(shuffle_index)
    x = x[shuffle_index, :]
    y = y[shuffle_index, :]
    while 1:
        index_batch = shuffle_index[batch_size*counter: batch_size*(counter+1)]
        x_batch = x[index_batch, :].todense()
        y_batch = y[index_batch, :].todense()
        counter+=1
        yield(x_batch, y_batch)
        if counter >= number_of_batches:
            np.random.shuffle(shuffle_index)
            counter = 0
            
NNPred = pd.DataFrame()
for gene in Y_colnames["Y_smo_colnames"]:

    K.clear_session()
    
    print(np.where(Y_colnames["Y_smo_colnames"] == gene))
    print(gene)
    
    
    model = Sequential()

    model.add(Dropout(0.6,  input_shape = (X_shape2, )))
    model.add(Dense(units=128, activation="relu"))
    model.add(Dropout(0.2))
    model.add(Dense(units=64, activation="relu"))
    model.add(Dropout(0.2))
    model.add(Dense(units=1, activation="relu"))


    model.compile(loss="mean_squared_error",
                  optimizer="sgd",
                  metrics="mean_squared_error")
    
    batch_size = 100
    epo = 25
    Y_train_genewise = Y_train[:, Y.col == gene].copy()
    train_steps =(X_train.shape[0]//batch_size)
    his=model.fit(batch_generator(X_train.asfptype(), Y_train_genewise, batch_size),epochs=epo,
                   steps_per_epoch=train_steps, verbose=0)
    
    pre = model.predict(X_test)
    NNPred[gene] = pre.reshape(-1)

    