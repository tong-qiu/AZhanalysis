import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import activations
import matplotlib.pyplot as plt
import time
import os
import sys

lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)
sys.path.append('../package')
from package.curveplot import curveplot

signal_all = pd.DataFrame()
signal_mass = [300, 400, 420, 440, 460, 500, 600, 700]#, 800, 900, 1000, 1200, 1400, 1600, 2000]
#signal_mass_train = [300, 400, 420, 440, 460, 600, 700]
for each in signal_mass:
    df_temp = pd.read_csv(str(each) + ".csv", index_col=0)
    df_temp["mass"] = each
    # df_temp["mLL"] *= 1000
    df_temp.drop(columns=["nTags", "MCChannelNumber", "mVHres"], inplace=True)
    signal_all = pd.concat([df_temp, signal_all], ignore_index=True)
background_all = pd.read_csv("background.csv", index_col=0)
#background_all["mass"] = np.random.rand(len(background_all)) * 2000
background_all["mass"] = np.random.choice(signal_mass, len(background_all))
background_all.drop(columns=["nTags", "MCChannelNumber", "mVHres"], inplace=True)

train_bkg, test_bkg = train_test_split(background_all, test_size=0.4, random_state=2)
train_signal, test_signal = train_test_split(signal_all, test_size=0.4, random_state=2)
train_y = len(train_bkg) * [0] + len(train_signal) * [1]
test_y = len(test_bkg) * [0] + len(test_signal) * [1]
train_x = pd.concat([train_bkg, train_signal], ignore_index=True)
test_x = pd.concat([test_bkg, test_signal], ignore_index=True)
train_weight = train_x["weight"].to_numpy()
train_x.drop(columns=["weight"], inplace=True)
test_weight = test_x["weight"].to_numpy()
test_x.drop(columns=["weight"], inplace=True)

scaler = StandardScaler()
train_x_before = train_x
train_x = scaler.fit_transform(train_x)
test_x = scaler.transform(test_x)

# train_y = np.array(train_y)[[each for each in train_x[train_x.mass != 500].index.values]]
# train_x = train_x.drop(train_x[train_x.mass == 500].index)

# start = time.time()
# with tf.device('/CPU:' + str(0)):
#     model = keras.models.Sequential()
#     model.add(keras.layers.Dense(100, input_dim=len(train_x_before.columns), activation="relu"))
#     model.add(keras.layers.Dense(100, activation="relu"))
#     model.add(keras.layers.Dense(100, activation="relu"))
#     model.add(keras.layers.Dense(1, activation=activations.sigmoid))
#     model.compile(loss="binary_crossentropy", optimizer=keras.optimizers.SGD(lr=0.005), metrics=["accuracy"])#0.00000000005
#     history = model.fit(train_x, np.array(train_y), sample_weight=train_weight, epochs=5, validation_data=(test_x, np.array(test_y), test_weight), shuffle=True, batch_size=70)
#     #history = model.fit(train_x, np.array(train_y), epochs=13, validation_data=(test_x, np.array(test_y)), shuffle=True, batch_size=70)
# model.save('testmodel') 
# end = time.time() - start
# print(end)
# pd.DataFrame(history.history).plot(figsize=(8, 5))
# plt.grid(True)
# plt.gca().set_ylim(0, 1)
# plt.show()


model = tf.keras.models.load_model('testmodel')

for each_mass in signal_mass:
    #test_bkg["mass"] = each_mass
    #test_bkg.loc[:, 'mass'] = each_mass
    test_bkg = test_bkg.assign(mass=each_mass)
    thissignal = test_signal.loc[test_signal['mass'] == each_mass]
    bkgresult = model.predict(scaler.transform(test_bkg.drop(columns=["weight"]).to_numpy()))
    sigresult = model.predict(scaler.transform(thissignal.drop(columns=["weight"]).to_numpy()))
    bkgresultweight = test_bkg["weight"]
    sigresultweight = thissignal["weight"]
    # print(len(sigresult), len(sigresultweight))
    # print(sigresult)

    
    bins = np.linspace(0, 1, 50)
    sighist, bin_edges = np.histogram(sigresult.flatten(), density=True, bins=bins, weights=sigresultweight.to_numpy())
    bkghist, bin_edges = np.histogram(bkgresult.flatten(), density=True, bins=bins, weights=bkgresultweight.to_numpy())
    curveplot([(bin_edges[0:-1] + bin_edges[1:])/2]*2, [sighist, bkghist], filename="output" + str(each_mass), ylimit=[0,20], labels=["sig", "bkg"], xlimit=[0, 1], 
              yshift=0.05, xshift=0.03, ylabel="arbitary unit", xlabel="NN output")