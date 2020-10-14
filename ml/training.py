import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import activations
import matplotlib.pyplot as plt
import time

signal_all = pd.DataFrame()
signal_mass = [300, 400, 420, 440, 460, 500, 600, 700]#, 800, 900, 1000, 1200, 1400, 1600, 2000]
for each in signal_mass:
    df_temp = pd.read_csv(str(each) + ".csv", index_col=0)
    df_temp["mass"] = each
    df_temp.drop(columns=["nTags", "MCChannelNumber"], inplace=True)
    signal_all = pd.concat([df_temp, signal_all], ignore_index=True)
background_all = pd.read_csv("background.csv", index_col=0)
#background_all["mass"] = np.random.rand(len(background_all)) * 2000
background_all["mass"] = np.random.choice(signal_mass, len(background_all))
background_all.drop(columns=["nTags", "MCChannelNumber"], inplace=True)

train_bkg, test_bkg = train_test_split(background_all, test_size=0.4, random_state=2)
train_signal, test_signal = train_test_split(signal_all, test_size=0.4, random_state=2)
train_y = len(train_bkg) * [1] + len(train_signal) * [0]
test_y = len(test_bkg) * [1] + len(test_signal) * [0]
train_x = pd.concat([train_bkg, train_signal], ignore_index=True)
test_x = pd.concat([test_bkg, test_signal], ignore_index=True)
train_weight = train_x["weight"].to_numpy()
train_x.drop(columns=["weight"], inplace=True)
test_weight = test_x["weight"].to_numpy()
test_x.drop(columns=["weight"], inplace=True)

scaler = StandardScaler()
scaler.fit_transform(train_x)
scaler.transform(test_x)

start = time.time()
with tf.device('/CPU:' + str(0)):
    model = keras.models.Sequential()
    #model.add(keras.layers.InputLayer(input_shape=input_shape))
    model.add(keras.layers.Dense(100, input_dim=len(train_x.columns), activation="relu"))
    model.add(keras.layers.Dense(100, activation="relu"))
    # model.add(keras.layers.Dense(50, activation="relu"))
    model.add(keras.layers.Dense(1, activation=activations.sigmoid))
    model.compile(loss="binary_crossentropy", optimizer=keras.optimizers.SGD(lr=0.00000005), metrics=["accuracy"])
    # model.add(keras.layers.Dense(2, activation="softmax"))
    # model.compile(loss="categorical_crossentropy", optimizer=keras.optimizers.SGD(lr=0.00000005), metrics=["accuracy"])
    #history = model.fit(train_x.to_numpy(), np.array(train_y), sample_weight=train_weight, epochs=10, validation_data=(test_x.to_numpy(), np.array(test_y), test_weight), shuffle=True)
    history = model.fit(train_x.to_numpy(), np.array(train_y), epochs=10, validation_data=(test_x.to_numpy(), np.array(test_y)), shuffle=True, batch_size=70)
model.save('testmodel') 
end = time.time() - start
print(end)
pd.DataFrame(history.history).plot(figsize=(8, 5))
plt.grid(True)
plt.gca().set_ylim(0, 1)
plt.show()


# new_model = tf.keras.models.load_model('testmodel')
# new_model.predict()
