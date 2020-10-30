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
import pickle
import math

# no sideband
# [2.0430678624883347, 1.43312542155723, 1.4121729284953546, 1.398337068647945, 1.3156597431937982, 1.2399629052816066, 1.2280585389802612, 1.1285337629528844]
# [1.9839632086466101, 1.4712313655503626, 1.4479517208453963, 1.4215206802627485, 1.3587278642095806, 1.283189124722305, 1.1736956185509486, 1.1112712251097387]
lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)
sys.path.append('../package')
from package.curveplot import curveplot

binning2tagresolvedsr = [0.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 
    440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 630.0, 660.0, 690.0, 720.0, 750.0, 
    780.0, 810.0, 850.0, 910.0, 970.0, 1030.0, 1090.0, 1150.0, 1210.0, 1270.0, 1330.0, 1390.0, 1450.0, 
    1510.0, 1570.0, 1900.0, 2300.0, 2700.0, 3100.0, 9000.0]

def pickleit(obj, path):
    outfile = open(path, 'wb')
    pickle.dump(obj, outfile)
    outfile.close()

def unpickleit(path):
    infile = open(path, 'rb')
    output = pickle.load(infile)
    infile.close()
    return output

def pre_selection(df):
    #dfout = df.loc[df["region"] == 1]
    dfout = df.loc[(df["region"] == 1) | (df["region"] == 2)]
    dfout = dfout.loc[dfout["regime"] == 1]
    dfout = dfout.loc[dfout["nTags"] == 2]
    #dfout.drop(columns=["nTags", "MCChannelNumber", "MV2c10B3", "pTJ3", "phiJ3", "etaJ3", "region", "regime", "phiB1", "phiB2"], inplace=True)
    dfout.drop(columns=["nTags", "MCChannelNumber", "region", "regime"], inplace=True)
    return dfout

def significance_binned(backgrounds, signal, logsig=True, portion=1):
    backgrounds_content = np.array(backgrounds)/portion
    signal_content = np.array(signal)/portion
    total = 0

    if not logsig:
        return sum(signal_content)/math.sqrt(sum(backgrounds_content))

    for each_b, each_s in zip(backgrounds_content, signal_content):
        if each_b > 0 and each_s > 0:
            total += 2 * ((each_s + each_b) * math.log(1 + each_s/each_b) - each_s)
    return math.sqrt(total)

def cal_significance(sig, bkg, signal_mass, bins):
    test_bkg = bkg["mVHres"].to_numpy()/1000.
    bkgresultweight = bkg["weight"]
    significance = []
    for each_mass in signal_mass:
        thissignal = sig.loc[sig['mass'] == each_mass]
        sigresultweight = thissignal["weight"]
        sigresult = thissignal["mVHres"].to_numpy()/1000.

        sighist, bin_edges = np.histogram(sigresult.flatten(), bins=bins, weights=sigresultweight.to_numpy())
        bkghist, bin_edges = np.histogram(test_bkg, bins=bins, weights=bkgresultweight.to_numpy())
        significance.append((each_mass, significance_binned(bkghist, sighist, logsig=True)))
    return dict(significance)

def main():
    # setting
    loadcsv = False
    dotraining = True
    domassplot = True

    signal_mass = [300, 400, 420, 440, 460, 500, 600, 700]#, 800, 900, 1000, 1200, 1400, 1600, 2000]
    if loadcsv:
        signal_all = pd.DataFrame()
        #signal_mass_train = [300, 400, 420, 440, 460, 600, 700]
        for each in signal_mass:
            df_temp = pd.read_csv(str(each) + ".csv", index_col=0)
            df_temp = pre_selection(df_temp)
            df_temp["mass"] = each
            signal_all = pd.concat([df_temp, signal_all], ignore_index=True)
        background_all = pd.read_csv("background.csv", index_col=0)
        background_all = pre_selection(background_all)
        #background_all["mass"] = np.random.rand(len(background_all)) * 2000
        background_all["mass"] = np.random.choice(signal_mass, len(background_all))
        pickleit(background_all, "bkg.pickle")
        pickleit(signal_all, "sig.pickle")
        exit(1)
    else:
        background_all = unpickleit("bkg.pickle")
        signal_all = unpickleit("sig.pickle")
        #result_2tag = cal_significance(signal_all, background_all, signal_mass, binning2tagresolvedsr)
        result_2tag = dict([(300, 5.636571301755268), (400, 15.682911909977934), (420, 17.402946494678936), (440, 19.59381319990806), (460, 22.61216334233906), (500, 28.083905499575355), (600, 40.36491380816585), (700, 50.96207372482288), (800, 60.86608130389973), (900, 66.48638123492984), (1000, 72.15575742621934), (1200, 77.40686627138746), (1400, 70.99760718165942), (1600, 60.77058118028304), (2000, 48.35965326790206)])
        # background_all.drop(columns=["mVHres"], inplace=True)
        # signal_all.drop(columns=["mVHres"], inplace=True)

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
    if dotraining:
        start = time.time()
        with tf.device('/CPU:' + str(0)):
            model = keras.models.Sequential()
            model.add(keras.layers.Dense(100, input_dim=len(train_x_before.columns), activation="selu", kernel_initializer="lecun_normal"))
            #model.add(keras.layers.LeakyReLU(alpha=0.2))
            model.add(keras.layers.Dense(100, activation="selu", kernel_initializer="lecun_normal"))
            #model.add(keras.layers.LeakyReLU(alpha=0.2))
            model.add(keras.layers.Dense(100, activation="selu", kernel_initializer="lecun_normal"))
            #model.add(keras.layers.LeakyReLU(alpha=0.2))
            model.add(keras.layers.Dense(100, activation="selu", kernel_initializer="lecun_normal"))
            model.add(keras.layers.Dense(1, activation=activations.sigmoid))
            opt1 = keras.optimizers.SGD(lr=0.0002)
            opt2 = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999)
            model.compile(loss="binary_crossentropy", optimizer=opt2, weighted_metrics=["accuracy"])#0.00000000005
            history = model.fit(train_x, np.array(train_y), sample_weight=train_weight, epochs=5, validation_data=(test_x, np.array(test_y), test_weight), shuffle=True, batch_size=100)
            #history = model.fit(train_x, np.array(train_y), epochs=13, validation_data=(test_x, np.array(test_y)), shuffle=True, batch_size=70)
    
        model.save('testmodel') 
        end = time.time() - start
        print(end)
        # pd.DataFrame(history.history).plot(figsize=(8, 5))
        # plt.grid(True)
        # plt.gca().set_ylim(0, 1)
        # plt.show()
    else:
        model = tf.keras.models.load_model('testmodel')

    if domassplot:
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
    



    significance = []
    for each_mass in signal_mass:
        test_bkg = test_bkg.assign(mass=each_mass)
        thissignal = test_signal.loc[test_signal['mass'] == each_mass]
        bkgresult = model.predict(scaler.transform(test_bkg.drop(columns=["weight"]).to_numpy()))
        sigresult = model.predict(scaler.transform(thissignal.drop(columns=["weight"]).to_numpy()))
        bkgresultweight = test_bkg["weight"]
        sigresultweight = thissignal["weight"]
        bins = np.linspace(0, 1, 50)
        sighist, bin_edges = np.histogram(sigresult.flatten(), bins=bins, weights=sigresultweight.to_numpy())
        bkghist, bin_edges = np.histogram(bkgresult.flatten(), bins=bins, weights=bkgresultweight.to_numpy())
        print(each_mass, significance_binned(bkghist, sighist, logsig=True, portion=0.4), result_2tag[each_mass])
        significance.append(significance_binned(bkghist, sighist, logsig=True, portion=0.4)/result_2tag[each_mass])
        # curveplot([(bin_edges[0:-1] + bin_edges[1:])/2]*2, [sighist/0.4, bkghist/0.4], filename="rawoutput" + str(each_mass), ylimit=[0,14000], labels=["sig", "bkg"], xlimit=[0, 1], 
        #            yshift=0.05, xshift=0.03, ylabel="arbitary unit", xlabel="NN output", log_y=True)
    print(significance)


if __name__ == "__main__":
    main()
