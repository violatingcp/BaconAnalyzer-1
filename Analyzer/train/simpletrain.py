#!/usr/local/bin/python2.7

from sys import exit 
from os import environ
environ['KERAS_BACKEND'] = 'tensorflow'

import numpy as np

import setGPU
from  convert import process_file
from keras.layers import Input, Dense, Dropout, Activation, Conv2D, MaxPooling2D, LSTM, Convolution1D, MaxPooling1D, MaxPooling1D
from keras.models import Model,Sequential 
from keras.callbacks import ModelCheckpoint 
from keras.optimizers import Adam 
from keras.utils import np_utils
from keras import backend as K
from keras.preprocessing import sequence
K.set_image_data_format('channels_last')


def reform(arr_x,arr_y, train_frac, val_frac, label):
    n = len(arr_x)
    ns = {}
    ns['train'] = (0,int(train_frac*n))
    ns['val']   = (ns['train'][1],ns['train'][1]+int(val_frac*n))
    ns['test']  = (ns['val'][1],n)
    weight_norm = 100. / n
    x = {}; y = {}; w = {}; extras = {}
    for subset in ['train','val','test']:
        n_ = ns[subset]
        x[subset] = arr_x[n_[0]:n_[1]]
        y[subset] = arr_y[n_[0]:n_[1]]
        w[subset] = weight_norm * np.ones(n_[1]-n_[0])
    return {'x':x,'y':y,'w':w}
    
def load_data(train_frac, val_frac):
    arr_x = process_file("/eos/cms/store/group/phys_exotica/dijet/qbert/qbertbits-v13-gen/QCD_Pt_1000to1400_13TeV_pythia8_ext.root",30,["part_pt","part_eta","part_phi"])
    arr_y = process_file("/eos/cms/store/group/phys_exotica/dijet/qbert/qbertbits-v13-gen/QCD_Pt_1000to1400_13TeV_pythia8_ext.root",1 ,["jet_pt"])
    arr1_y = []
    for pArr in arr_y:
        arr1_y.append(pArr[0][0])
    arr1_y = np.array(arr1_y)
    data = reform(arr_x,arr1_y,train_frac,val_frac,0)
    x = {}; y = {}; w = {}
    for subset in ['train','val','test']:
        w[subset] = data['w'][subset]
        x[subset] = data['x'][subset]
        y[subset] = data['y'][subset]
    return x,y,w

x,y,w = load_data(0.5, 0.25)
dim0 = x['train'].shape[1]
dim1 = x['train'].shape[2]
print "Test:",dim0,dim1
inputs = Input(shape=(dim0, dim1), name='input')
conv = Convolution1D(32, 32, padding='valid', activation='relu', input_shape=(dim0,dim1))(inputs)
conv = Convolution1D(16, 1, padding='valid', activation='relu')(conv)
conv = Convolution1D(8, 1, padding='valid', activation='relu')(conv)
conv = Convolution1D(4, 1, padding='valid', activation='relu')(conv)
#l = Dense(32, activation='relu')(conv)
#l = Dense(64, activation='relu')(l)
#l = Dense(64, activation='relu')(l)
lstm = LSTM(100)(conv)
output = Dense(1)(l)
model = Model(inputs=inputs, outputs=output)

model.compile(optimizer=Adam(),loss='mean_squared_error',metrics=['accuracy'])

#model = Sequential()
#model.add(LSTM(1,input_shape=(2,10)))
#model.add(Dense(1))
#model.compile(loss='mean_squared_error', optimizer='adam')

model.fit(x['train'], y['train'],batch_size=50, epochs=30, verbose=2)
print model.summary()

y_pred = model.predict(x['test'])
print y_pred
print y['test']
print y['train']

score = model.evaluate(x['test'], y['test'], batch_size=32, verbose=1, sample_weight=w['test'])
print '' 
print 'NN score =',score
