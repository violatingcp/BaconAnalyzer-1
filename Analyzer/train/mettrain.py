#!/usr/local/bin/python2.7

from sys import exit 
from os import environ
environ['KERAS_BACKEND'] = 'tensorflow'

import numpy as np

#import setGPU
from  convert_flat import process_file
from keras.layers import Input, Dense, Dropout, Activation, Conv2D, MaxPooling2D, LSTM, Convolution1D, MaxPooling1D, MaxPooling1D,BatchNormalization,Reshape
from keras.models import Model,Sequential 
from keras.callbacks import ModelCheckpoint 
from keras.optimizers import Adam 
from keras.utils import np_utils
from keras import backend as K
from keras.preprocessing import sequence
K.set_image_data_format('channels_last')


def reform(arr_x,arr_y,arr_w,train_frac, val_frac, label):
    n = len(arr_x)
    ns = {}
    ns['train'] = (0,int(train_frac*n))
    ns['val']   = (ns['train'][1],ns['train'][1]+int(val_frac*n))
    ns['test']  = (ns['val'][1],n)
    x = {}; y = {}; w = {}; extras = {}
    for subset in ['train','val','test']:
        n_ = ns[subset]
        x[subset] = arr_x[n_[0]:n_[1]]
        y[subset] = arr_y[n_[0]:n_[1]]
        w[subset] = np.ones(n_[1]-n_[0])
        #w[subset] = arr_w[n_[0]:n_[1]]
    return {'x':x,'y':y,'w':w}
    
def load_data(iFile,train_frac, val_frac):
    arr_x = process_file(iFile,1,False,["puppiu","puppiuphi","pfu","pfuphi"])
    arr_y = process_file(iFile,1,True,["zpt","zphi"])
    arr_w = process_file(iFile,1,True,["scale"])
    data = reform(arr_x,arr_y,arr_w,train_frac,val_frac,0)
    x = {}; y = {}; w = {}
    for subset in ['train','val','test']:
        w[subset] = data['w'][subset]
        x[subset] = data['x'][subset]
        y[subset] = data['y'][subset]
    return x,y,w

def save_model(signal=None, frame=None):
    print 'Saving model'
    model.save('model3_pf100.h5')
    print 'Saved model'
    flog.close()
    exit(1)

#def get_adv_loss():
#    def adv_loss(y_true, y_pred,yd_true,yd_pred):
#        pResidDa=fAdDa.predict(y_pred)
#        pResidMc=fAdMc.predict(y_pred)
#        lLossDa =K.mean(K.square(pResidDa-pResidMc),axis=-1)
#        return K.mean(K.square(y_pred - y_true), axis=-1)+lLossDa
#    return adv_loss

#def reformat(iBase,iPred):
    

#def train_adv(iModel,iAdvMC,iAdvDa,x,y,xdata,ydata):
#    d_loss  = iModel.train_on_batch(x, y)
#    ymc     = iModel.predict(x)
#    ydat    = iModel.predict(xdata)
#    ydat    = reformat(ydata,ydat)
#    ymc     = reformat(y,    ymc)
#    am_loss = iAdvMC  .train_on_batch(x,ymc)
#    ad_loss = iAdvDa  .train_on_batch(x,y)

x,y,w    = load_data("ZJets.root"     ,0.5, 0.25)
#xd,yd,wd = load_data("ZJets_data.root",0.5, 0.25)
dim1 = x['train'].shape[1]
print "Test:",x['train'].shape,y['train'].shape
#inputs = Input(shape=(dim1), name='input')
#norm = BatchNormalization(momentum=0.6, name='input_inclusive_bnorm')(inputs)
#conv = Convolution1D(32, 1, padding='valid', activation='relu', input_shape=(dim1,dim2))(norm)
#conv = Convolution1D(16, 1, padding='valid', activation='relu')(conv)
#conv = Convolution1D(8, 1, padding='valid', activation='relu')(conv)
#donv = Convolution1D(4, 1, padding='valid', activation='relu')(inputs)
#conv = Convolution1D(2, 1, padding='valid', activation='relu')(conv)
#l = Dense(32, input_shape=(dim1,1),activation='relu')(inputs)
#l = Dense(64, activation='relu')(l)
#l = Dense(64, activation='relu')(l)
#output = Dense(1,activation='linear',init='normal')(l)
#model = Model(inputs=inputs, outputs=output)
#model.compile(optimizer=Adam(),loss='mean_squared_error',metrics=['accuracy'])

model = Sequential()                                                                                                                                                                                     #model.add(LSTM(1,input_shape=(dim0,dim1)))  
model.add(Dense(32,input_shape=(dim1,)))
model.add(Dense(64))
model.add(Dense(2)) 
model.compile(loss='mean_squared_error', optimizer='adam')   

#advd = Sequential()
#advd.add(self.generator())
#advd.add(self.discriminator())
#advd.compile(loss='mean_squared_error', optimizer=optimizer,metrics=['accuracy'])

#advm = Sequential()
#advm.add(self.generator())
#advm.add(self.discriminator())
#advm.compile(loss='mean_squared_error', optimizer=optimizer,metrics=['accuracy'])

flog = open('train3.log','w')
model.fit(x['train'], y['train'],sample_weight=w['train'],
          #validation_data=(x['val'],y['val'],w['val']),
          batch_size=50, epochs=30, verbose=2)
print model.summary()

y_pred = model.predict(x['test'])
print y_pred[:5]
print y['test'][:5]
print y['train'][:5]

score = model.evaluate(x['test'], y['test'], batch_size=32, verbose=1, sample_weight=w['test'])
print '' 
print 'NN score =',score

save_model()
