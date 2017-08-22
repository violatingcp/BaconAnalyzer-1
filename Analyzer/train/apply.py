#!/usr/local/bin/python2.7

from sys import exit 
from os import environ, system
import ROOT as r
environ['KERAS_BACKEND'] = 'tensorflow'

import numpy as np

from keras.layers import Input, Dense, Dropout, Activation, concatenate, Conv2D, MaxPooling2D, Conv2DTranspose, LSTM, Convolution1D, MaxPooling1D, MaxPooling1D
from keras.models import Model, load_model 
from keras.callbacks import ModelCheckpoint 
from keras.optimizers import Adam 
from keras.utils import np_utils
from array import array

model = load_model('model3_pf100.h5')
from keras import backend as K
K.set_image_data_format('channels_last')

def setupAddress(iVars,iTree,iArr):
    for pArr in iArr:
        pVar = array( 'f', [ 0. ] )
        iVars.append(pVar)
        iTree.SetBranchAddress(pArr,iVars[len(iVars)-1])

lFile  = r.TFile("ZJets.root")
lTree  = lFile.Get("Events")
lArr   = ["puppiu","puppiuphi","pfu","pfuphi"]
lVars  = []
lCheck = ["zpt","zphi"]
lCVars = []
lPred  = ["u","phi"]
lPVars = []

setupAddress(lVars, lTree,lArr)
setupAddress(lCVars,lTree,lCheck)

lOFile = r.TFile("ZJets_out.root","RECREATE")
lFTree = r.TTree("Events","Events")
for i0 in range(len(lArr)):
    lFTree.Branch(lArr[i0],lVars[i0],lArr[i0]+"/F")

for i0 in range(len(lCheck)):
    lFTree.Branch(lCheck[i0],lCVars[i0],lCheck[i0]+"/F")

for pPred in lPred:
    pVar = array( 'f', [ 0. ] )
    lPVars.append(pVar)
    lFTree.Branch("dnn_"+pPred,lPVars[len(lPVars)-1],"dnn_"+pPred+"/F")

for i0 in range(lTree.GetEntriesFast()):
    lTree.GetEntry(i0)
    if i0 % 100 != 0:
        continue
    if i0 % 10000 == 0:
        print i0,float(i0)/float(lTree.GetEntriesFast())
    pVar =[]
    pVars=[]
    for i1 in range(len(lVars)):
        pVar.append(lVars[i1][0])
    pVars.append(pVar)
    pVars=np.array(pVars)
    pred = model.predict(pVars)
    for i1 in range(len(pred[0])):
        lPVars[i1][0] = pred[0][i1]
    lFTree.Fill()

lOFile.cd()
lFTree.Write()


