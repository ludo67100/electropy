# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 14:14:18 2018

This class loads HdF5 recordings from MCS acquisition system as matrices of shape ((channel,data))

Allows to load Raw signals, spike-filtered signals or LFP signals (one function for each)
+ associated time vectors
+ associated sampling rates

All in Volts and Seconds 

Hope it will work 

Then all you have to do is to load HdF5IO from eletroPy package; init class with smthg = HdF5IO(filepath)

After that u can load every instance with associated function, they are all described bellow. 

@author: ludovic.spaeth
"""

class HdF5IO:
    
    def __init__(self,filepath):
        import h5py as h5
        file_ = h5.File(filepath)
        
        self.file = file_['Data'] #Loads first node 

#----------RAW RECORDINGS------------------------------------------------------        
    def raw_record(self): #Gets Raw Records as matrix ((channel,data))
        
        raw = self.file['Recording_0']['AnalogStream']['Stream_1']['ChannelData']
        
        import numpy as np 
        raw_record = np.zeros((raw.shape[0],raw.shape[1]))
        raw_conv = float(self.file['Recording_0']['AnalogStream']['Stream_1']['InfoChannel'][0][9]) #Gain 
        
        for i in range(raw.shape[0]): #Stores data in new matrix 
            raw_record[i,:] = raw[i,:]/raw_conv #From pV to V
    
        return raw_record
    
    def raw_time(self): #Gets time vector for raw records 
        import numpy as np
        raw_tick = int(self.file['Recording_0']['AnalogStream']['Stream_1']['InfoChannel'][0][8])/1000000.0 #exp6 to pass from us to s
        raw_length = len(self.file['Recording_0']['AnalogStream']['Stream_1']['ChannelData'][0])        
        raw_time = np.arange(0,raw_length*raw_tick,raw_tick)        
        return raw_time
        
    def raw_sampling_rate(self): #Gets sampling rate
        
        raw_tick = float(self.file['Recording_0']['AnalogStream']['Stream_1']['InfoChannel'][0][8])/1000000.0
        
        return 1.0/raw_tick #In Hz
    
    
#---------Spike filtered recordings--------------------------------------------
    def filt_record(self): #Gets filt Records as matrix (channel,data)
        
        filt = self.file['Recording_0']['AnalogStream']['Stream_0']['ChannelData'] 
        
        import numpy as np 
        filt_record = np.zeros((filt.shape[0],filt.shape[1]))
        filt_conv = float(self.file['Recording_0']['AnalogStream']['Stream_0']['InfoChannel'][0][9]) #Gain 

        
        for i in range(filt.shape[0]):
            filt_record[i,:] = filt[i,:]/filt_conv #From pV to V
        
        return filt_record
    
    def filt_time(self): #Gets time vector for raw records 
        import numpy as np
        filt_tick = int(self.file['Recording_0']['AnalogStream']['Stream_0']['InfoChannel'][0][8])/1000000.0 #exp6 to pass from us to s
        filt_length = len(self.file['Recording_0']['AnalogStream']['Stream_0']['ChannelData'][0])        
        filt_time = np.arange(0,filt_length*filt_tick,filt_tick)        
        return filt_time
        
    def filt_sampling_rate(self):
        
        filt_tick = float(self.file['Recording_0']['AnalogStream']['Stream_0']['InfoChannel'][0][8])/1000000.0
        
        return 1.0/filt_tick #In Hz

    
#---------LFP filtered recordings--------------------------------------------
    def LFP_record(self): #Gets filt Records as matrix (channel,data)
        
        lfp = self.file['Recording_0']['AnalogStream']['Stream_2']['ChannelData'] 
        
        import numpy as np 
        lfp_record = np.zeros((lfp.shape[0],lfp.shape[1]))
        lfp_conv = float(self.file['Recording_0']['AnalogStream']['Stream_2']['InfoChannel'][0][9]) #Gain 

        
        for i in range(lfp.shape[0]):
            lfp_record[i,:] = lfp[i,:]/lfp_conv #From pV to V
        
        return lfp_record
    
    def LFP_time(self): #Gets time vector for raw records 
        import numpy as np
        lfp_tick = int(self.file['Recording_0']['AnalogStream']['Stream_2']['InfoChannel'][0][8])/1000000.0 #exp6 to pass from us to s
        lfp_length = len(self.file['Recording_0']['AnalogStream']['Stream_2']['ChannelData'][0])        
        lfp_time = np.arange(0,lfp_length*lfp_tick,lfp_tick)        
        return lfp_time
        
    def LFP_sampling_rate(self):
        
        lfp_tick = float(self.file['Recording_0']['AnalogStream']['Stream_2']['InfoChannel'][0][8])/1000000.0
        
        return 1.0/lfp_tick #In Hz
        
        


##----------------TRIAL---------------------------------------------
#        
#    
#filepath = r"H:\Federica\2018-08-10T13-54-41McsRecording.h5"
#
#f = HdF5IO(filepath)
#
#
#from matplotlib import pyplot as plt 
#
#fig, ax = plt.subplots(16,1)
#
#for i in range(16):
#    
#    ax[i].plot(f.LFP_time(),f.LFP_record()[i])
    

