#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 14:16:27 2020

@author: ted
"""
import numpy as np
from scipy.io.wavfile import read, write
from pydub import AudioSegment
from pydub.playback import play

def triangle_IC(L, N, pluck_x, pluck_y, dimensions = 2):
    r0 = np.zeros((N,dimensions))
    r0[:,0] = np.linspace(0,L,N)
    for j in range(N):
        ratio = j/(N-1)
        if ratio < pluck_x:
            r0[j,1] = pluck_y/pluck_x * L * ratio
        else: 
            r0[j,1] = pluck_y/(1-pluck_x) * L * (1-ratio)
    v0 = np.zeros_like(r0)
    return r0, v0
    
def scale_data(data):
    M = np.max(np.abs(data))
    return (data/M).astype('float32')
def write_wav(fname, wav_rate, r1, r2 = None):
    if r2 is None:
        r2 = r1
    write(fname, wav_rate, np.array([scale_data(r1), scale_data(r2)]).T )
def play_wav(fname):
    sound = AudioSegment.from_wav(fname)
    play(sound)
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    xy = triangle_IC(10,30,.2,.1)[0]
    plt.scatter(xy[:,0], xy[:,1])

