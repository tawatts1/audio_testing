#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 13:56:19 2020

@author: ted
"""
from gs_methods import triangle_IC, write_wav, play_wav
from integrators import verlet
from energy_analysis import KE_sys_tot
import numpy as np
import matplotlib.pyplot as plt

class GString:
    def __init__(self, 
                 spring_k, lamb, L, N,
                 F_sys, F_args, 
                 dimensions = 2):
        self.k          = spring_k
        self.lamb       = lamb
        self.L          = L
        self.N          = N
        self.F_sys      = F_sys
        self.F_args     = F_args
        self.dimensions = dimensions
    def calculate(self, 
                  r0v0, wav_rate = 48000, over_calc = 1, T = 5,
                  integrator = verlet, delete_v = True):
        #r0, v0 = r0v0
        print('calculating...')
        self.wav_rate = wav_rate
        self.over_calc = over_calc
        rate = wav_rate*over_calc
        dt = 1/rate
        m = self.lamb * self.L/(self.N-1)
        m_array = np.array([[m for _ in range(self.N)]]*self.dimensions).T
        shape = (int(rate*T), self.N, self.dimensions)
        self.shape = shape
        r = np.zeros(shape) # not r = v = np.zeros(shape) for some reason
        v = np.zeros_like(r)
        r[0], v[0] = r0v0
        r,v = integrator(r, v, dt, self.F_sys, self.F_args, m_array)
        self.r = r
        self.time = np.linspace(0,T,len(self.r))
        plt.scatter(self.time, KE_sys_tot(v,m, particle_axis=1), s = .1)
        plt.title('Total KE over indeces')
        plt.show()
        if delete_v:
            self.v = v[-1]
        else:
            self.v = v
        print('done calculating')
        
    def gen_wav(self, 
                 fname, listening_ratios = [.4,.6], position_indeces = [1,1], 
                 test_mode = False):
        self.fname = fname
        #left index, right index, used for left and right speakers
        li,ri = [int(a*(self.N-1)) for a in listening_ratios] 
        write_wav(fname, self.wav_rate, 
                  self.r[::self.over_calc, li, position_indeces[0]], 
                  self.r[::self.over_calc, ri, position_indeces[1]])
        if test_mode:
            plot_i = [0,1,2,3,4]
        else:
            plot_i = [0,self.shape[0]//2, 3*self.shape[0]//4, -1]
        
        for i in plot_i:
            plt.plot(self.r[i,:,0],self.r[i,:,1], label = f'i= {i}')
        plt.legend()
        plt.show()
    def play_wav(self):
        play_wav(self.fname)
    def plot(self, i_list):
        for i in i_list:
            plt.plot(self.r[i,:,0],self.r[i,:,1], label = f'i= {i}')
        plt.legend()
        plt.show()
        
        
        
        
if __name__ == '__main__':
    from system_forces import F_linear
    # IC check
    k1 = 5e6
    lamb = 1.0
    L = 1
    N = 20
    pluck_x = .2
    pluck_y = 2e-3
    x0_ratio = 0.5
    x0 = x0_ratio*L/(N-1)
    T = 1
    
    IC = triangle_IC(L, N, pluck_x, pluck_y) #shape check
    plt.plot(IC[0][:,0], IC[0][:,1])
    plt.show()
    print(IC[1]==0) # velocity check
    
    
    gs1 = GString(k1, lamb, L, N, 
                  F_linear, [k1, x0])
    gs1.calculate(IC, T = T)
    gs1.gen_wav('test_new_gs.wav', listening_ratios=[0.45, 0.55])
    gs1.play_wav()
    
