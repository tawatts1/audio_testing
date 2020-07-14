#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 19:15:14 2020

@author: ted
"""
import numpy as np

def verlet(r,v,dt, F_sys, F_args, m):
    '''
    Parameters
    ----------
    r : numpy array
        3d array: (ith time, jth particle, kth dimension)
        Make sure the 0th time is already initialized when passed and the 
        rest are zeros. 
        
    v : numpy array
        3d array of velocity, like r. 
    dt : float
        timestep
    F_sys : Function
        Force as a function of position, where the position array is the 
        first argument
    F_args : list-like
        other arguments for F_sys
    m : 1-d array of masses
        

    Returns
    -------
    r : Numpy array of history of positions
    v : Numpy array of history of velocities
    '''
    shape = np.shape(r)
    F = F_sys(r[0,:,:], *F_args)
    for i in range(shape[0] - 1):
        vhalf = v[i] +0.5*dt*F / m
        r[i+1,:,:] = r[i,:,:] + dt * vhalf
        F = F_sys(r[i+1,:,:], *F_args)
        v[i+1,:,:] = vhalf + dt * F / m
    return r,v

if __name__ == '__main__':
    
    
    pass


