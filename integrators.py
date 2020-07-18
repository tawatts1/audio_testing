#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 19:15:14 2020

@author: ted
"""
import numpy as np

def verlet(r,v,dt, F_sys, F_args, m):
    '''
    Velocity verlet integrator (good at conserving energy, order h**2)
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
        vhalf =       v[i] + 0.5*dt * F / m
        
        r[i+1,:,:] =  r[i,:,:] + dt * vhalf
        
        F = F_sys(r[i+1,:,:], *F_args)
        
        v[i+1,:,:] = vhalf + 0.5*dt * F / m
    return r,v

def rk4(r, v, dt, F_sys, F_args, m):
    shape = np.shape(r)
    #convert r, v to one array, 
    #integrate with doctored forces, 
    #convert back
    #np.concatenate((x1, x2), axis = 1)
    for i in range(shape[0] - 1): #for every timestep
        drs = []
        dvs = []
        rdot = v[i]
        vdot = F_sys(r[i,:,:], *F_args) / m
        drs.append(rdot)
        dvs.append(vdot)
        for h in [0.5, 0.5, 1]:
            temp = rdot
            rdot = v[i] + h*dt*vdot #at a slightly updated velocity
            vdot = F_sys(r[i,:,:] + h*dt*temp, *F_args) / m
            drs.append(rdot)
            dvs.append(vdot)
        r[i+1:,:] = r[i,:,:] + dt/6 * (drs[0] + 2*drs[1] + 2*drs[2] + drs[3])
        v[i+1:,:] = v[i,:,:] + dt/6 * (dvs[0] + 2*dvs[1] + 2*dvs[2] + dvs[3])
    return r, v
def rk4_f(r, v, dt, F_sys, F_args, m):
    shape = np.shape(r)
    #convert r, v to one array, 
    #integrate with doctored forces, 
    #convert back
    #np.concatenate((x1, x2), axis = 1)
    for i in range(shape[0] - 1): #for every timestep
        drdt = np.zeros((4, shape[1], shape[2]))
        dvdt = np.zeros_like(drdt)
        #rdot = v[i]
        vdot = F_sys(r[i,:,:], *F_args) / m
        drdt[0] = v[i]
        dvdt[0] = vdot
        for j, h in enumerate([0.5, 0.5, 1]):
            #temp = rdot
            drdt[j+1] = v[i] + h*dt*dvdt[j] #at a slightly updated velocity
            dvdt[j+1] = F_sys(r[i,:,:] + h*dt*drdt[j], *F_args) / m
            #drs.append(rdot)
            #dvs.append(vdot)
        r[i+1:,:] = r[i,:,:] + dt/6 * (drdt[0] + 2*drdt[1] + 2*drdt[2] + drdt[3])
        v[i+1:,:] = v[i,:,:] + dt/6 * (dvdt[0] + 2*dvdt[1] + 2*dvdt[2] + dvdt[3])
    return r, v



def euler(r, v, dt, F_sys, F_args, m):
    '''
    Simple Euler integrator
    
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
    for i in range(shape[0] - 1):
        r[i+1] = r[i] + v[i]*dt
        v[i+1] = v[i] + dt*F_sys(r[i,:,:], *F_args) / m
    return r,v
def enhanced_euler(r, v, dt, F_sys, F_args, m):
    '''
    Enhanced Euler integrator
    
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
    for i in range(shape[0] - 1):
        r[i+1] = r[i] + v[i]*dt
        v[i+1] = v[i] + dt*F_sys(r[i+1,:,:], *F_args) / m
    return r,v
if __name__ == '__main__':
    
    
    pass


