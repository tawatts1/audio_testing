#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 17:35:45 2020

@author: ted
"""
import numpy as np

#def distance(r1, r2):
#    rdiff = r2 - r1
#    return np.sqrt(np.sum( (rdiff)**2 ) )
def sys_diff(rs):
    return rs[1:,:] - rs[:-1,:]
def length_ax1(rdiff, axis = 1):
    return np.sqrt(np.sum(rdiff**2, axis = axis))
def unit_vectors(rdiff):
    lengths = length_ax1(rdiff)
    divisor = np.array([lengths for _ in range(np.shape(rdiff)[1])]).T
    return rdiff/divisor
def F_linear(rs, k, x0):
    rdiff = sys_diff(rs)
    rdiff_shape = np.shape(rdiff)
    rdiff_hat = unit_vectors(rdiff)
    F = np.zeros((rdiff_shape[0]+1, rdiff_shape[1]))
    F[1:,:]  -= (rdiff - rdiff_hat*x0)
    F[:-1,:] += (rdiff - rdiff_hat*x0)
    #Fnorm = np.zeros_like(F)
    #Fnorm[:,0] = Fnorm[:,1] = length_ax1(F)
    F = F*k
    F[0,:] = F[-1,:] = np.zeros((rdiff_shape[1]))
    return F

if __name__ == '__main__':
    
    
    pass


