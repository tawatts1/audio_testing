#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 08:19:45 2020

@author: ted
"""
from gs_class import GString
from gs_methods import triangle_IC
from system_forces import F_linear, F_quadratic, F_constant
from integrators import verlet, euler, enhanced_euler, rk4, rk4_f

def make_and_play(k, lamb, L,
                  N, pluck_x, pluck_y,
                  x0_ratio, T, F, 
                  fname, overcalc = 1, integr = verlet,
                  testing = False):
    x0 = x0_ratio*L/(N-1)
    IC = triangle_IC(L, N, pluck_x, pluck_y) 
    gs1 = GString(k, lamb, L, N, 
                  F, [k, x0])
    gs1.calculate(IC, T = T, over_calc = overcalc, integrator = integr)
    gs1.gen_wav(fname+
                f'{k:.1e}_{lamb}_{L}_{N}_{pluck_x}_{pluck_y}_{x0_ratio}_{T}_{overcalc}.wav',
                listening_ratios=[0.45, 0.55], test_mode = testing)
    gs1.play_wav()
if __name__ == '__main__':
    '''
    make_and_play(5e6, 1, 1,
                  30, .3, .01, 
                  0.5, 1, F_linear,
                  'sounds/linear_EEuler_', overcalc = 1, integr = enhanced_euler, testing = True) # makes a beautiful KE plot
    '''
    make_and_play(5e6, 1, 1,
                  30, .3, .01, 
                  0.5, 1, F_linear,
                  'sounds/linear_verlet_', overcalc = 1, integr = verlet, testing = True) # makes a beautiful KE plot
    
    
    '''
    make_and_play(5e6, 1, 1,
                  15, .3, .1, 
                  0, 10, F_quadratic,
                  'sounds/quadratic_', overcalc = 5)
    
    make_and_play(5e6, 1, 1,
                  15, .1, .1, 
                  0, 10, F_quadratic,
                  'sounds/quadratic_', overcalc = 1)
    
    make_and_play(5e7, 1, 1,
                  10, .4, .01, 
                  0, 1, F_constant,
                  'sounds/constant_', overcalc = 10, testing  = False)
    '''
    
    
