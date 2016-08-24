#! /home/jaggu/anaconda/bin/python2.7

from __future__ import division
from scipy.optimize import curve_fit
import numpy as np


def one_exponential(ydata,norm=False):

    def _func(x,a,b,c):
        return a*np.exp(-b*x)+c

    x_data = np.arange( 0,len(ydata),1)
    if norm:
        a0, b0, c0 = 1, 0.02, 1
        y_data = ydata/np.amax(ydata)
    else:
        a0, b0, c0 = 1000,0.02,1000
        y_data = ydata
    
    param_opt,pconv = curve_fit(_func,x_data,y_data,p0=(a0,b0,c0))
    _a,_b,_c = param_opt


    y_fit = _func(x_data,_a,_b,_c)
    rmse = np.sqrt(np.mean((y_data-y_fit)**2))
    sse = np.sum((y_data-y_fit)**2)
    sst = np.sum((y_data-np.mean(y_data))**2)
    r_2 = 1 - sse/sst

    return param_opt,rmse,r_2, y_fit

def two_exponential(ydata,norm=False):

    def _dfunc(x,a,b,d,e):
        return a*np.exp(-b*x)+a*np.exp(-d*x)+e

    x_data = np.arange( 0,len(ydata),1)
    if norm:
        a0,b0,d0,e0 = 0.5,0.01,0.01,0.2
        y_data = ydata/np.amax(ydata)
    else:
        a0,b0,d0,e0 = 500,0.001,0.01,500
        y_data = ydata

    try:
        param_opt,pconv = curve_fit(_dfunc,x_data,y_data,p0=(a0,b0,d0,e0))
    except RuntimeError:
        print "Error: curve fit failed"
        param_opt = (a0,b0,d0,e0)
    _a,_b,_d,_e = param_opt
    
    y_fit = _dfunc(x_data,_a,_b,_d,_e)
    sse = np.sum((y_data-y_fit)**2)
    sst = np.sum((y_data-np.mean(y_data))**2)
    try:
        r_2 = 1 - sse/sst
    except RuntimeWarning:
        r_2 = 0

    rmse = np.sqrt(np.mean((y_data-y_fit)**2))

    return param_opt,rmse,r_2,y_fit

def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

def fit_exp_linear(t, y, C=0):
    y = np.array(y,dtype=float)
    t = np.array(t,dtype=int)
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    y_fit = model_func(t,A,K,C)
    return A, K, y_fit





