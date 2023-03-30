# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 13:56:10 2021

@author: vfolomeev
"""

import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

fig, ax = plt.subplots()

u=pd.read_csv('C:/Lanit/Рабочий стол/LBM/Sogol/0u.csv',delimiter=' ')
X=u.values[:,0]
Y=u.values[:,1]
U=u.values[:,2]
V=u.values[:,3]

print(X.size)

print(u)
print(X)
plt.quiver(X,Y,U,V)
#plt.scatter(Y,U)

