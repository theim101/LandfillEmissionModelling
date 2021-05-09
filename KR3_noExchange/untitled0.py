#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:56:03 2020

@author: theimovaara
"""
import numpy as np
#import scipy.integrate as spi
import scipy.stats as stats
import pandas as pd
import scipy.special as spspec

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


b1 = 10
b2 = 12

c = 20

V =np.arange(100)

f1 = b1 * (V>=c) + (b2 + (b1-b2)/c * V)*(V<c)



L = V * f1


plt.plot(V,L)