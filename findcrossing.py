# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 15:30:22 2021

@author: esoria
"""
def matching(x,xp, s_vector):
    from scipy.interpolate import interp1d
    import numpy as np
    from scipy.interpolate import InterpolatedUnivariateSpline
    not_equal = np.array(np.where(np.diff(x) != 0))
    #not_equal = not_equal[0]
   # indexes = np.concatenate(not_equal, not_equal[-1]+1)
    indexes=not_equal
    #u = interp1d(x[indexes],xp[indexes],kind='cubic',fill_value="extrapolate")(s_vector)
    s= InterpolatedUnivariateSpline(x[indexes],xp[indexes],k=3)
    u= s(s_vector)
    return u 
'''
    from scipy.interpolate import interp1d
    import numpy as np
    from scipy.interpolate import InterpolatedUnivariateSpline
    not_equal = np.array(np.where(np.diff(x) != 0),)
   # indexes=np.append(not_equal, not_equal[-1]+1)
    if np.size(not_equal)> np.size(x)-1:
       indexes=not_equal
    else:
        indexes=not_equal[1:]
       # indexes = np.concatenate((not_equal[0,:], [not_equal[0,-1]+1]),axis=0, out=None)
    #u = interp1d(x[indexes],xp[indexes],kind='cubic',fill_value="extrapolate")(s_vector)
    s= InterpolatedUnivariateSpline(x[indexes],xp[indexes],k=3)
    u= s(s_vector)
    return u 
'''
