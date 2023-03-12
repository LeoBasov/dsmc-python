import numba
import time
import numpy as np

def func(v):
    N = len(v)
    
    for i in range(N):
        v[i] += v[i]
        
    return v

@numba.njit(cache=True)
def func_njit(v):
    N = len(v)
    
    for i in range(N):
        v[i] += v[i]
        
    return v

@numba.njit(numba.float64[:, :](numba.float64[:, :]), cache=True)
def func_njit_sig(v):
    N = len(v)
    
    for i in range(N):
        v[i] += v[i]
        
    return v

def time_f(f, v):
    t = time.time()
    f(v)
    return time.time() - t

if __name__ == '__main__':   
    N = 100000
    v = np.random.random((N, 3))
    
    t = time_f(func, v)
    t_njit = time_f(func_njit, v)
    t_njit_sig = time_f(func_njit_sig, v)
    
    print("base line:                    {:.3e} ".format(t/t))
    print("njit / base line:             {:.3e} ".format(t_njit/t))
    print("njit + signature / base line: {:.3e} ".format(t_njit_sig/t))
    print('done')