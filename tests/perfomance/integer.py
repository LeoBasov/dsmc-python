import numba
import time
import numpy as np

def func(N):
    res = 0
    
    for i in range(N):
        res += i
        
    return res

@numba.njit(cache=True)
def func_njit(N):
    res = 0
    
    for i in range(N):
        res += i
        
    return res

@numba.njit(numba.int32(numba.int32), cache=True)
def func_njit_sig(N):
    res = 0
    
    for i in range(N):
        res += i
        
    return res

def time_f(f, N):
    t = time.time()
    f(N)
    return time.time() - t

if __name__ == '__main__':    
    N = 10000000
    
    t = time_f(func, N)
    t_njit = time_f(func_njit, N)
    t_njit_sig = time_f(func_njit_sig, N)
    
    print("base line:                    {:.3e} ".format(t/t))
    print("njit / base line:             {:.3e} ".format(t_njit/t))
    print("njit + signature / base line: {:.3e} ".format(t_njit_sig/t))
    print('done')