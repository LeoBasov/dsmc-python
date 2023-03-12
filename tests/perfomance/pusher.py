import numba
import numpy as np
import time

def _push(velocities, positions, dt):
    old_positions = np.copy(positions)
    for p in numba.prange(len(positions)):
        positions[p] = positions[p] + velocities[p]*dt
    return (velocities, positions, old_positions)

def push(velocities, positions, dt):
    N = len(positions)
    for p in numba.prange(N):
        positions[p] = positions[p] + velocities[p]*dt
    return positions

@numba.njit(cache=True)
def push_njit(velocities, positions, dt):
    N = len(positions)
    for p in numba.prange(N):
        positions[p] = positions[p] + velocities[p]*dt
    return positions

@numba.njit(numba.float64[:, :](numba.float64[:, :], numba.float64[:, :], numba.float64), cache=True)
def push_njit_sig(velocities, positions, dt):
    N = len(positions)
    for p in numba.prange(N):
        positions[p] = positions[p] + velocities[p]*dt
    return positions

@numba.njit(numba.float64[:, :](numba.float64[:, :], numba.float64[:, :], numba.float64), cache=True, parallel=True)
def push_njit_sig_par(velocities, positions, dt):
    N = len(positions)
    for p in numba.prange(N):
        positions[p] += velocities[p]*dt
    return positions

def time_f(f, v, p, dt):
    t = time.time()
    f(v, p, dt)
    return time.time() - t

if __name__ == '__main__':   
    N = 100000
    v = np.random.random((N, 3))
    p = np.random.random((N, 3))
    dt = 1e-6
    
    t = time_f(_push, v, p, dt)
    t_new = time_f(push, v, p, dt)
    t_njit = time_f(push_njit, v, p, dt)
    t_njit_sig = time_f(push_njit_sig, v, p, dt)
    t_njit_sig_par = time_f(push_njit_sig_par, v, p, dt)
    
    print("base line:                          {:.3e} ".format(t/t))
    print("new / base line:                    {:.3e} ".format(t_new/t))
    print("new / base line:                    {:.3e} ".format(t_njit/t))
    print("njit + signature / base line:       {:.3e} ".format(t_njit_sig/t))
    print("njit + signature + par / base line: {:.3e} ".format(t_njit_sig_par/t))
    print('done')