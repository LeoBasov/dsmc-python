import matplotlib.pyplot as plt
import csv
import numpy as np
import argparse

def read_data(file_name):
    with open(file_name) as file:
        reader = csv.reader(file, delimiter=',')
        
        for line in reader:
            l = [m for m in line if m]
            data = [float(l[i]) for i in range(len(l))]
            
    return data

def read_ref_data(file_name):
    x = []
    val = []
    with open(file_name) as file:
        reader = csv.reader(file, delimiter=',')
    
        for line in reader:
            x.append(float(line[0]))
            val.append(float(line[1]))
            
    return x, val

if __name__ == '__main__':
    n = read_data("n.csv")
    p = read_data("p.csv")
    T = read_data("T.csv")
    
    x_n, n_ref = read_ref_data("ref_data/sod_n.csv")
    p_n, p_ref = read_ref_data("ref_data/sod_p.csv")
    T_n, T_ref = read_ref_data("ref_data/sod_T.csv")
            
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))
    
    ax0.plot(np.linspace(0, 0.1, len(n)), n)
    ax0.plot(x_n, n_ref)
    ax0.set_xlim([0, 0.1])
    ax0.set_ylim([0, 3e+22])
    ax0.set_title("number denisty")
    ax0.set_xlabel("L / m")
    ax0.set_ylabel("n / m^(-3)")
    
    ax1.plot(np.linspace(0, 0.1, len(n)), T)
    ax1.plot(T_n, T_ref)
    ax1.set_xlim([0, 0.1])
    ax1.set_ylim([150, 500])
    ax1.set_title("temperature")
    ax1.set_xlabel("L / m")
    ax1.set_ylabel("T / K")
    
    ax2.plot(np.linspace(0, 0.1, len(n)), p)
    ax2.plot(p_n, p_ref)
    ax2.set_xlim([0, 0.1])
    ax2.set_ylim([0, 120])
    ax2.set_title("pressure")
    ax2.set_xlabel("L / m")
    ax2.set_ylabel("p / Pa")
    
    plt.show()
