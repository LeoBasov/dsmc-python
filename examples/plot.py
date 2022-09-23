import matplotlib.pyplot as plt
import csv
import numpy as np

if __name__ == '__main__':
    res = []

    with open('test.csv') as file:
        reader = csv.reader(file, delimiter=',')
        
        for line in reader:
            l = [m for m in line if m]
            n = [float(l[i]) for i in range(1, len(l))]
            
        plt.plot(n)
        plt.show()
