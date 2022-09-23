import matplotlib.pyplot as plt
import csv
import numpy as np

if __name__ == '__main__':
    res = []

    with open('test.csv') as file:
        reader = csv.reader(file, delimiter=',')
        res = []
        
        for line in reader:
            l = [m for m in line if m]
            data = (np.array(l)).astype(float)
            data = np.sort(data)
            
            N = 100
            sor = np.zeros((N, ))
            x = np.zeros((N, ))
            dx = 0.05/N
            x[0] = dx
            q = 0
            
            for i in range(len(data)):
                while data[i] > x[q]:
                    q += 1
                    x[q] = x[q - 1] + dx
                    
                sor[q] += 1
                    
            
            x.resize(q)
            sor.resize(q)
            
            res.append((x, sor))
            
    for i in range(len(res)):
        if i%100 == 0:
            plt.plot(res[i][0], res[i][1])
            plt.show()
