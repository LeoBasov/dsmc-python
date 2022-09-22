import matplotlib.pyplot as plt
import csv
import numpy as np

if __name__ == '__main__':
    res = []

    with open('test.csv') as file:
        reader = csv.reader(file, delimiter=',')
        
        for line in reader:
            l = [m for m in line if m]
            data = (np.array(l)).astype(float)
            data = np.sort(data)
            
            N = 100
            sor = np.zeros((N, ))
            x = np.zeros((N, ))
            step = data[-1] / N
            q = 0
            print(step, data.shape)
            
            for i in range(len(data)):
                if  data[i] < step:
                    x[q] = data[i]
                    sor[q] += 1
                else:
                    q += 1
                    step += step
                    sor[q] += 1
                    x[q] = data[i]
                    
            
            x.resize(q)
            sor.resize(q)
            
            plt.plot(x, sor)
            plt.show()
            plt.hist(data)
            plt.show()
            break

