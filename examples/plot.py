import matplotlib.pyplot as plt
import csv
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('file_name', type=str)

    args = parser.parse_args()

    with open(args.file_name) as file:
        reader = csv.reader(file, delimiter=',')
        
        for line in reader:
            l = [m for m in line if m]
            n = [float(l[i]) for i in range(len(l))]
            
        plt.plot(n)
        plt.show()
