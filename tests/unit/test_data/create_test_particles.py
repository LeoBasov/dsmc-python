import numpy as np

if __name__ == "__main__":
    N = 100
    print("writing {} particles".format(N))
    with open("particles.csv", "w") as file:
        for i in range(N):
            pos = np.random.random(3)*2.0 - np.ones(3)
            file.write("{},{},{}\n".format(pos[0], pos[1], pos[2]))
            
    print("done")