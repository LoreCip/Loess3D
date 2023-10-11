import numpy as np

np.random.seed(00000)

n, m, l = [30]*3

x, y, z = np.linspace(0, 1, n), np.linspace(0, 1, m), np.linspace(0, 1, l)
X, Y, Z = np.meshgrid(x, y, z)

f = X**2 + Y**2 + Z**2
f_ran = np.random.normal(f, 0.1)

with open('data.dat', 'w') as f:

    f.write(f"{n}\t{m}\t{l}\n")
    for arr in [X, Y, Z, f_ran]:
        arr = arr.flatten()
        for i in range(len(arr)):
            f.write(str(arr[i]) + "\t")

        f.write("\n")
