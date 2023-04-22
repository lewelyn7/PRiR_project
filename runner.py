import os

sizes = [2**i for i in range(10,33, 2)]

for i in sizes:
    for j in range(3):
        print(f"executing with size: {i} attempt: {j+1}")
        os.system(f"./main {i}")
    