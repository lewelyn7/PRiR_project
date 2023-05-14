import os

sizes = [2**i for i in range(10,29, 2)]

for size in sizes:
    for attempt in range(5):
        for thread in [2,4,6,8]:
            print(f"executing with size: {size} attempt: {attempt+1}")
            os.system(f"./main {size} {thread}")
    