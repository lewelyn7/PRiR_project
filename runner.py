import os
import signal
import sys
stop = False
def signal_handler(sig, frame):
    global stop
    stop = True
    print('You pressed Ctrl+C!')

signal.signal(signal.SIGINT, signal_handler)

sizes = [2**i for i in range(10,29, 2)]

for size in sizes:
    for attempt in range(5):
        for thread in [2,4,6,8]:
            if stop:
                print("exiting")
                exit()
            print(f"executing with size: {size} attempt: {attempt+1}")
            os.system(f"mpirun -n {thread}  --oversubscribe main_new {size} {thread} 1234")
    