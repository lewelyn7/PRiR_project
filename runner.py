import os
import time
import signal
import sys
import enum

class Modes(enum.Enum):
    OMP = 0
    MPI = 1
    HYBRID = 2
    
mode = Modes.OMP

stop = False
def signal_handler(sig, frame):
    global stop
    stop = True
    print('You pressed Ctrl+C!')

signal.signal(signal.SIGINT, signal_handler)

sizes = [2**i for i in range(12,29, 2)]
sizes.append(2**29)
time.sleep(3)

for size in sizes:
    for attempt in range(5):
        for thread in [2,4,6,8]:
            if stop:
                print("exiting")
                exit()
            print(f"executing with size: {size} attempt: {attempt+1}")
            if mode == Modes.MPI: 
                os.system(f"mpirun -n {thread}  --oversubscribe main_mpi {size} {thread} 4096")
            elif mode == Modes.OMP:
                os.system(f"./main {size} {thread} 345")
                
    