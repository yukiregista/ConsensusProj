import psutil
import threading
import time

def monitor_cpu(initial_time,event):
    print("START monitor_cpu")
    while not event.wait(0.5):
        elapsed_time = time.time() - initial_time
        cpu_percent = psutil.cpu_percent(percpu=True)
        mem = psutil.virtual_memory() 

        cpu_percent = '\t'.join(["{:10.4f}".format(v) for v in cpu_percent])
        print("time:", int(elapsed_time))
        print(f"Used Memory: {mem.used / (1024 ** 3):.2f} GB")
        print(f"Available Memory: {mem.available / (1024 ** 3):.2f} GB")
        print(f"Total Memory: {mem.total / (1024 ** 3):.2f} GB")
        print(f"Memory percent: {mem.percent}")
    print("END monitor_cpu")
    


