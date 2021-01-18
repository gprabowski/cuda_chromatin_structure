import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

chromosomes = ["chr1", "chr14", "chr21"]

def test():
    for chromosome in chromosomes:
        for i in range(50):
            os.system(f"./3dnome_gpu -s ./stg_gpu.ini -c {chromosome} -o ./{chromosome}/")
            # os.system(f"./3dnome_cpu -s ./stg.ini -c {chromosome} -o ./{chromosome}/")

        os.replace("times_gpu.txt", f"{chromosome}/times_gpu.txt")
        os.replace("scores_gpu.txt", f"{chromosome}/scores_gpu.txt")
        # os.replace("times_cpu.txt", f"{chromosome}/times_cpu.txt")
        # os.replace("scores_cpu.txt", f"{chromosome}/scores_cpu.txt")


def stats(file_name: str):
    with open(file_name) as f:
        values = [float(row) for row in f]
        _min = min(values)
        _max = max(values)
        avg = np.mean(values)
        stdev = np.std(values)

        print(f"Min: {_min}, Max: {_max}, Mean: {avg}, Std dev: {stdev}")

def bar_chart():
    c = {
        "cpu": [1308.28172, 44.971254, 2.3483246],
        "gpu": [4.3018849999999995, 0.8185148200000001, 0.28969110000000003]
    }

    j = { x: c[x] for x in ["cpu", "gpu"] }
    df = pd.DataFrame(j)

    # index = np.arange(len(chromosomes))
    index = np.arange(1)
    bar_width = 0.35

    fig, ax = plt.subplots()

    cpu = ax.bar(index, df["cpu"][0], bar_width, label="CPU")
    gpu = ax.bar(index+bar_width, df["gpu"][0], bar_width, label="GPU")

    # ax.set_xlabel('')
    ax.set_ylabel("Execution Time (s)")
    ax.set_title('CPU vs GPU performance comparison')
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(["Chromosome 1"])
    ax.legend()

    def autolabel(rects,data):
        """
        Attach a text label above each bar displaying its height
        """
        c = 0
        initial = 0.091
        offset = 0.205
        use_global_coordinate = False

        if use_global_coordinate:
            for i in data:        
                ax.text(initial+offset*c, 0.05, str(i), horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes,fontsize=8)
                c=c+1
        else:
            for rect,i in zip(rects,data):
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()/2., height, secondsToMinutes(i) ,ha='center', va='bottom')
                
    autolabel(cpu, [df["cpu"][0]])
    autolabel(gpu, [df["gpu"][0]])

    plt.show()

def secondsToMinutes(sec: float):
    mins = math.floor(sec / 60)
    secs = int(sec) % 60
    ms = sec - float(int(sec))

    return "00:{:02d}:{:02d}.{}".format(mins, secs, int(ms * 1000))


def rmsd(file_1: str, file_2: str):
    with open(file_1) as f1:
        with open(file_2) as f2:
            mean = np.mean(np.array([float(row) for row in f1]))
            values = np.array([float(row) for row in f2])
            rmsd = math.sqrt(sum(np.square([x - mean for x in values]) / len(values)))
            print(f"RMSD: {rmsd}")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        test()
    elif len(sys.argv) == 2:
        if(sys.argv[1] == "bar"):
            bar_chart()
        else:
            stats(sys.argv[1])
    else:
        rmsd(sys.argv[1], sys.argv[2])