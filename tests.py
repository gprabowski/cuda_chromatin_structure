import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def test():
    # for chromosome in ["chr1", "chr14", "chr21"]:
    for c in range(9,23):
        if c == 0 or c == 13 or c == 20 or c == 1:
            continue
        for i in range(10):
            if not c == 22:
                os.system(f"./3dnome_gpu -s ./stg_gpu.ini -c chr{c+1} -o ./chr{c+1}/")
                os.system(f"./3dnome_cpu -s ./stg.ini -c chr{c+1} -o ./chr{c+1}/")
            else:
                os.system(f"./3dnome_gpu -s ./stg_gpu.ini -c chrX -o ./chrX/")
                os.system(f"./3dnome_cpu -s ./stg.ini -c chrX -o ./chrX/")

        os.system(f"mv ./*.txt chr{c+1}/")

def aggregate(path: str):
    times = []
    scores = []
    temp = 0.0

    for fi in sorted(os.listdir(path)):
        if ".txt" in fi:
            with open(path + '/' + fi) as f:
                values = [float(row) for row in f]
                avg = str(np.mean(values))

                print(f'Processing {fi}')
                
                if fi[-6:-4] == '16':
                    temp = avg
                    print("was 16, waiting")
                elif fi[-5:-4] == '8':
                    if 'times' in fi:
                        times.append(avg)
                        times.append(temp)
                    else:
                        scores.append(avg)
                        scores.append(temp)
                    
                    print("was 8, writing both 8 and 16")
                else:
                    if 'times' in fi:
                        times.append(avg)
                    else:
                        scores.append(avg)
                    
                    print("was 4")

    f_times = open(path + "/times_turing.txt", 'w+')
    f_scores = open(path + "/scores_turing.txt", 'w+')

    header = '128_4,128_8,128_16,256_4,256_8,256_16,512_4,512_8,512_16\n'

    f_times.write(header)
    f_scores.write(header)

    f_times.write(','.join(times))
    f_scores.write(','.join(scores))

    f_times.close()
    f_scores.close()

def stats(path: str):
    for fi in os.listdir(path):
        if ".txt" in fi:
            with open(path + '/' + fi) as f:
                values = [float(row) for row in f]
                _min = min(values)
                _max = max(values)
                avg = np.mean(values)
                stdev = np.std(values)

                print(f"====== {fi} ========")
                print(f"Min: {_min}, Max: {_max}, Mean: {avg}, Std dev: {stdev}")
                print()

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
            # aggregate(sys.argv[1])
    else:
        rmsd(sys.argv[1], sys.argv[2])