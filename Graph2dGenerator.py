import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import numpy as np
import pandas as pd

def size_vs_speedup():

    sizes = [
        (181, 1), (194, 2), (162, 3), (119, 4), (120, 5), (129, 6), (127, 7), (110, 8), (90, 9),
        (108, 10), (96, 11), (110, 12), (59, 13), (58, 14), (65, 15), (59, 16), (72, 17), (55, 18), (53, 19), (48, 20), (23, 21), (37, 22), (97, 23)
    ]

    sorted_sizes = sorted(sizes, key=lambda x: x[0], reverse=True)
    cpu_times = []
    gpu_times = []

    for i in range(23):
        if i == 0:
            cpu_times.append(1308.28172)
            gpu_times.append(4.17)
            continue
        elif i == 13:
            cpu_times.append(44.971254)
            gpu_times.append(0.64)
            continue
        elif i == 20:
            cpu_times.append(2.3483246)
            gpu_times.append(0.13)
            continue


        for device in ['cpu', 'gpu']:                
            if i == 22:
                with open(f'chrX/times_{device}.txt') as cpu:
                    mean = np.mean([float(val) for val in cpu.readlines() if not val == ''])
                    if device == 'cpu':
                        cpu_times.append(mean)
                    else:
                        gpu_times.append(mean)
            else:
                with open(f'chr{i+1}/times_{device}.txt') as cpu:
                    mean = np.mean([float(val) for val in cpu.readlines() if not val == ''])
                    if device == 'cpu':
                        cpu_times.append(mean)
                    else:
                        gpu_times.append(mean)
            
    cpu = np.array(cpu_times)
    gpu = np.array(gpu_times)
    plt.plot([val[0] for val in sorted_sizes], (cpu / gpu) * 100, 'ro')
    plt.ylabel('Speedup in %')
    plt.xlabel('Chromosome size')
    plt.title('Increase in performance with respect to chromosome size')
    plt.show()

def group_bar_chart_cpu_gpu():
    chromosomes = ['chr1', 'chr14', 'chr21']
    labels = ["cpu", "gpu1", "gpu2"]

    data = {
        "cpu": [1308.28172, 44.971254, 2.3483246], 
        "gpu1": [4.17, 0.64, 0.13], 
        "gpu2": [4.08, 0.46, 0.21]
    }

    scores = {
        "cpu": [1833.4105999999997, 142.53581999999997, 15.27742],
        "gpu1": [1774.88, 137.91, 15.58], 
        "gpu2": [1787.7, 159.09, 17.36]
    }

    j = { x: scores[x] for x in labels }
    df = pd.DataFrame(j)

    bar_width = 0.2
    index = np.arange(len(chromosomes))
    fig, ax = plt.subplots()
    colors = plt.cm.BuPu(np.linspace(0.5, 0.9, len(labels)))
    cell_text = []

    cpu = ax.bar(index, df["cpu"], bar_width, label="CPU", color=colors[0])
    gpu1 = ax.bar(index+1*bar_width, df["gpu1"], bar_width, label="GPU1 (Turing)", color=colors[1])
    gpu2 = ax.bar(index+2*bar_width, df["gpu2"], bar_width, label="GPU2 (Ampere)", color=colors[2])

    y_offset = df["cpu"]
    cell_text.append(['%1.2f' % (x) for x in y_offset])
    y_offset = df["gpu1"]
    cell_text.append(['%1.2f' % (x) for x in y_offset])
    y_offset = df["gpu2"]
    cell_text.append(['%1.2f' % (x) for x in y_offset])

    # Add a table at the bottom of the axes
    the_table = ax.table(cellText=cell_text,
                    rowLabels=labels,
                    rowColours=colors,
                    colLabels=chromosomes,
                    loc='bottom')

    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.2, bottom=0.2)

    # Initialize the vertical-offset for the stacked bar chart.
    y_offset = np.zeros(len(chromosomes))

    # ax.set_xlabel('')
    ax.set_yscale('log')
    ax.set_ylabel("Score")
    ax.set_title('CPU vs GPU score comparison')
    ax.set_xticks([])
    # ax.set_xticks(index + len(chromosomes) * bar_width / 2)
    # ax.set_xticklabels(chromosomes)
    ax.legend()

    plt.show()

# accepts list of file names in order: chr1, chr14, chr21
def group_bar_chart_gpu(files: list):
    chromosomes = ['chr1', 'chr14', 'chr21']
    data = dict()
    labels = []

    for i in range(len(files)):
        with open(files[i]) as f:
            labels = [x.split('\n')[0] for x in f.readline().split(',') if not x == '']
            d = [float(x) for x in f.readline().split(',') if not x == '']

            for j in range(len(labels)):
                if i == 0: data[labels[j]] = []
                data[labels[j]].append(d[j])
            
    j = { x: data[x] for x in labels }
    df = pd.DataFrame(j)
    print(df)

    bar_width = 0.1
    index = np.arange(len(chromosomes))
    fig, ax = plt.subplots()
    colors = plt.cm.BuPu(np.linspace(0.5, 0.9, len(labels)))
    cell_text = []

    # Initialize the vertical-offset for the stacked bar chart.
    y_offset = np.zeros(len(chromosomes))

    for i in range(len(labels)):
        ax.bar(index+i*bar_width, df[labels[i]], bar_width, label=labels[i], color=colors[i])
        y_offset = df[labels[i]]
        cell_text.append(['%1.2f' % (x) for x in y_offset])

    # Add a table at the bottom of the axes
    the_table = ax.table(cellText=cell_text,
                    rowLabels=labels,
                    rowColours=colors,
                    colLabels=chromosomes,
                    loc='bottom')

    
    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.2, bottom=0.2)

    # ax.set_xlabel('')
    ax.set_ylabel("Score")
    ax.set_title('Turing GPU Thread-Block configuration scores')
    # ax.set_xticks(index + len(labels) * bar_width / 2)
    ax.set_xticks([])
    # ax.set_xticklabels(chromosomes)
    ax.legend()

    plt.show()

if __name__ == "__main__":

    if len(sys.argv) == 1:
        group_bar_chart_cpu_gpu()
    elif len(sys.argv) == 3:
        pass
    elif len(sys.argv) == 4:
        group_bar_chart_gpu([sys.argv[1], sys.argv[2], sys.argv[3]])