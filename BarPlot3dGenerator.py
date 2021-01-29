import os
import sys
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm


def prepare_data_for_3dplot():
    options = ['128_4','128_8','128_16','256_4','256_8','256_16','512_4','512_8','512_16']

    for s in ['3', '5']:
        for c in ["chr1", "chr14", "chr21"]:
            for st in ['scores', 'times']:
                with open(f'{os.path.expanduser("~")}/Desktop/ampere_raw/{c}_{st}_ampere_{s}.txt', 'w+') as fw:
                    fw.write(','.join(options) + '\n')
                    for o in options:
                        with open(f'{os.path.expanduser("~")}/Desktop/ampere_raw/{s}/{c}/{st}_gpu_{o}.txt') as f:
                            l = np.array([float(val) for val in f])
                            fw.write('{:.2f},'.format(np.mean(l)))



def generate_plots_3d(file_path: str):
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(projection='3d')

    _x = np.arange(3)
    _y = np.arange(3)

    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    top = []
    with open(file_path) as f:   
        f.readline()
        top = [float(val) for val in f.readline().split(',') if not val == '']
        if not len(top) == 9:
            del top[9]

    bottom = np.zeros_like(top)
    width = depth = 0.5
    cmap = cm.get_cmap('turbo')

    max_height = np.max(top)   # get range of colorbars so we can normalize
    min_height = np.min(top)
    # scale each z to [0,1], and get their rgb values
    rgba = [cmap((k-min_height)/max_height) for k in top]

    ax1.bar3d(x, y, bottom, width, depth, top, shade=True, color=rgba)

    if 'chr21_' in file_path:
        ax1.set_title(f'Chromosome 21')
    elif 'chr14_' in file_path:
        ax1.set_title(f'Chromosome 14')
    else:
        ax1.set_title(f'Chromosome 1')

    ax1.set_ylabel("threads per block")
    ax1.set_xlabel("block multiplier")

    if '_times_' in file_path:
        ax1.set_zlabel("Time (s)")
    else:
        ax1.set_zlabel("Score")

    ax1.set_xticks(_x + width / 2)
    ax1.set_yticks(_y + width / 2)
    ax1.set_yticklabels(["128", "256", "512"])
    ax1.set_xticklabels(["x4", "x8", "x16"])

    plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        generate_plots_3d()
    elif len(sys.argv) == 2:
        generate_plots_3d(sys.argv[1])