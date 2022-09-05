import os
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def read_first_line(file):
    with open(file, 'rt') as fd:
        first_line = fd.readline()
        return first_line


def get_dictionary_from_data(all_files):
    my_dict = {}
    for files in all_files:
        if files.endswith(".txt"):
            my_dict[int(re.findall(r'\d+', files)[0])] = float(read_first_line(files))
    return sorted(my_dict.items())


def plot_data(my_dict):
    x = []
    y = []
    for i in my_dict:
        x.append(i[0])
        y.append(i[1])
    fontsize = 20
    plt.plot(x, y, "-o", markersize=15, linewidth=5)
    plt.xlabel('Number of points', fontsize=fontsize)
    plt.ylabel('$u_{max}$ [mm]', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    # plt.xticks(np.arange(10000, 40000, 5000))
    plt.xlim([10000, 35000])
    plt.yticks(fontsize=fontsize)
    # plt.yticks(np.arange(10000, 40000, 5000))
    plt.ylim([0.85, 1.05])
    plt.grid()
    plt.show()


def main():
    all_files = os.listdir("./")
    my_dict = get_dictionary_from_data(all_files)
    plot_data(my_dict)


if __name__ == '__main__':
    main()
