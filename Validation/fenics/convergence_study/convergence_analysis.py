import os
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def read_first_line(file):
    with open(file, 'rt') as fd:
        first_line = fd.readline()
        return first_line


def get_dictionary_from_data(all_files, dir):
    my_dict = {}
    for files in all_files:
        if files.endswith(".txt"):
            my_dict[int(re.findall(r'\d+', files)[0])] = float(read_first_line(dir + files))
    return sorted(my_dict.items())


def plot_data(my_dict_liver, my_dict_beam):
    x_liver = []
    y_liver = []
    for i in my_dict_liver:
        x_liver.append(i[0])
        y_liver.append(i[1] * 1000)
    x_beam = []
    y_beam = []
    for i in my_dict_beam:
        x_beam.append(i[0])
        y_beam.append(i[1] * 1000)
    fontsize = 20
    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    ax.plot(x_liver, y_liver, "-o", markersize=15, linewidth=5, label="Liver", color="blue")
    ax2.plot(x_beam, y_beam, "-o", markersize=15, linewidth=5, label="Beam", color="orange")
    ax.set_xlabel('Number of DOFs', fontsize=fontsize)
    # plt.ylabel('$u_{max}$ [mm]', fontsize=fontsize)
    ax.set_ylabel('$u^{liver}_{max}$ [mm]', color='blue', fontsize=fontsize)
    ax2.set_ylabel('$u^{beam}_{max}$ [mm]', color='orange', fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.xaxis.set_ticks(np.arange(0, 560000, 60000))
    # plt.xticks(fontsize=fontsize)
    ax.set_xticks(np.arange(0, 560000, 60000))
    ax.set_xlim([20000, 580000])
    plt.yticks(fontsize=fontsize)
    # plt.yticks(np.arange(10000, 40000, 5000))
    ax2.set_ylim([-27.75, -25.5])
    ax.set_ylim([-4.5, -1])
    ax.set_yticks(np.arange(-4.5, -0.95, 0.39))
    ax.grid()
    fig.legend(prop={'size': fontsize})

    # secax = ax.secondary_yaxis('right')
    # # secax.xaxis.set_minor_locator(AutoMinorLocator())
    # secax.set_ylabel('$X_{other}$')
    plt.show()


def main():
    dir_liver = "./liver/"
    all_files_liver = os.listdir(dir_liver)
    my_dict_liver = get_dictionary_from_data(all_files_liver, dir_liver)
    dir_beam = "./beam/"
    all_files_beam = os.listdir(dir_beam)
    my_dict_beam = get_dictionary_from_data(all_files_beam, dir_beam)
    plot_data(my_dict_liver, my_dict_beam)


if __name__ == '__main__':
    main()
