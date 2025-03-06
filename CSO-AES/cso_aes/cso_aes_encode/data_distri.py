import numpy as np
import matplotlib.pyplot as plt
import os
from .mlp_encoding_extract import decode
import multiprocessing

def Freedman_Diaconis_bins(data):
    # 计算四分位数
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1

    # 计算Freedman-Diaconis规则的条柱宽度
    n = len(data)
    bin_width = 2 * IQR / (n ** (1 / 3))

    # 计算组数
    data_range = np.ptp(data)
    number_of_bins = int(np.ceil(data_range / bin_width))
    return number_of_bins

def scott(data):
    std_dev = np.std(data)
    # 计算样本数量
    n = len(data)

    # 使用Scott法则计算带宽
    bin_width_scott = 3.5 * std_dev / (n ** (1 / 3))

    # 计算组数
    data_range = np.ptp(data)  # 数据的范围
    number_of_bins = int(np.ceil(data_range / bin_width_scott))
    return number_of_bins

def distribution(array_data,bw,method,fig_title, body, plot_model):
    D = len(array_data[0])
    zero_freq_intervals_list = []
    fig, axes = plt.subplots(D, 1, figsize=(8, D * 4))  # 宽度为8英寸，高度为D * 4英寸
    fig.suptitle(f'body_{body} distribution of type_{fig_title} (method_{method})')
    max_min = []
    bins = []
    array_data = np.array(array_data)

    for i in range(D):
        new_data = array_data[:, i]
        max_min.append([max(new_data),min(new_data)])
        if method =='Freedman_Diaconis':
            bin = Freedman_Diaconis_bins(new_data)
        elif method =='self_input':
            bin = int(np.ceil((max(new_data)-min(new_data))/bw))
        elif method == 'scott':
            bin = scott(new_data)
        elif method == 'std':
            bin = int(np.ceil((max(new_data)-min(new_data))/(np.std(new_data)/10)))
        else:
            raise ValueError("method no exsit!")
        bins.append(bin)
        frequencies, bin_edges, patches = axes[i].hist(new_data, bins=bin, alpha=0.6, color='g',edgecolor='black')  # 绘制直方图
        zero_freq_intervals = [[bin_edges[i], bin_edges[i + 1]] for i in range(len(bin_edges) - 1) if frequencies[i] == 0]
        zero_freq_intervals_list.append(zero_freq_intervals)
        axes[i].set_title(f'Dimension {i}')
        axes[i].set_xlabel('Value')
        axes[i].set_ylabel('Frequency')
        # print(frequencies)
        #print(len(bins),len(frequencies))
    plt.tight_layout()  # 调整子图间距
    if plot_model==True:
        plt.savefig(f'body_{body} distribution of type_{fig_title}',dpi=300)
    plt.close()
    #plt.show()
    return zero_freq_intervals_list, max_min, bins

def worker(args):
    type_atoms, bw, method, fig_title, body, plot_model = args
    stru_temp = [atom[:-1] for atom in type_atoms]
    tt = np.array(stru_temp)
    zero_freq_intervals_list, max_min, bins = distribution(tt, bw, method, fig_title, body, plot_model)
    return zero_freq_intervals_list, max_min, bins

def data_base_distribution(data_base_data, bw, method, body, plot_model):
    train_data = decode(data_base_data)
    large_zero_freq_intervals_list = []
    large_max_min = []
    large_bins = []

    params = [(type_atoms, bw, method, str(type), body, plot_model) for type, type_atoms in enumerate(train_data)]

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    results = pool.map(worker, params)

    # 关闭进程池，不再接受新的任务
    pool.close()

    # 等待所有进程完成
    pool.join()

    # 将结果收集到列表中
    for zero_freq_intervals_list, max_min, bins in results:
        large_zero_freq_intervals_list.append(zero_freq_intervals_list)
        large_max_min.append(max_min)
        large_bins.append(bins)

    return large_zero_freq_intervals_list, large_max_min, large_bins

