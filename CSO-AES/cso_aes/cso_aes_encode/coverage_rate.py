import time
from .mlp_encoding_extract import decode,des_out2pkl
import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
import bisect
from collections import Counter
from matplotlib.colors import LinearSegmentedColormap

def is_in_interval(x, interval):
    """判断 x 是否在区间 interval 内"""
    return interval[0] <= x < interval[1]

def find_interval(x, intervals):
    """使用二分查找找到 x 所在的区间，如果不在任何区间内，返回 True(在零频区间返回False)"""
    # 提取区间的起始点
    starts = [interval[0] for interval in intervals]

    # 使用 bisect_left 找到 x 应该插入的位置
    index = bisect.bisect_left(starts, x)

    max_value = intervals[0][1]
    min_value = intervals[-1][0]

    if x == max_value or x == min_value:
        return True
    else:
        # 检查 x 是否在找到的区间内
        if index < len(intervals) and is_in_interval(x, intervals[index]):
            # return index
            #print(x, intervals[index])
            return False
        elif index > 0 and is_in_interval(x, intervals[index - 1]):
            #print(x, intervals[index - 1])
            # return index - 1
            return False
        else:
            # return -1
            return True

def return_in_zero_freq_intervals_index(data_list, zero_freq_intervals, max_min):
    all_intervals = [[-float('inf'), max_min[1]], *zero_freq_intervals, [max_min[0], float('inf')]]
    label_list = [find_interval(x, all_intervals) for x in data_list]
    #print(sum(label_list))
    return label_list

def md_sub_extract(array_data, type_zero_freq_intervals_list, max_min, coverage_rate_method):
    coverage_rate_list = []

    D = len(array_data[0])
    new_lable_array = np.zeros(len(array_data), dtype=int)
    for i in range(D):
        new_data = array_data[:, i]  # 直接使用NumPy的切片操作
        label_list = return_in_zero_freq_intervals_index(new_data, type_zero_freq_intervals_list[i], max_min[i])
        new_lable_array += np.array(label_list, dtype=int)
        D_coverage_rate = np.sum(label_list) / len(label_list) * 100
        coverage_rate_list.append(D_coverage_rate)

    if coverage_rate_method == 'mean':
        coverage_rate = sum(coverage_rate_list)/len(coverage_rate_list)
    elif coverage_rate_method == 'min':
        coverage_rate = min(coverage_rate_list)
    else:
        raise ValueError("coverage_rate_method has only mean and min!")
    return new_lable_array, coverage_rate

def coverage_rate(data, large_zero_freq_intervals_list, large_max_min, body, plot_out,coverage_rate_method,plot_model):
    md_data = decode(data)
    type_coverage_rate = []
    type_coverage_rate_100 = []
    type_coverage_rate_index = []

    # 使用ProcessPoolExecutor来并行执行md_sub_extract
    with ProcessPoolExecutor() as executor:
        # 准备参数列表，供进程池中的函数使用
        tasks = []
        for type, (type_atoms, type_zero_freq_intervals_list, max_min) in enumerate(
                zip(md_data, large_zero_freq_intervals_list, large_max_min)):
            stru_temp = [(atom[:-1]) for atom in type_atoms]
            tt = np.array(stru_temp)
            type_coverage_rate_100.append(100)
            if len(tt) != 0:
                # 将每个md_sub_extract的参数作为一个元组添加到任务列表中
                tasks.append((tt, type_zero_freq_intervals_list, max_min, coverage_rate_method))
                # 记录数据中有元素的type
                type_coverage_rate_index.append(type)

        # 使用executor.map并行执行md_sub_extract，并收集结果
        results = list(executor.map(md_sub_extract, *zip(*tasks)))

    # 处理结果并绘图
    for type, (result, task) in enumerate(zip(results, tasks)):
        lable_array, coverage_rate = result
        tt, type_zero_freq_intervals_list, max_min, coverage_rate_method = task
        D = len(tt[0])

        plot_data = [tt[:, i] for i in range(tt.shape[1])]

        # # 定义自定义颜色映射
        # colors = ['red'] + [plt.cm.viridis(i / float(D)) for i in range(D)][::-1]  # 将红色作为第一个颜色，后面跟着viridis颜色映射
        # custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, D + 1)
        # scatter = plt.scatter(plot_data[0], plot_data[1], s=0.1,
        #                       c=lable_array, cmap=custom_cmap, vmin=-0.5, vmax=D + 0.5)

        cmap = plt.cm.get_cmap('viridis', D+1)  # 使用viridis颜色映射，并设置为8个等级
        scatter = plt.scatter(plot_data[0], plot_data[1], s=0.2,
                              c=lable_array, cmap=cmap, vmin=-0.5, vmax=D + 0.5)
        # 添加颜色条
        element_counts = Counter(lable_array)
        sorted_counts = sorted(element_counts.items(), key=lambda item: item[0])
        count_str = "\n".join(f"{elem}: {count}" for elem, count in sorted_counts)
        #count_str = count_str + '\n' + str(sum(lable_array) / len(lable_array) / D)

        cbar = plt.colorbar(scatter, ticks=np.arange(0, D + 1, 1))  # 设置颜色条的刻度位置

        # 设置颜色条的刻度标签
        cbar.ax.set_yticklabels(np.arange(0, D + 1, 1))  # 设置刻度标签

        plt.title(f'{body}_body_type_{str(type)} mean_coverage_rate:{round(coverage_rate, 5)}%')
        plt.xlabel('Dimension_0')
        plt.ylabel('Dimension_1')
        plt.text(0.95, 0.95, count_str, transform=plt.gca().transAxes, verticalalignment='top',
                 horizontalalignment='right')
        if plot_model:
            plt.savefig(os.path.join(plot_out, f'{body}_body_type_{str(type)}.png'), dpi=300)
            plt.close()
        type_coverage_rate.append(coverage_rate)

    for a, b in zip(type_coverage_rate_index, type_coverage_rate):
        type_coverage_rate_100[a] = b
    return type_coverage_rate_100


if __name__ == '__main__':
    data_base_data = r"database_two_body_coding_zlib.pkl"
    md_data = decode(r"md_two_body_coding_zlib.pkl")
    type_coverage_rate = coverage_rate(r"md_two_body_coding_zlib.pkl",large_zero_freq_intervals_list, max_min, body, plot_out,coverage_rate_method)
