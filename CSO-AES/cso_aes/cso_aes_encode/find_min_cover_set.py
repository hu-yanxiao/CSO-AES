from multiprocessing import Pool
from collections import defaultdict
from collections import Counter
import time
import random
import math

def cover_count(sets,elem):
    cover_count = sum(1 for s in sets if elem in s)
    return cover_count

def sort_by_count_with_counter(lists):
    temp_lists = []
    for a in lists:
        temp_lists += a
    tt = Counter(temp_lists)
    return sorted(tt, key=tt.get, reverse=True)

def list2intersection(list1,list2):
    # 将列表转换为集合
    set1 = set(list1)
    set2 = set(list2)
    # 计算两个集合的交集
    intersection = set1.intersection(set2)
    return intersection

def set_intersection(list1,list2):
    # 将列表转换为集合
    set1 = set(list1)
    set2 = set(list2)
    # 计算两个集合的交集
    intersection = set1.intersection(set2)
    # 计算交集集合的元素个数
    count = len(intersection)
    return count

def find_min_cover_set(lists):
    cover_set = set()
    sets = [set(a) for a in lists]

    # 使用字典来建立元素到其覆盖集合的映射
    element_to_sets = defaultdict(set)
    for index, s in enumerate(sets):
        for elem in s:
            element_to_sets[elem].add(index)
    sorted_elements = sort_by_count_with_counter(lists)
    element_to_sets = {a: element_to_sets[a] for a in sorted_elements}

    # 临时集合，用于存放所有未覆盖的集合的索引
    remaining_sets = set(range(len(sets)))

    while remaining_sets:
        # 找到覆盖最多集合的元素
        best_elem = None
        best_cover_count = 0
        for elem, covered_sets in element_to_sets.items():
            # 统计这个元素覆盖的集合数
            cover_count = len(covered_sets & remaining_sets)
            if cover_count > best_cover_count:
                best_cover_count = cover_count
                best_elem = elem

        # 更新覆盖集合
        cover_set.add(best_elem)

        # 移除已被覆盖的集合
        remaining_sets -= element_to_sets[best_elem]

        # 从字典中移除已覆盖的元素
        for elem in list(element_to_sets.keys()):
            if remaining_sets.isdisjoint(element_to_sets[elem]):
                del element_to_sets[elem]

    return list(cover_set)

# def fwss(lists, min_cover_index, num):
#     # 使用字典来统计每个元素的出现次数，避免多次遍历all_elements
#     element_count = {}
#     for s in lists:
#         for elem in s:
#             element_count[elem] = element_count.get(elem, 0) + 1
#     # 使用字典来存储元素的计数，避免在列表中查找元素的索引
#     count_list = [element_count[a] for a in min_cover_index]
#     # 根据元素的计数对min_cover_index进行排序
#     sorted_min_cover_index = sorted(min_cover_index, key=lambda x: count_list[min_cover_index.index(x)], reverse=True)
#     # 选择前num个元素
#     select_sorted_min_cover_index = sorted_min_cover_index[:num]
#     # 根据元素的计数对计数列表进行排序
#     sorted_count_list = sorted(count_list, reverse=True)
#     return sorted_min_cover_index, sorted_count_list, select_sorted_min_cover_index

# frequency_weight_selection_structure
def fwss(lists, min_cover_index, num):
    all_elements = []
    count_list = []
    for s in lists:
        all_elements.extend(s)
    for a in min_cover_index:
        count = all_elements.count(a)
        count_list.append(count)
    sorted_min_cover_index = sorted(min_cover_index, key=lambda x: count_list[min_cover_index.index(x)], reverse=True)
    aa=num//5  #  top 20%
    bb=num-aa  # remaining 80%
    if bb==0:
        bb=1
    cc=(len(sorted_min_cover_index)-aa)//bb # step
    if cc==0:
        cc=1
    print(f"a={aa};b={bb};c={cc}")
    select_sorted_min_cover_index = sorted_min_cover_index[:aa]+sorted_min_cover_index[aa::cc]
    if num >=len(sorted_min_cover_index):
        select_sorted_min_cover_index = sorted_min_cover_index[:num]
    return sorted_min_cover_index, sorted(count_list, reverse=True), select_sorted_min_cover_index

# frequency_weight_selection_structure + mean_select_index
def fwss_plus_mean_select_index(lists, min_cover_index, num, mean_select_index,mean_coverage_rate,logger):
    all_elements = []
    count_list = []
    for s in lists:
        all_elements.extend(s)
    for a in min_cover_index:
        count = all_elements.count(a)
        count_list.append(count)
    sorted_min_cover_index = sorted(min_cover_index, key=lambda x: count_list[min_cover_index.index(x)], reverse=True)
    temp_sorted_min_cover_index = sorted_min_cover_index
    sorted_min_cover_index = mean_select_index + [item for item in sorted_min_cover_index if item not in mean_select_index]

    set_mean_coverage_rate_value = 90

    if mean_coverage_rate > set_mean_coverage_rate_value:
        aa = num // 5  # top 20%
        bb = num - aa  # remaining 80%
        if bb == 0:
            bb = 1
        cc = (len(sorted_min_cover_index) - aa) // bb  # step
        if cc == 0:
            cc = 1
        print(f"a={aa};b={bb};c={cc}")
        select_sorted_min_cover_index = sorted_min_cover_index[:aa] + sorted_min_cover_index[aa::cc]
        if num >= len(sorted_min_cover_index):
            select_sorted_min_cover_index = sorted_min_cover_index[:num]
    else:
        if len(sorted_min_cover_index) * 0.2<num:
            select_num = math.floor(len(sorted_min_cover_index) * 0.8)
            select_sorted_min_cover_index = sorted_min_cover_index[select_num:]
        else:
            select_sorted_min_cover_index = sorted_min_cover_index[len(sorted_min_cover_index)-num:]
        logger.info(f'mean_coverage_rate is less than {set_mean_coverage_rate_value}%, select the last 20% structure')

    return set_intersection(select_sorted_min_cover_index,mean_select_index), set_intersection(select_sorted_min_cover_index,temp_sorted_min_cover_index), \
           set_intersection(list2intersection(select_sorted_min_cover_index,mean_select_index),list2intersection(select_sorted_min_cover_index,temp_sorted_min_cover_index)),select_sorted_min_cover_index

# 随机生成超大的集合组成的整数列表
def generate_random_sets(num_sets, min_set_size, max_set_size, value_range):
    sets = []
    for _ in range(num_sets):
        set_size = random.randint(min_set_size, max_set_size)
        # 生成一个集合，避免重复元素
        random_set = set(random.sample(range(value_range), set_size))
        sets.append(random_set)
    return sets

if __name__ == "__main__":
    # 设置参数
    num_sets = 2000  # 集合的数量
    min_set_size = 1  # 每个集合的最小大小
    max_set_size = 5000  # 每个集合的最大大小
    value_range = 50000  # 元素的范围

    random.seed(42)
    # 生成随机集合
    class_sets = generate_random_sets(num_sets, min_set_size, max_set_size, value_range)
    # print(random_sets)
    # class_sets = [
    #     [1, 2, 2, 2, 2, 2, 2],
    #     [2, 3, 4, 4, 4, 4],
    #     [3, 4, 5, 4, 5, 6],
    #     [6, 7, 5, 5],
    #     [8, 9, 8, 8, 8, 8],
    #     [9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 8, 8, 8],
    #     [11, 12, 13, 13],
    #     [12, 13, 13],
    # ]
    start = time.time()
    result = find_min_cover_set(class_sets)
    end = time.time()
    print(f"时间:{end - start}，最小覆盖集合:{result}")

    start = time.time()
    print(fwss(class_sets, result,10))
    end = time.time()
    print("时间:", end - start)



