import numpy as np
import pickle
from joblib import dump, load
import zlib
import copy
import os

def decode(data_pkl):
    with open(data_pkl, 'rb') as f:
        compressed_data = f.read()
        decompressed_data = zlib.decompress(compressed_data)
        loaded_list = pickle.loads(decompressed_data)
    return loaded_list

#获取two_body,three_body......的对应的列表
def mtp_many_body_list(mtp_type):
    if mtp_type == 'l2k2.mtp':
        dic = {'<0>':[0,10],'<11>':[27,28,29],'<22>': [30,31,32],
               '<211>':[72,73,80,81],'<222>':[91,92,99,100]}
    elif mtp_type == 'l2k3.mtp':
        dic = {'<0>': [0, 10, 20], '<11>': [46, 47, 48, 49, 50, 51], '<22>': [52, 53, 54, 55, 56, 57],
               '<211>': [178,179,180,187,188,189,196,197,198], '<222>': [211,212,213,220,221,222,229,230,231]}
    else:
        raise ValueError("mtp_type does not exist! If you want to add, modify the program!")

    two_body = []
    three_body = []
    four_body = []

    for key,value in dic.items():
        if len(key)-1 == 2:
            two_body += value
        if len(key)-1 == 3:
            three_body += value
        if len(key)-1 == 4:
            four_body += value
    return two_body,three_body,four_body

def alpha_moment_mapping(hyx_mtp_path):
    with open(hyx_mtp_path,'r') as f:
        lines = f.readlines()
    alpha_moment_mapping_list = []
    for line in lines:
        if 'alpha_moment_mapping' in line:
            alpha_moment_mapping_str = line
            start = alpha_moment_mapping_str.index('{') + 1
            end = alpha_moment_mapping_str.index('}')
            content = alpha_moment_mapping_str[start:end]
            # 将内容转换为列表
            alpha_moment_mapping_list = [int(num.strip()) for num in content.split(',')]
    if len(alpha_moment_mapping_list) == 0:
        raise ValueError("len(alpha_moment_mapping_list) = 0!")
    return alpha_moment_mapping_list

def extract_mtp_many_body_index(mtp_type,hyx_mtp_path):
    two_body,three_body,four_body = mtp_many_body_list(mtp_type)
    tt = alpha_moment_mapping(hyx_mtp_path)
    two_body = [tt.index(a) for a in two_body]
    three_body = [tt.index(a) for a in three_body]
    four_body = [tt.index(a) for a in four_body]
    return two_body,three_body,four_body

def des_out2pkl(des_out_path, prefix, num_ele, mtp_type, hyx_mlp_path, body_name_list,out_path):
    length = num_ele
    type_list = [a for a in range(length)]
    total_list = [[] for _ in range(length)]

    two_body_list = copy.deepcopy(total_list)
    three_body_list = copy.deepcopy(total_list)
    four_body_list = copy.deepcopy(total_list)

    two_body, three_body, four_body = extract_mtp_many_body_index(mtp_type, hyx_mlp_path)

    string = '#start'
    with open(des_out_path,'r') as f:
        lines = f.readlines()
    atom_index = 0
    for index,line in enumerate(lines):
        if string in line:
            atom_num = int(line.split()[1])
            #print(atom_num)
            #print(index)
            atom_lines = lines[index+1:index+atom_num+1]
            for atom_line in atom_lines:
                a = atom_line.split()[1:]
                a = [float(i) for i in a]
                a.append(atom_index)
                type = int(atom_line.split()[0])
                for index, type_ori in enumerate(type_list):
                    if type == type_ori:
                        total_list[index].append(a)
                        two_body_temp = [a[i] for i in two_body]
                        two_body_temp.append(a[-1])
                        three_body_temp = [a[i] for i in three_body]
                        three_body_temp.append(a[-1])
                        four_body_temp = [a[i] for i in four_body]
                        four_body_temp.append(a[-1])

                        two_body_list[index].append(two_body_temp)
                        three_body_list[index].append(three_body_temp)
                        four_body_list[index].append(four_body_temp)
            atom_index += 1

    # with open('coding_zlib.pkl', 'wb') as f:
    #     f.write(compressed_data)

    body_list = [two_body_list,three_body_list,four_body_list]
    body_name = [prefix+'_two_body_',prefix+'_three_body_',prefix+'_four_body_']
    dic = {'two':0,'three':1,'four':2}
    body_index = [dic[a] for a in body_name_list]

    for index, (body,name) in enumerate(zip(body_list,body_name)):
        serialized_data = pickle.dumps(body)
        compressed_data = zlib.compress(serialized_data)
        # print(total_list[0])
        if index in body_index:
            with open(os.path.join(out_path, name + 'coding_zlib.pkl'), 'wb') as f:
                f.write(compressed_data)

if __name__ == '__main__':
    hyx_mtp_path = 'hyx.mtp'
    mtp_type = 'l2k2.mtp'
    two_body, three_body, four_body = extract_mtp_many_body_index(mtp_type, hyx_mtp_path)
    print(two_body,three_body,four_body)
    des_out_path = 'md.out'
    prefix = 'md'
    ele = ['O','1','2']
    body_list = ['two']
    out_path = os.getcwd()
    des_out2pkl(des_out_path, prefix, len(ele), mtp_type, hyx_mtp_path, body_list,out_path)


