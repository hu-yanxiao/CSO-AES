import os
from cso_aes_das.gen_while_loop import check_run_position,gen_while_loop,mkdir
from cso_aes_das.file_conversion import cfg2xyz,sort_xyz2cfg
from cso_aes_das.other import touch
import yaml
import glob
from ase.io import iread,write

def check(main_loop_npt,main_loop_nvt):
    if main_loop_npt is not None and main_loop_nvt is not None:
        return main_loop_npt
    elif main_loop_npt is not None and main_loop_nvt is None:
        return main_loop_npt
    elif main_loop_npt is None and main_loop_nvt is not None:
        return main_loop_nvt
    else:
        print('error')

def xyz_para_cfg(path,xyz_input, mtp_type):
    out = os.path.join(path, 'init','original.cfg')
    mtp = os.path.join(path, 'init', mtp_type)
    ori_mtp = os.path.join(path, 'init','original.mtp')
    old_para = os.path.join(path, 'init','old_parameter.yaml')
    para = os.path.join(path, 'init','parameter.yaml')

    ele = sort_xyz2cfg(xyz_input,out)
    with open(old_para, 'rb') as file:
        data = yaml.safe_load(file)
    data['ele'] = ele
    data['ele_model'] = 1
    data['mtp_type'] = mtp_type

    with open(para, 'w') as file:
        yaml.safe_dump(data, file, default_flow_style=False)

    with open(mtp, 'r') as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if 'species_count' in line:
                lines[index] = f'species_count = {len(ele)}' + '\n'
    with open(ori_mtp, "w") as file:
        for line in lines:
            file.write(line)
    #touch(os.path.join(path,'init'),'__ok__')

def main():
    main_loop_npt = None
    main_loop_nvt = None

    main_loop_npt = [[200],[300]]
    #main_loop_nvt = [[50, 100, 200], [300, 400, 500], [600, 700, 800],[900,1000,1100]]

    '''xyz转cfg,同时改变para参数(ele,ele_model,mtp_type)和mtp元素个数'''
    pwd = os.getcwd()
    path = pwd
    xyz_input = 'abacus_LiF_npt.xyz' #绝对路径
    mtp_type = 'l2k2.mtp'
    all_label = 'LiF'
    xyz_para_cfg(path,xyz_input, mtp_type)

    '''cfg最后导出xyz文件'''
    output2xyz = True
    xyz_name = 'all_sample_data.xyz'
    yaml_file = os.path.join(pwd, 'init', 'parameter.yaml')
    with open(yaml_file, 'rb') as file:
        yaml_data = yaml.safe_load(file)
    ele = yaml_data['ele']
    ele_model = yaml_data['ele_model']

    init_threshold = 0.1
    threshold_coff = 1.2

    tt = check(main_loop_npt, main_loop_nvt)

    start_position,gen_num,main_num = check_run_position(pwd,tt)
    sleep_time = 10
    max_gen = 10

    model = 'full-automatic' #semi-automatics

    if model == 'semi-automatic':
        for a in range(main_num,len(tt)):
            try:
                npt = main_loop_npt[a]
            except:
                npt = None
            try:
                nvt = main_loop_nvt[a]
            except:
                nvt = None
            print(npt, nvt, start_position, gen_num, init_threshold, threshold_coff, sleep_time, max_gen)
            gen_while_loop(pwd, npt, nvt, start_position, gen_num, init_threshold, threshold_coff, sleep_time,max_gen)
    elif model == 'full-automatic':
        for a in range(main_num,len(tt)):
            try:
                npt = main_loop_npt[a]
            except:
                npt = None
            try:
                nvt = main_loop_nvt[a]
            except:
                nvt = None
            start_position, gen_num, new_main_num = check_run_position(pwd, tt)
            gen_while_loop(pwd, npt, nvt, start_position, gen_num, init_threshold, threshold_coff, sleep_time,max_gen)
            if new_main_num<len(tt)-1:
                temp = new_main_num + 1
                path = os.path.join(pwd,'main_'+str(temp),'gen_0')
                mkdir(path)

        if output2xyz == True:
            main_list = glob.glob(os.path.join(pwd,'main_*'))
            main_list = [os.path.basename(a) for a in main_list]
            main_num_list = [int(a.replace('main_','')) for a in main_list]
            main_sorted_pairs = sorted(zip(main_num_list, main_list), key=lambda pair: pair[0])

            main_list = [t[1] for t in main_sorted_pairs]
            last_main = main_list[-1]
            gen_list = glob.glob(os.path.join(pwd,last_main,'gen_*'))
            gen_list = [os.path.basename(a) for a in gen_list]
            gen_num_list = [int(a.replace('gen_','')) for a in gen_list]
            gen_sorted_pairs = sorted(zip(gen_num_list, gen_list), key=lambda pair: pair[0])
            gen_list = [t[1] for t in gen_sorted_pairs]
            last_gen = gen_list[-1]

            out = os.path.join(pwd,xyz_name)
            if os.path.exists(out):
                os.remove(out)
            train_cfg = os.path.join(pwd,last_main,last_gen,'train_mlp','train.cfg')
            cfg2xyz(ele, ele_model, train_cfg, out)
            temp = []
            for a in list(iread(out)):
                a.info['label'] = all_label
                temp.append(a)
            write(out,temp,format='extxyz')
            length = len(list(iread(out)))

            with open(os.path.join(pwd,'app.log'),'a') as f:
                f.write(f'The numbers of all structures is {length}')

if __name__ == '__main__':
    main()
