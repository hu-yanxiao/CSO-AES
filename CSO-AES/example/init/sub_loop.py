import os
import yaml
from path import Path
import shutil
from tqdm import tqdm
from cso_aes_das.mkdir import mkdir_vasp
from cso_aes_das.work_dir import work_deepest_dir, submit_lammps_task, scf_dir, bsub_dir, check_finish,check_scf, check_filter_xyz_0,deepest_dir,delete_dump
from cso_aes_das.calc_ensemble_ambiguity import get_force_ambiguity, ambiguity_extract,check_lmp_error
from cso_aes_das.logger import setup_logger
from cso_aes_das.sample_xyz import sample_main
from cso_aes_das.scf_lmp_data import scf_lammps_data
from ase.io import write,iread
from cso_aes_das.main_calc import main_calc,check_and_modify_calc_dir

from cso_aes_das.other import remove,end_yaml,touch,mkdir
from cso_aes_das.train_mlp import pre_train_mlp,start_train
from cso_aes_das.das_update_ambiguity import das_update_ambiguity,record_yaml,af_limit_update,af_limit_record
import glob
import json
import sys
from cso_aes_encode.mlp_encode_sample_flow import main_sample_flow
import time
logger = setup_logger()

def main():
    pwd = os.getcwd()

    ######导入计算引擎提取函数#####
    path = os.path.dirname(os.path.dirname(pwd))
    sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join(path, 'init')))
    import bsub_script
    scf_cal_engine = bsub_script.scf_cal_engine
    if bsub_script.scf_cal_engine == 'abacus':
        from cso_aes_das.abacus_main_xyz import abacus_main_xyz as scf2xyz
    elif bsub_script.scf_cal_engine == 'cp2k':
        from cso_aes_das.cp2k_main_xyz import cp2k_main_xyz as scf2xyz
    elif bsub_script.scf_cal_engine == 'vasp':
        from sus_das.vasp_main_xyz import vasp_main_xyz as scf2xyz
    else:
        raise ValueError(f'{lmp_in.scf_cal_engine} no exist! but you can add new_scf_cal_engine by yourself.')

    try:
        end_threshold_low, end_threshold_high, end_n, end_cluster_threshold_init, end_k = end_yaml('end.yaml')
    except:
        end_threshold_low, end_threshold_high, end_n, end_cluster_threshold_init, end_k = None, None, None, None,None

    with open('parameter.yaml', 'rb') as file:
        yaml_data = yaml.safe_load(file)

    mlp_MD = yaml_data['mlp_MD']
    ele = yaml_data['ele']
    size = tuple(eval(yaml_data['size']))
    task_submission_method = yaml_data['task_submission_method']
    mlp_nums = yaml_data['mlp_nums']
    ele_model = yaml_data['ele_model']
    original_COMMAND = bsub_script.original_COMMAND
    subsequent_COMMAND = bsub_script.subsequent_COMMAND

    nvt_lattice_scaling_factor = yaml_data['nvt_lattice_scaling_factor']

    das_ambiguity = yaml_data['das_ambiguity']
    af_default = yaml_data['af_default']
    #af_limit = yaml_data['af_limit']   #后面会根据 af_limit_list更新af_limit的值 update注意不要运行两遍,否则可能增加两个True
    af_failed = yaml_data['af_failed']
    over_fitting_factor = yaml_data['over_fitting_factor']

    threshold_low = yaml_data['threshold_low']
    threshold_high = yaml_data['threshold_high']
    end = yaml_data['end']
    num_elements = yaml_data['num_elements']

    sample = yaml_data['sample']
    n = sample['n']
    cluster_threshold_init = sample['cluster_threshold_init']
    k = sample['k']
    clustering_by_ambiguity = sample['clustering_by_ambiguity']

    mlp_encode_model = yaml_data['mlp_encode_model']
    bw_method = yaml_data['bw_method']  # scott, Freedman_Diaconis, std(std/10), self_input
    bw = yaml_data['bw']
    body_list = yaml_data['body_list']  # 同时考虑的多体编码
    mtp_type = yaml_data['mtp_type']
    stru_num = yaml_data['stru_num']
    coverage_rate_threshold = yaml_data['coverage_rate_threshold']
    coverage_rate_method = yaml_data['coverage_rate_method']

    dft = yaml_data['dft']
    calc_dir_num = dft['calc_dir_num']
    force_threshold = dft['force_threshold']

    #强行置为1
    if mlp_encode_model:
        mlp_nums = 1
    else:
        if mlp_nums<3:
            print("During the use of DAS, mlp_nums should be greater than or equal to 3")
            logger.info("During the use of DAS, mlp_nums should be greater than or equal to 3")
            sys.exit(1)
    ###########训练势函数#############
    label = pre_train_mlp(pwd, mlp_nums,ele,ele_model,logger,original_COMMAND,subsequent_COMMAND)
    if label:
        start_train(pwd,task_submission_method,mlp_nums,logger)

    ######检查训练的势函数是否报错######
    gen_num = int(os.path.basename(pwd).replace('gen_', ''))
    train_mlp = deepest_dir(pwd, 'train_mlp')
    if gen_num != 0 or (gen_num == 0 and len(train_mlp) != 1):
        for a in deepest_dir(pwd, 'train_mlp'):
            log_file = os.path.join(a, 'logout')
            with open(log_file, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if 'Killed' in line:
                        print(f"Scaling_mlp training error, please check! : {log_file}")
                        logger.info(f'Scaling_mlp training error, please check! : {log_file}')
                        touch(pwd, '__error__')
                        sys.exit(1)
                end_lines = lines[-25:]
                for a in end_lines:
                    if 'nan' in a:
                        print(f"Scaling_mlp training error, please check! : {log_file}")
                        logger.info(f'Scaling_mlp training error, please check! : {log_file}')
                        touch(pwd, '__error__')
                        sys.exit(1)

    # ##########建立目录和文件###########
    mkdir_vasp(pwd, mlp_MD, ele, size, mlp_nums, ele_model,
               nvt_lattice_scaling_factor,  mlp_encode_model)

    dirs_1 = work_deepest_dir(pwd)
    dirs_2 = bsub_dir(pwd)

    ###开始遍历计算目录提交任务########
    submit_lammps_task(pwd,logger,task_submission_method)

    ############检查任务是否全部计算完成######
    check_finish(dirs_1,logger,'All MD calculations have been completed')

    if mlp_encode_model == False:
        # ##########计算分歧，挑出分歧大的结构################
        last_gen = 'gen_' + str((int(gen_num) - 1))
        last_gen_path = os.path.join(os.path.dirname(pwd), last_gen)
        mtp_path = os.path.join(pwd, 'mtp')
        last_scf_filter_xyz = os.path.join(last_gen_path, 'scf_lammps_data', 'scf_filter.xyz')

        if end_threshold_low != threshold_low or end_threshold_high != threshold_high:
            for dir in dirs_2:
                os.chdir(dir)
                files_to_delete = glob.glob('*filter*')
                for file_name in files_to_delete:
                    remove(file_name)

            af_adaptive = None

            if das_ambiguity:
                if gen_num == 0:
                    xyz = None
                    model_fns = None
                else:
                    xyz = last_scf_filter_xyz
                    model_fns = [os.path.join(mtp_path,a) for a in glob.glob(os.path.join(mtp_path,'current*'))]
                yaml_file = os.path.join(pwd, 'parameter.yaml')
                new_af_limit = af_limit_update(pwd, yaml_file)
                af_adaptive,label = das_update_ambiguity(ele=ele, ele_model=ele_model, af_default=af_default, af_limit=new_af_limit,
                                     af_failed=af_failed, over_fitting_factor=over_fitting_factor, logger=logger).run(xyz,
                                                                                                                      model_fns,
                                                                                                                      gen_num)
                af_limit_record(pwd, label)

                threshold_low = float(af_adaptive)
                threshold_high = af_failed

            select_stru_num = 0

            for a in tqdm(dirs_1):
                dir = os.path.join(pwd, a)
                total_stru = get_force_ambiguity(dir)  # 计算分歧阈值
                num, stru,interval,hist= ambiguity_extract(dir, 'force.0.dump', 'af.out', threshold_low,threshold_high, ele, ele_model, end, num_elements)  # 挑出分歧阈值大的结构
                path_parts = os.path.normpath(dir).split(os.sep)
                final_path = os.sep.join(path_parts[-3:])
                filter_path = os.path.join((os.sep.join(path_parts[:-2])), 'filter.xyz')
                write(filter_path, stru, format='extxyz', append=True)  # 结构写入xyz
                logger.info(f'{final_path}: According to the ambiguity in the {round(threshold_low, 3)}-{threshold_high} range , {num} structures are selected from {total_stru} structures. Interval:{interval} Statistical number:{hist}')
                select_stru_num += num
                error, message = check_lmp_error(dir)
                if error:
                    logger.warning(f'{message}')
            yaml_file = os.path.join(pwd, 'parameter.yaml')
            record_yaml(yaml_file, float(af_adaptive), int(select_stru_num))
        else:
            logger.info(f'end_threshold equals threshold: skip to select the structure by ambiguity')

        if end_threshold_low == threshold_low and end_threshold_high == threshold_high and end_n == n and end_cluster_threshold_init == cluster_threshold_init and end_k == k:
            logger.info(f'(threshold, n, cluster_threshold_init, k) parameters are equal: skip to select the structure by MBTR+Brich')
        else:
            for dir in dirs_2:
                os.chdir(dir)
                if os.path.getsize('filter.xyz') != 0:
                    num = len(list(iread('filter.xyz')))
                    if n*k <= num:
                        select, total = sample_main(os.getcwd(), n=n, threshold_init=cluster_threshold_init, k=k, clustering_by_ambiguity=clustering_by_ambiguity)  # 采样结构生成pkl等文件
                        name = os.path.basename(dir)
                        logger.info(f'{name}: selected {select} structures from {total} structures in data by MBTR+Brich.')
                    elif n*k > num:
                        shutil.copy('filter.xyz',f'{num}_sample_filter.xyz')
    elif mlp_encode_model == True:
        sample_xyz_list = glob.glob(os.path.join(pwd,'work','*_sample_filter.xyz'))
        ll = len(sample_xyz_list)
        if ll == 0:
            main_sample_flow(pwd, dirs_1, bw, bw_method, body_list, ele, ele_model, mtp_type,
                                                    stru_num, coverage_rate_threshold, coverage_rate_method, logger)
        elif ll == 1:
            logger.info(
                f'*_sample_filter.xyz already exists.({sample_xyz_list[0]})')
        else:
            raise ValueError('Multiple *_sample_filter.xyz, Please delete!')
    else:
        raise ValueError('mlp_encode_model must be a bool value')

    # #############yaml保存##################
    init = os.path.join(pwd, 'parameter.yaml')
    end = os.path.join(pwd, 'end.yaml')
    shutil.copy(init, end)

    if mlp_encode_model == False:
        if check_filter_xyz_0(pwd):
            #########把文件汇总至scf_lammps_data目录，开始批量DFT计算############
            if end_threshold_low == threshold_low and end_threshold_high == threshold_high and end_n == n and end_cluster_threshold_init == cluster_threshold_init and end_k == k and calc_dir_num:
                logger.info('(calc_dir_num, threshold, n, cluster_threshold_init, k) parameters are equal: skip scf calculations')
            else:
                delete_filter_select = 'no'  # 把挑选出来剩余的xyz删除
                delete_scf_file = 'yes'  # 把之前的scf下所有文件全部删除，(如果calc_dir_num改变，对delete_scf_file = 'no'会影响)

                total_atom_list = scf_lammps_data(delete_filter_select, pwd)
                scf_path = os.path.join(pwd, 'scf_lammps_data', 'scf')
                os.chdir(scf_path)

                if delete_scf_file == 'yes':
                    for i in [a for a in os.listdir() if a != 'total_sample_filter.xyz']:
                        remove(i)
                calc_dir_num = check_and_modify_calc_dir(pwd, total_atom_list, calc_dir_num)
                num = main_calc(total_atom_list, calc_dir_num, pwd)
                if not check_scf(pwd):
                    os.system('python start_calc.py')
                logger.info(f'The {num} structures are divided into {calc_dir_num} dft calculation tasks to be submitted')

            #####检查SCF是否完成#####
            logger.info(f'In the process of checking whether the SCF calculation is complete......')
            check_finish(scf_dir(pwd), logger, 'All scf calculations have been completed')

            #####提取已经计算好的且SCF收敛的结构,不收敛的生成busb.lsf文件#####
            current = os.path.join(pwd,'scf_lammps_data','scf','filter')
            out_name = os.path.join(pwd,'scf_lammps_data','scf_filter.xyz')
            ori_out_name = os.path.join(pwd, 'scf_lammps_data', 'ori_scf_filter.xyz')
            remove(out_name)  #这里默认追加写入，删掉是为了防止，重复写入
            ok_count,len_count,no_success_path, force_count,_ = scf2xyz(current, out_name, ori_out_name, force_threshold)

            with open(os.path.join(pwd,'no_success_path.json'), 'w') as json_file:
                json.dump(no_success_path, json_file)

            gen = os.path.basename(pwd)
            logger.info(f'Active learning continues: {gen} | {scf_cal_engine}_completed_number:{ok_count} | '
                        f'Successful_collection_structure/scf_convergent_number:{len_count} | '
                        f'force_threshold_number({force_threshold}):{force_count}')
        else:
            logger.info(f'The active learning loop ends')
            touch(pwd, '__end__')
    elif mlp_encode_model == True:
        sample_xyz = glob.glob(os.path.join(pwd,'work','*_sample_filter.xyz'))[0]
        if os.path.getsize(sample_xyz) != 0:
            total_atom_list = list(iread(sample_xyz))
            scf_lammps_data_path = os.path.join(pwd, 'scf_lammps_data')
            scf_path = os.path.join(scf_lammps_data_path, 'scf')
            mkdir(scf_lammps_data_path)
            mkdir(scf_path)
            os.chdir(scf_path)
            calc_dir_num = check_and_modify_calc_dir(pwd, total_atom_list, calc_dir_num)
            num = main_calc(total_atom_list, calc_dir_num, pwd)
            if not check_scf(pwd):
                os.system('python start_calc.py')
            logger.info(f'The {num} structures are divided into {calc_dir_num} dft calculation tasks to be submitted')

            #####检查SCF是否完成#####
            logger.info(f'In the process of checking whether the SCF calculation is complete......')
            check_finish(scf_dir(pwd), logger, 'All scf calculations have been completed')

            #####提取已经计算好的且SCF收敛的结构,不收敛的生成busb.lsf文件#####
            current = os.path.join(pwd, 'scf_lammps_data', 'scf', 'filter')
            out_name = os.path.join(pwd, 'scf_lammps_data', 'scf_filter.xyz')
            ori_out_name = os.path.join(pwd, 'scf_lammps_data', 'ori_scf_filter.xyz')
            remove(out_name)  # 这里默认追加写入，删掉是为了防止，重复写入
            ok_count,len_count,no_success_path, force_count, force_of_force_count_0 = scf2xyz(current, out_name, ori_out_name, force_threshold)

            with open(os.path.join(pwd, 'no_success_path.json'), 'w') as json_file:
                json.dump(no_success_path, json_file)

            gen = os.path.basename(pwd)

            if force_count !=0:
                logger.info(f'Active learning continues: {gen} | {scf_cal_engine}_completed_number:{ok_count} | '
                            f'Successful_collection_structure/scf_convergent_number:{len_count} | '
                            f'force_threshold_number({force_threshold}):{force_count}') #force_count：小于force_threshold的数
            else:
                logger.info(f'Active learning continues: {gen} | {scf_cal_engine}_completed_number:{ok_count} | '
                            f'Successful_collection_structure/scf_convergent_number:{len_count} | '
                            f'force_threshold_number({force_threshold}):1 (Warning force_count==0, The structure with the minimum atomic force(force:{force_of_force_count_0}) is chosen.)')  # force_count：小于force_threshold的数
        else:
            logger.info(f'The active learning loop ends')
            touch(pwd, '__end__')
    delete_dump(dirs_1)
    touch(pwd,'__ok__')

if __name__ == '__main__':
    main()
