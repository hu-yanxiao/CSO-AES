import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from ase.io import iread
import time
from sklearn.metrics import r2_score
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import MinMaxScaler
from joblib import Parallel, delayed
import sys
import os
import matplotlib

#matplotlib.rcParams['font.family'] = 'Times New Roman'

def select_by_absolute_value(num1, num2):
    if abs(num1) > abs(num2):
        return num1
    else:
        return num2

def process_batch(x_batch, y_batch):
    xy_batch = np.vstack([x_batch, y_batch])
    kde_batch = gaussian_kde(xy_batch)
    z_batch = kde_batch(xy_batch)
    return z_batch

def process_with_joblib(x, y, batch_size=10000, n_jobs=-1):
    # 将数据分成多个批次
    num_batches = (len(x) + batch_size - 1) // batch_size
    batches = [(x[i * batch_size:(i + 1) * batch_size], y[i * batch_size:(i + 1) * batch_size]) for i in range(num_batches)]

    # 使用 Joblib 进行并行处理
    results = Parallel(n_jobs=n_jobs)(delayed(process_batch)(x_batch, y_batch) for x_batch, y_batch in batches)

    # 合并结果
    return np.concatenate(results)

def compute_forces(atoms_dft, atoms_sus):
    """
    计算单个 DFT 和 SUS 原子对象的力的范数和向量。
    """
    f_dft = np.linalg.norm(atoms_dft.get_forces(), axis=1).tolist()
    f_sus = np.linalg.norm(atoms_sus.get_forces(), axis=1).tolist()
    f_dft_vec = atoms_dft.get_forces().tolist()
    f_sus_vec = atoms_sus.get_forces().tolist()
    return f_dft, f_sus, f_dft_vec, f_sus_vec

def compute_energy(atoms_dft, atoms_sus):
    """
    计算单个 DFT 和 SUS 原子对象的平均能量。
    """
    e_dft = atoms_dft.get_potential_energy() / atoms_dft.get_global_number_of_atoms()
    e_sus = atoms_sus.get_potential_energy() / atoms_sus.get_global_number_of_atoms()
    return e_dft, e_sus

def compute_stress_Gpa(atoms_dft, atoms_sus):
    """
    计算单个 DFT 和 SUS 原子对象的平均能量。
    """
    ev_A32GPa = 160.21766208

    dft_s = (atoms_dft.info['virial'] / atoms_dft.get_volume() * ev_A32GPa).reshape(3, 3) #atoms_dft.info['stress'] 这里是kbar
    sus_s = (atoms_sus.info['virial'] / atoms_dft.get_volume() * ev_A32GPa).reshape(3, 3) #atoms_dft.info['stress'] ase计算的与vasp符号相反
    # diff = dft_s - sus_s
    #
    # # 计算Frobenius范数
    # frobenius_norm = np.linalg.norm(diff, 'fro')
    #
    # # 计算MAE（Frobenius范数除以9）
    # mae = frobenius_norm / 9

    return dft_s,sus_s

def dft_sus_e(dft_data, sus_data, num_processes=None,save_txt=False):
    # 读取数据
    start = time.time()
    dft = list(iread(dft_data))
    sus = list(iread(sus_data))
    end = time.time()

    dft[0].get_volume()

    print(f'read_xyz_time: {(end - start) / 60} min')

    # 并行计算
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_energy, zip(dft, sus))

    # 合并结果
    e_dft, e_sus = zip(*results)  # 解包结果
    e_dft = np.array(e_dft)
    e_sus = np.array(e_sus)

    t_save = time.time()
    if save_txt:
        os.makedirs('txt_files', exist_ok=True)
        data = np.column_stack((e_dft, e_sus))
        np.savetxt('txt_files/energy.txt', data,
                   delimiter=',', fmt='%.6f',
                   header='dft(eV/atom),sus2_pre(eV/atom)', comments='')
        print(f'Text file saved as: txt_files/force.txt')

    if save_txt:
        print(f'Data save time: {time.time() - t_save:.2f}s')

    # 绘图
    plt.scatter(e_dft, e_sus)
    # data = np.column_stack((e_dft, e_sus))
    # np.savetxt('energy.txt', data, delimiter=',', fmt='%f', header='dft,pre', comments='')

    x_lim = plt.gca().get_xlim()
    y_lim = plt.gca().get_ylim()
    line_x = min(x_lim[0], y_lim[0])
    line_y = max(x_lim[1], y_lim[1])
    plt.plot([line_x, line_y], [line_x, line_y], color='red', linestyle='--', label='y = x')

    # 计算 MAE、RMSE 和 R²
    mae = np.mean(np.abs(e_dft - e_sus)) * 1000
    rmse = np.sqrt(np.mean((e_dft - e_sus) ** 2)) * 1000
    r2 = r2_score(e_dft, e_sus)

    print(f'energy_MAE: {mae} meV, energy_RMSE: {rmse} meV, R2: {r2}')
    plt.text(0.05, 0.95, f'MAE: {mae:.3f} (meV/atom)\nRMSE: {rmse:.3f} (meV/atom)\nR\u00B2: {r2:.3f}', fontsize=12, ha='left', va='top',
             transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

    plt.ylabel('SUS² Energy (eV/atom)', fontsize=14, fontweight='bold')
    plt.xlabel('DFT Energy (eV/atom)', fontsize=14, fontweight='bold')
    plt.savefig('energy.png', dpi=300)
    plt.close()

def dft_sus_f(dft_data, sus_data, num_processes=None):
    # 读取数据
    start = time.time()
    dft = list(iread(dft_data))
    sus = list(iread(sus_data))
    end = time.time()

    print(f'read_xyz_time: {(end - start) / 60} min')

    # 并行计算
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_forces, zip(dft, sus))

    # 合并结果
    f_dft = []
    f_sus = []
    f_dft_vec = []
    f_sus_vec = []
    for result in results:
        f_dft.extend(result[0])
        f_sus.extend(result[1])
        f_dft_vec.extend(result[2])
        f_sus_vec.extend(result[3])

    # 转换为 NumPy 数组
    f_dft = np.array(f_dft)
    f_sus = np.array(f_sus)
    f_dft_vec = np.array(f_dft_vec)
    f_sus_vec = np.array(f_sus_vec)


    print(f'dft_20_num:{np.sum(f_dft > 20)}, sus_20_num:{np.sum(f_sus > 20)}, total_atomic_num:{len(f_dft)}')


    # # 保存到文件
    # data = np.column_stack((f_dft, f_sus))
    # np.savetxt('force.txt', data, delimiter=',', fmt='%f', header='dft,pre', comments='')

    # 绘图
    plt.scatter(f_dft, f_sus)
    x_lim = plt.gca().get_xlim()
    y_lim = plt.gca().get_ylim()
    line_x = min(x_lim[0], y_lim[0])
    line_y = max(x_lim[1], y_lim[1])
    plt.plot([line_x, line_y], [line_x, line_y], color='red', linestyle='--', label='y = x')

    # 计算 MAE 和 RMSE
    mae = np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) / 3) * 1000
    rmse = np.sqrt(np.mean((np.linalg.norm((f_dft_vec - f_sus_vec), axis=1)) ** 2) / 3) * 1000
    r2 = r2_score(f_dft_vec, f_sus_vec)

    print(f'force_MAE: {mae} (meV/Å), force_RMSE: {rmse} (meV/Å), R2: {r2}')
    plt.text(0.05, 0.95, f'MAE: {mae:.3f} (meV/Å)\nRMSE: {rmse:.3f} (meV/Å)\nR\u00B2: {r2:.3f}', fontsize=12, ha='left', va='top',
             transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

    plt.ylabel('SUS² Force (eV/Å)', fontsize=14, fontweight='bold')
    plt.xlabel('DFT Force (eV/Å)', fontsize=14, fontweight='bold')
    plt.savefig('force.png', dpi=300)
    plt.close()

    return f_dft, f_sus


def dft_sus_stress(dft_data, sus_data, num_processes=None,save_txt=False):
    # 读取数据
    start = time.time()
    dft = list(iread(dft_data))
    sus = list(iread(sus_data))
    end = time.time()

    print(f'read_xyz_time: {(end - start) / 60} min')

    # 并行计算
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_stress_Gpa, zip(dft, sus))
    # stress_dft = []
    # stress_sus = []
    # for result in results:
    #     stress_dft.extend(result[0])
    #     stress_sus.extend(result[1])

        # 收集所有应力数据
    dft_list = [result[0] for result in results]
    sus_list = [result[1] for result in results]

    N = len(dft_list)
    total_frobenius_mae = 0.0
    total_frobenius_squared = 0.0

    # 计算真实应力的均值矩阵
    mean_true = np.mean(dft_list, axis=0)

    # 初始化R²计算所需的平方和
    SS_res = 0.0  # 残差平方和
    SS_tot = 0.0  # 总平方和

    for dft, sus in zip(dft_list, sus_list):
        diff = sus - dft  # 预测值 - 真实值

        # MAE计算（基于Frobenius范数/9）
        frob_norm = np.linalg.norm(diff, 'fro')
        total_frobenius_mae += frob_norm / 9

        # RMSE计算准备（基于Frobenius范数平方）
        total_frobenius_squared += (frob_norm ** 2) / 9

        # R²计算准备
        SS_res += frob_norm ** 2
        SS_tot += np.linalg.norm(dft - mean_true, 'fro') ** 2

    # 计算最终指标
    MAE = total_frobenius_mae / N
    RMSE = np.sqrt(total_frobenius_squared / N)
    R2 = 1 - (SS_res / SS_tot) if SS_tot != 0 else float('nan')
    print(f'stress_MAE:{MAE} (GPa), stress_RMSE:{RMSE} (GPa), R2:{R2}')

    def plot_stress_components(dft_list, sus_list, mae, rmse, r2, save_txt):
        """
        将所有应力分量的散点图合并到一个图中

        参数:
            dft_list: DFT应力矩阵列表
            sus_list: 预测应力矩阵列表
            mae: 整体MAE
            rmse: 整体RMSE
            r2: 整体R²
            save_txt: 是否保存数据到文本文件
        """
        # 提取6个独立分量 (xx, yy, zz, xy, xz, yz)
        dft_components = []
        sus_components = []

        component_names = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']

        for dft, sus in zip(dft_list, sus_list):
            # 提取独立分量
            dft_flat = [dft[0, 0], dft[1, 1], dft[2, 2], dft[0, 1], dft[0, 2], dft[1, 2]]
            sus_flat = [sus[0, 0], sus[1, 1], sus[2, 2], sus[0, 1], sus[0, 2], sus[1, 2]]

            dft_components.append(dft_flat)
            sus_components.append(sus_flat)

        dft_components = np.array(dft_components)
        sus_components = np.array(sus_components)

        # 保存数据到文本文件（可选）
        if save_txt:
            os.makedirs('txt_files', exist_ok=True)
            all_data = np.column_stack((dft_components, sus_components))
            np.savetxt(f'txt_files/stress_all.txt', all_data,
                       delimiter=',', fmt='%.6f',
                       header='dft_xx(GPa),dft_yy(GPa),dft_zz(GPa),dft_xy(GPa),dft_xz(GPa),dft_yz(GPa),'
                              'sus2_pre_xx(GPa),sus2_pre_yy(GPa),sus2_pre_zz(GPa),sus2_pre_xy(GPa),sus2_pre_xz(GPa),sus2_pre_yz(GPa)', comments='')
            print(f'Text file saved as: txt_files/stress_all.txt')


        # 创建一个大图
        plt.figure(figsize=(10, 8))

        # 定义不同的颜色和标记
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
        markers = ['o', 's', '^', 'D', 'v', '<']

        # 绘制所有分量的散点图
        for i, (comp_name, color, marker) in enumerate(zip(component_names, colors, markers)):
            dft_comp = dft_components[:, i]
            sus_comp = sus_components[:, i]

            # 绘制散点图
            plt.scatter(dft_comp, sus_comp, s=40,
                        c=color, marker=marker, label=comp_name)

        # 添加理想线
        all_dft = dft_components.flatten()
        all_sus = sus_components.flatten()
        min_val = min(np.min(all_dft), np.min(all_sus))
        max_val = max(np.max(all_dft), np.max(all_sus))
        margin = 0.1 * (max_val - min_val)

        plt.plot([min_val - margin, max_val + margin],
                 [min_val - margin, max_val + margin], 'k--', lw=2)

        # 设置图表属性
        plt.xlabel('DFT Stress (GPa)',fontsize=14, fontweight='bold')
        plt.ylabel('Predicted Stress (GPa)',fontsize=14, fontweight='bold')
        plt.legend()

        plt.text(0.1, 0.95, f'MAE: {mae:.3f} (GPa)\nRMSE: {rmse:.3f} (GPa)\nR\u00B2: {r2:.3f}', fontsize=12,
                 ha='left',
                 va='top', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

        plt.tight_layout()

        plt.savefig('stress_components_comparison.png', dpi=300, bbox_inches='tight')
        #plt.show()
    plot_stress_components(dft_list, sus_list, MAE, RMSE, R2, save_txt)

def dft_sus_f_three(dft_data, sus_data, num_processes=None):
    # 读取数据
    start = time.time()
    dft = list(iread(dft_data))
    sus = list(iread(sus_data))
    end = time.time()

    print(f'read_xyz_time: {(end-start)/60} min')

    # 并行计算
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_forces, zip(dft, sus))

    # 合并结果
    f_dft = []
    f_sus = []
    f_dft_vec = []
    f_sus_vec = []
    for result in results:
        f_dft.extend(result[0])
        f_sus.extend(result[1])
        f_dft_vec.extend(result[2])
        f_sus_vec.extend(result[3])
    #print(f_dft,f_dft_vec)
    # for result in results:
    #     print(result)

    # 转换为分量形式
    three_f_dft_vec = [list(t) for t in zip(*f_dft_vec)]
    three_f_sus_vec = [list(t) for t in zip(*f_sus_vec)]

    # 绘图
    names = ['Fx', 'Fy', 'Fz']
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    for ax, x, y, name in zip(axs, three_f_dft_vec, three_f_sus_vec, names):
        x = np.array(x)
        y = np.array(y)
        ax.scatter(x, y)
        ax.set_title(f'{name}')
        ax.set_ylabel('SUS² Force (eV/Å)', fontsize=14, fontweight='bold')
        ax.set_xlabel('DFT Force (eV/Å)', fontsize=14, fontweight='bold')
        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()
        line_x = min(x_lim[0], y_lim[0])
        line_y = max(x_lim[1], y_lim[1])
        ax.plot([line_x, line_y], [line_x, line_y], color='red', linestyle='--', label='y = x')

        mae = np.mean(np.abs(x - y)) * 1000
        rmse = np.sqrt(np.mean((x - y) ** 2)) * 1000

        r2 = r2_score(x, y)

        print(f'{name}: Force_MAE: {mae} meV, Force_RMSE: {rmse} meV')
        ax.text(0.05, 0.95, f'MAE: {mae:.3f} (meV/Å)\nRMSE: {rmse:.3f} (meV/Å)\nR\u00B2: {r2:.3f}', fontsize=12, ha='left',
                va='top', transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))

    plt.tight_layout()
    plt.savefig('force_three.png', dpi=300)
    plt.close()

def plot_f(dft_data,sus_data,num_processes=None):
    # 设置全局字体和字号
    # plt.rcParams['font.family'] = 'serif'
    # plt.rcParams['font.serif'] = 'Times New Roman'
    # plt.rcParams['font.size'] = 20  # 全局字体大小设置为20

    # 读取数据
    start = time.time()
    dft = list(iread(dft_data))
    sus = list(iread(sus_data))
    end = time.time()

    print(f'read_xyz_time: {(end - start) / 60} min')

    time_1 = time.time()
    # 并行计算
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_forces, zip(dft, sus))

    # 合并结果
    f_dft, f_sus, f_dft_vec, f_sus_vec = [], [], [], []
    for result in results:
        f_dft.extend(result[0])
        f_sus.extend(result[1])
        f_dft_vec.extend(result[2])
        f_sus_vec.extend(result[3])

    # 转换为 NumPy 数组
    f_dft = np.array(f_dft)
    f_sus = np.array(f_sus)
    f_dft_vec = np.array(f_dft_vec)
    f_sus_vec = np.array(f_sus_vec)


    x = f_dft
    y = f_sus

    mae = np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) / 3) * 1000
    rmse = np.sqrt(np.mean((np.linalg.norm((f_dft_vec - f_sus_vec), axis=1)) ** 2) / 3) * 1000
    r2 = r2_score(f_dft_vec, f_sus_vec)

    print(f'force_MAE:{mae} (meV/Å), force_RMSE:{rmse} (meV/Å), R2:{r2}')

    x_sampled, y_sampled = x, y

    time_2 = time.time()
    print(time_2-time_1)

    # 计算密度
    z = process_with_joblib(x_sampled,y_sampled)
    # xy = np.vstack([x_sampled, y_sampled])
    # z = gaussian_kde(xy)(xy)

    # 按密度排序，确保密度高的点在上层显示
    idx = z.argsort()
    x_sorted, y_sorted, z_sorted = x_sampled[idx], y_sampled[idx], z[idx]

    # 将颜色映射的值归一化到[0, 1]范围
    z_normalized = (z_sorted - z_sorted.min()) / (z_sorted.max() - z_sorted.min())

    time_3 = time.time()
    print(time_3 - time_2)

    # 绘图
    fig, ax = plt.subplots()


    # 绘制散点图，颜色表示密度
    scatter = ax.scatter(
        x_sorted,
        y_sorted,
        marker='o',
        c=z_normalized,
        edgecolors='none',
        s=10,
        cmap='Spectral_r'
    )

    # 添加颜色条
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(scatter, cax=cax, label='density')

    # 设置颜色条标签字号和刻度字号
    cbar.ax.set_ylabel('density', fontsize=20)
    cbar.ax.tick_params()
    #cbar.ax.tick_params(labelsize=20)

    # 设置坐标轴标签

    ax.set_ylabel('SUS\u00B2 Force (eV/Å)', fontsize=14, fontweight='bold')
    ax.set_xlabel('DFT Force (eV/Å)', fontsize=14, fontweight='bold')

    # 调整坐标轴刻度标签的字号
    #ax.tick_params(axis='both', which='major', labelsize=20)

    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()
    print('x_axis_range:',x_lim)
    print('y_axis_range:',y_lim)
    line_x = select_by_absolute_value(x_lim[0], y_lim[0])
    line_y = select_by_absolute_value(x_lim[1], y_lim[1])
    print(line_x,line_y)
    ax.plot([line_x,line_y], [line_x,line_y], color='red', linestyle='--', label='y = x')

    ax.text(0.05, 0.95, f'MAE: {mae:.3f} (meV/Å)\nRMSE: {rmse:.3f} (meV/Å)\nR\u00B2: {r2:.3f}', fontsize=12, ha='left', va='top',
             transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))

    # 调整布局以防止元素重叠
    plt.tight_layout()

    # 显示图形
    plt.savefig('density', dpi=300)
    #plt.show()
    plt.close()


def opt_compute_forces(atoms_dft, atoms_sus):
    """Optimized force computation that returns flattened force arrays"""
    return (atoms_dft.get_forces().ravel(),
            atoms_sus.get_forces().ravel())

def dft_sus_f_combined(dft_data, sus_data, num_processes=None, save_txt=False, save_npy=False):
    # Read data
    t0 = time.time()
    dft = list(iread(dft_data))
    sus = list(iread(sus_data))
    print(f'Data reading time: {time.time() - t0:.2f}s')

    # Parallel force computation
    t1 = time.time()
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(opt_compute_forces, zip(dft, sus))

    # Stack all forces into single arrays (flattened components)
    forces_dft = np.concatenate([r[0] for r in results])
    forces_sus = np.concatenate([r[1] for r in results])
    print(f'Force computation time: {time.time() - t1:.2f}s')

    # Save data if requested
    t_save = time.time()
    if save_npy:
        os.makedirs('npy_files', exist_ok=True)
        np.save('npy_files/forces_dft.npy', forces_dft)
        np.save('npy_files/forces_sus.npy', forces_sus)
        print(f'Numpy files saved in: npy_files/')

    if save_txt:
        os.makedirs('txt_files', exist_ok=True)
        data = np.column_stack((forces_dft, forces_sus))
        np.savetxt('txt_files/force.txt', data,
                   delimiter=',', fmt='%.6f',
                   header='dft(eV/Å),sus2_pre(eV/Å)', comments='')
        print(f'Text file saved as: txt_files/force.txt')

    if save_npy or save_txt:
        print(f'Data save time: {time.time() - t_save:.2f}s')

    # Calculate metrics
    t2 = time.time()
    diff = forces_dft - forces_sus
    mae = np.mean(np.abs(diff)) * 1000  # meV/Å
    rmse = np.sqrt(np.mean(diff ** 2)) * 1000  # meV/Å
    r2 = r2_score(forces_dft, forces_sus)

    # Plot all components together
    plt.figure(figsize=(8, 8))
    plt.scatter(forces_dft, forces_sus, alpha=0.3, s=10,
                c='blue', edgecolors='none')

    # Set equal limits and reference line
    lim_min = min(forces_dft.min(), forces_sus.min())
    lim_max = max(forces_dft.max(), forces_sus.max())
    plt.plot([lim_min, lim_max], [lim_min, lim_max],
             'r--', lw=1, alpha=0.7)

    plt.xlabel('DFT Force (eV/Å)', fontsize=16,fontweight='bold')
    plt.ylabel('SUS\u00B2 Force (eV/Å)', fontsize=16,fontweight='bold')
    #plt.title('Force Components Comparison', fontsize=16)

    # Annotate with metrics
    metrics_text = (f'MAE = {mae:.2f} meV/Å\n'
                    f'RMSE = {rmse:.2f} meV/Å\n'
                    f'R² = {r2:.3f}')
    plt.text(0.05, 0.95, metrics_text,
             transform=plt.gca().transAxes,
             ha='left', va='top',
             bbox=dict(facecolor='white', alpha=0.8))

    #plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('force_combined.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f'Metrics and plotting time: {time.time() - t2:.2f}s')
    print(f'Total execution time: {time.time() - t0:.2f}s')
    print(f'\nFinal Metrics (All force components):')
    print(f'MAE: {mae:.2f} meV/Å')
    print(f'RMSE: {rmse:.2f} meV/Å')
    print(f'R²: {r2:.3f}')


# Usage examples:
# dft_sus_f_combined('dft.xyz', 'sus.xyz', num_processes=8)  # No saving
# dft_sus_f_combined('dft.xyz', 'sus.xyz', num_processes=8, save_txt=True)  # Save text only
# dft_sus_f_combined('dft.xyz', 'sus.xyz', num_processes=8, save_npy=True)  # Save numpy only
# dft_sus_f_combined('dft.xyz', 'sus.xyz', num_processes=8, save_txt=True, save_npy=True)  # Save both

if __name__ == "__main__":
    dft_data = sys.argv[1]
    sus_data = sys.argv[2]
    print(f'dft_data:{dft_data} sus_data:{sus_data}')

    num_processes = 64

    start = time.time()
    dft_sus_e(dft_data, sus_data, num_processes=num_processes,save_txt=True)
    end = time.time()
    print(f'run_time: {(end - start) / 60} min')
    #
    start = time.time()
    dft_sus_f(dft_data, sus_data, num_processes=num_processes)
    end = time.time()
    print(f'run_time: {(end - start) / 60} min')


    start = time.time()
    dft_sus_f_three(dft_data, sus_data, num_processes=num_processes)
    end = time.time()
    print(f'run_time: {(end - start) / 60} min')

    start = time.time()
    dft_sus_f_combined(dft_data, sus_data, num_processes=num_processes,save_txt=True)
    end = time.time()
    print(f'run_time: {(end - start) / 60} min')

    start = time.time()
    dft_sus_stress(dft_data, sus_data, num_processes=num_processes, save_txt=True)
    end = time.time()
    print(f'run_time: {(end - start) / 60} min')

    #
    # start = time.time()
    # plot_f(dft_data, sus_data, num_processes=num_processes)
    # end = time.time()
    # print(f'run_time: {(end - start) / 60} min')
