#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare DFT and SUS² results: energy, |F| and Fx/Fy/Fz.
Scatter plots are density-coloured; optional CSV export for Origin.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from multiprocessing import Pool
from ase.io import iread
import time, sys, os
from sklearn.metrics import r2_score
from joblib import Parallel, delayed
from matplotlib import rcParams
import matplotlib.ticker as ticker

# ========= 开关区 =========
SAVE_DATA = False  # True → 额外保存 3 个 CSV
num_processes = 10  # 并行进程数
cmap_name = 'Spectral_r'
scatter_size = 10  # 每个点半径

# 设置全局字体为Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams.update({
    'font.size': 30,
    'axes.titlesize': 30,
    'axes.titleweight': 'bold',
    'axes.labelsize': 30,
    'axes.labelweight': 'bold',
    'xtick.labelsize': 30,
    'ytick.labelsize': 30,
    'legend.fontsize': 20,
    'legend.title_fontsize': 20,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'axes.linewidth': 4,
    'xtick.major.width': 4,
    'ytick.major.width': 4,
    'xtick.major.size': 14,
    'ytick.major.size': 14,
    'xtick.minor.width': 4,
    'ytick.minor.width': 4,
    'xtick.minor.size': 7,
    'ytick.minor.size': 7,
})


# ========= 辅助函数 =========
def hist2d_density(x, y, bins=120):
    """用 2D 直方图快速估算每个点的密度（对大数据更快）。"""
    h, xe, ye = np.histogram2d(x, y, bins=bins)
    ix = np.clip(np.digitize(x, xe) - 1, 0, bins - 1)
    iy = np.clip(np.digitize(y, ye) - 1, 0, bins - 1)
    return h[ix, iy]


def density_scatter(ax, x, y, bins=120):
    """在指定 ax 上绘制密度散点并返回 scatter 句柄。"""
    z = hist2d_density(x, y, bins)
    order = z.argsort()  # 低密度点先画
    sc = ax.scatter(x[order], y[order],
                    c=z[order], s=scatter_size, cmap=cmap_name,
                    norm=mcolors.LogNorm(), edgecolors='none')
    return sc


def compute_energy(a_dft, a_sus):
    e_dft = a_dft.get_potential_energy() / a_dft.get_global_number_of_atoms()
    e_sus = a_sus.get_potential_energy() / a_sus.get_global_number_of_atoms()
    return e_dft, e_sus


def compute_forces(a_dft, a_sus):
    return a_dft.get_forces(), a_sus.get_forces()


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


def select_abs(a, b):  # 保留原作者的小辅助
    return a if abs(a) > abs(b) else b


# ========= 能量 =========
def dft_sus_e(dft_xyz, sus_xyz, num_processes=None):
    t0 = time.time()
    dft = list(iread(dft_xyz))
    sus = list(iread(sus_xyz))
    print(f'read_xyz_time: {time.time() - t0:.1f}s  frames={len(dft)}')

    with Pool(processes=num_processes) as pool:
        res = pool.starmap(compute_energy, zip(dft, sus))
    e_dft, e_sus = map(np.array, zip(*res))

    # ===== 绘图 =====
    fig, ax = plt.subplots(figsize=(10, 8))  # 添加这行
    sc = density_scatter(ax, e_dft, e_sus, bins=100)
    cbar = fig.colorbar(sc, ax=ax, label='density')
    cbar.outline.set_linewidth(2)
    cbar.set_label('Density', fontsize=25, fontweight='bold')
    cbar.ax.tick_params(which='major', labelsize=20, width=2, length=10, direction='in')
    cbar.ax.tick_params(which='minor', labelsize=12, width=2, length=6, direction='in')

    lim_low = min(ax.get_xlim()[0], ax.get_ylim()[0])
    lim_high = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lim_low, lim_high], [lim_low, lim_high], '--')  # 灰阶

    mae = np.mean(np.abs(e_dft - e_sus)) * 1000
    rmse = np.sqrt(np.mean((e_dft - e_sus) ** 2)) * 1000
    r2 = r2_score(e_dft, e_sus)
    ax.text(0.05, 0.95,
            f'MAE: {mae:.3f} (meV/atom)\n'
            f'RMSE: {rmse:.3f} (meV/atom)\n'
            f'R\u00B2: {r2:.3f}',
            transform=ax.transAxes, va='top',
            bbox=dict(facecolor='white', alpha=.5), fontsize=20)

    ax.set_xlabel('DFT Energy (eV/atom)', fontweight='bold')
    ax.set_ylabel('SUS2 Energy (eV/atom)', fontweight='bold')
    plt.tight_layout()
    plt.savefig('energy.jpg', dpi=300)
    plt.close()

    if SAVE_DATA:
        np.savetxt('energy_data.csv',
                   np.column_stack((e_dft, e_sus)),
                   delimiter=',', header='DFT_E,SUS2_E', comments='')


# ========= 总力 =========
def dft_sus_f(dft_xyz, sus_xyz, num_processes=None):
    t0 = time.time()
    dft = list(iread(dft_xyz))
    sus = list(iread(sus_xyz))
    print(f'read_xyz_time: {time.time() - t0:.1f}s')

    with Pool(processes=num_processes) as pool:
        res = pool.starmap(compute_forces, zip(dft, sus))

    f_dft_vec = np.vstack([r[0] for r in res])
    f_sus_vec = np.vstack([r[1] for r in res])
    f_dft = np.linalg.norm(f_dft_vec, axis=1)
    f_sus = np.linalg.norm(f_sus_vec, axis=1)

    # ===== 绘图 =====
    fig, ax = plt.subplots(figsize=(10, 8))  # 添加这行
    sc = density_scatter(ax, f_dft, f_sus, bins=100)
    cbar = fig.colorbar(sc, ax=ax, label='density')
    cbar.outline.set_linewidth(2)
    cbar.set_label('Density', fontsize=25, fontweight='bold')
    cbar.ax.tick_params(which='major', labelsize=20, width=2, length=10, direction='in')
    cbar.ax.tick_params(which='minor', labelsize=12, width=2, length=6, direction='in')

    lim_low = min(ax.get_xlim()[0], ax.get_ylim()[0])
    lim_high = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lim_low, lim_high], [lim_low, lim_high], '--')

    mae = np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) / 3) * 1000
    rmse = np.sqrt(np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) ** 2) / 3) * 1000
    r2 = r2_score(f_dft_vec, f_sus_vec)
    ax.text(0.05, 0.95,
            f'MAE: {mae:.3f} (meV/Å)\n'
            f'RMSE: {rmse:.3f} (meV/Å)\n'
            f'R\u00B2: {r2:.3f}',
            transform=ax.transAxes, va='top',
            bbox=dict(facecolor='white', alpha=.5))

    ax.set_xlabel('DFT Force (eV/Å)', fontweight='bold')
    ax.set_ylabel('SUS2 Force (eV/Å)', fontweight='bold')
    plt.tight_layout()
    plt.savefig('force.jpg', dpi=300)
    plt.close()

    if SAVE_DATA:
        np.savetxt('force_data.csv',
                   np.column_stack((f_dft, f_sus)),
                   delimiter=',', header='|F_DFT|,|F_SUS2|', comments='')

    return f_dft_vec, f_sus_vec


# ========= Fx/Fy/Fz =========
def dft_sus_f_three(dft_xyz, sus_xyz, num_processes=None):
    t0 = time.time()
    dft = list(iread(dft_xyz))
    sus = list(iread(sus_xyz))
    print(f'read_xyz_time: {time.time() - t0:.1f}s')

    with Pool(processes=num_processes) as pool:
        res = pool.starmap(compute_forces, zip(dft, sus))

    f_dft_vec = np.vstack([r[0] for r in res])
    f_sus_vec = np.vstack([r[1] for r in res])

    names = ['Fx', 'Fy', 'Fz']
    fig, axs = plt.subplots(1, 3, figsize=(30, 8))

    for i, ax in enumerate(axs):
        x = f_dft_vec[:, i]
        y = f_sus_vec[:, i]
        sc = density_scatter(ax, x, y, bins=80)
        lim_low = min(ax.get_xlim()[0], ax.get_ylim()[0])
        lim_high = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lim_low, lim_high], [lim_low, lim_high], '--')

        mae = np.mean(np.abs(x - y)) * 1000
        rmse = np.sqrt(np.mean((x - y) ** 2)) * 1000
        r2 = r2_score(x, y)
        ax.text(0.05, 0.95,
                f'MAE: {mae:.1f}\nRMSE: {rmse:.1f}\nR²: {r2:.3f}',
                transform=ax.transAxes, va='top',
                bbox=dict(facecolor='white', alpha=.5))
        ax.set_xlabel('DFT Force (eV/Å)')
        ax.set_ylabel('SUS2 Force (eV/Å)')
        ax.set_title(names[i])

        # 只在中间子图加一次 colorbar
        if i == 1:
            cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label='density')
            cbar.outline.set_linewidth(2)
            cbar.set_label('Density', fontsize=25, fontweight='bold')
            cbar.ax.tick_params(which='major', labelsize=20, width=2, length=10, direction='in')
            cbar.ax.tick_params(which='minor', labelsize=12, width=2, length=6, direction='in')

    plt.tight_layout()
    plt.savefig('force_three.jpg', dpi=300)
    plt.close()

    if SAVE_DATA:
        np.savetxt('force_components.csv',
                   np.column_stack((f_dft_vec, f_sus_vec)),
                   delimiter=',',
                   header='Fx_DFT,Fy_DFT,Fz_DFT,Fx_SUS2,Fy_SUS2,Fz_SUS2',
                   comments='')


# ========= 应力分析 =========
def dft_sus_stress(dft_xyz, sus_xyz, num_processes=None):
    """应力分析函数"""
    t0 = time.time()
    dft = list(iread(dft_xyz))
    sus = list(iread(sus_xyz))
    print(f'read_xyz_time: {time.time() - t0:.1f}s  frames={len(dft)}')

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_stress_Gpa, zip(dft, sus))

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
    print(f'stress_MAE: {MAE:.3f} (GPa), stress_RMSE: {RMSE:.3f} (GPa), R2: {R2:.3f}')

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
    if SAVE_DATA:
        os.makedirs('txt_files', exist_ok=True)
        all_data = np.column_stack((dft_components, sus_components))
        np.savetxt('txt_files/stress_all.txt', all_data,
                   delimiter=',', fmt='%.6f',
                   header='dft_xx(GPa),dft_yy(GPa),dft_zz(GPa),dft_xy(GPa),dft_xz(GPa),dft_yz(GPa),'
                          'sus2_pre_xx(GPa),sus2_pre_yy(GPa),sus2_pre_zz(GPa),sus2_pre_xy(GPa),sus2_pre_xz(GPa),sus2_pre_yz(GPa)',
                   comments='')
        print('Text file saved as: txt_files/stress_all.txt')

    # ===== 图1: 所有分量合并的密度散点图 =====
    fig1, ax1 = plt.subplots(figsize=(10, 8))

    # 将所有分量合并为单个数组
    dft_all = dft_components.flatten()
    sus_all = sus_components.flatten()

    # 使用密度散点图
    sc1 = density_scatter(ax1, dft_all, sus_all, bins=100)
    cbar1 = fig1.colorbar(sc1, ax=ax1, label='density')
    cbar1.outline.set_linewidth(2)
    cbar1.set_label('Density', fontsize=25, fontweight='bold')
    cbar1.ax.tick_params(which='major', labelsize=20, width=2, length=10, direction='in')
    cbar1.ax.tick_params(which='minor', labelsize=12, width=2, length=6, direction='in')

    # 添加理想线
    lim_low = min(ax1.get_xlim()[0], ax1.get_ylim()[0])
    lim_high = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    ax1.plot([lim_low, lim_high], [lim_low, lim_high], '--', color='gray', linewidth=2)

    # 计算所有分量的统计指标 自己测试得指标
    # mae_all = np.mean(np.abs(dft_all - sus_all))
    # rmse_all = np.sqrt(np.mean((dft_all - sus_all) ** 2))
    # r2_all = r2_score(dft_all, sus_all)

    ax1.text(0.05, 0.95,
             f'MAE: {MAE:.3f} (GPa)\n'
             f'RMSE: {RMSE:.3f} (GPa)\n'
             f'R\u00B2: {R2:.3f}',
             transform=ax1.transAxes, va='top',
             bbox=dict(facecolor='white', alpha=.5), fontsize=20)

    ax1.set_xlabel('DFT Stress (GPa)', fontweight='bold')
    ax1.set_ylabel('SUS2 Stress (GPa)', fontweight='bold')
    plt.tight_layout()
    plt.savefig('stress_all_components_density.jpg', dpi=300)
    plt.close()

    # ===== 图2: 各分量分开的散点图 =====
    fig2, ax2 = plt.subplots(figsize=(10, 8))

    # 定义不同的颜色和标记
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
    markers = ['o', 's', '^', 'D', 'v', '<']

    # 绘制所有分量的散点图
    for i, (comp_name, color, marker) in enumerate(zip(component_names, colors, markers)):
        dft_comp = dft_components[:, i]
        sus_comp = sus_components[:, i]

        # 绘制散点图
        ax2.scatter(dft_comp, sus_comp, s=40,
                    c=color, marker=marker, label=comp_name, alpha=0.7)

    # 添加理想线
    all_dft = dft_components.flatten()
    all_sus = sus_components.flatten()
    min_val = min(np.min(all_dft), np.min(all_sus))
    max_val = max(np.max(all_dft), np.max(all_sus))
    margin = 0.1 * (max_val - min_val)

    ax2.plot([min_val - margin, max_val + margin],
             [min_val - margin, max_val + margin], 'k--', lw=2)

    # 设置图表属性
    ax2.set_xlabel('DFT Stress (GPa)',  fontweight='bold')
    ax2.set_ylabel('Predicted Stress (GPa)',fontweight='bold')
    ax2.legend()

    ax2.text(0.2, 0.95, f'MAE: {MAE:.3f} (GPa)\nRMSE: {RMSE:.3f} (GPa)\nR\u00B2: {R2:.3f}', fontsize=12,
             ha='left', va='top', transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))

    plt.tight_layout()
    plt.savefig('stress_components_comparison.jpg', dpi=300, bbox_inches='tight')
    plt.close()

    # ===== 图3: 各分量的密度散点图（子图形式） =====
    fig3, axs = plt.subplots(2, 3, figsize=(24, 16))
    axs = axs.flatten()

    for i, (ax, comp_name) in enumerate(zip(axs, component_names)):
        dft_comp = dft_components[:, i]
        sus_comp = sus_components[:, i]

        # 使用密度散点图
        sc = density_scatter(ax, dft_comp, sus_comp, bins=80)

        # 添加理想线
        lim_low = min(ax.get_xlim()[0], ax.get_ylim()[0])
        lim_high = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lim_low, lim_high], [lim_low, lim_high], '--', color='gray', linewidth=2)

        # 计算单个分量的统计指标
        mae_comp = np.mean(np.abs(dft_comp - sus_comp))
        rmse_comp = np.sqrt(np.mean((dft_comp - sus_comp) ** 2))
        r2_comp = r2_score(dft_comp, sus_comp)

        ax.text(0.05, 0.95,
                f'MAE: {mae_comp:.3f}\nRMSE: {rmse_comp:.3f}\nR²: {r2_comp:.3f}',
                transform=ax.transAxes, va='top', fontsize=14,
                bbox=dict(facecolor='white', alpha=.5))

        ax.set_xlabel('DFT Stress (GPa)', fontweight='bold')
        ax.set_ylabel('SUS2 Stress (GPa)', fontweight='bold')
        ax.set_title(f'Stress {comp_name}', fontweight='bold')

    # 调整子图布局，为colorbar留出空间
    plt.tight_layout()
    fig3.subplots_adjust(right=0.88)  # 为colorbar留出右边空间

    # 添加colorbar到最右边
    cbar_ax = fig3.add_axes([0.90, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    cbar = fig3.colorbar(sc, cax=cbar_ax, label='density')
    cbar.outline.set_linewidth(2)
    cbar.set_label('Density', fontsize=25, fontweight='bold')
    cbar.ax.tick_params(which='major', labelsize=20, width=2, length=10, direction='in')
    cbar.ax.tick_params(which='minor', labelsize=12, width=2, length=6, direction='in')

    plt.savefig('stress_components_density_subplots.jpg', dpi=300, bbox_inches='tight')
    plt.close()

    print("Stress analysis completed: 3 plots generated")


# ========= 综合图：能量、合力、应力合并 =========
def create_combined_plot(dft_xyz, sus_xyz, num_processes=None,R2=False):
    """创建一行三列的综合图：能量、合力、应力"""
    print("Creating combined plot...")

    # 计算能量数据
    dft = list(iread(dft_xyz))
    sus = list(iread(sus_xyz))

    # 能量数据
    with Pool(processes=num_processes) as pool:
        energy_res = pool.starmap(compute_energy, zip(dft, sus))
    e_dft, e_sus = map(np.array, zip(*energy_res))

    # 力数据
    with Pool(processes=num_processes) as pool:
        force_res = pool.starmap(compute_forces, zip(dft, sus))
    f_dft_vec = np.vstack([r[0] for r in force_res])
    f_sus_vec = np.vstack([r[1] for r in force_res])
    f_dft = np.linalg.norm(f_dft_vec, axis=1)
    f_sus = np.linalg.norm(f_sus_vec, axis=1)

    # 应力数据
    with Pool(processes=num_processes) as pool:
        stress_res = pool.starmap(compute_stress_Gpa, zip(dft, sus))
    dft_list = [result[0] for result in stress_res]
    sus_list = [result[1] for result in stress_res]

    # 提取应力分量
    dft_components = []
    sus_components = []
    for dft, sus in zip(dft_list, sus_list):
        dft_flat = [dft[0, 0], dft[1, 1], dft[2, 2], dft[0, 1], dft[0, 2], dft[1, 2]]
        sus_flat = [sus[0, 0], sus[1, 1], sus[2, 2], sus[0, 1], sus[0, 2], sus[1, 2]]
        dft_components.append(dft_flat)
        sus_components.append(sus_flat)

    dft_components = np.array(dft_components)
    sus_components = np.array(sus_components)
    dft_all = dft_components.flatten()
    sus_all = sus_components.flatten()

    # ===== 创建综合图 =====
    fig, axs = plt.subplots(1, 3, figsize=(30, 10))

    # 子图 (a): 能量
    # 方法1：在绘图前计算统一的范围
    sc1 = density_scatter(axs[0], e_dft, e_sus, bins=100)
    xlim = axs[0].get_xlim()
    ylim = axs[0].get_ylim()
    lim_low = min(xlim[0], ylim[0])
    lim_high = max(xlim[1], ylim[1])

    axs[0].plot([lim_low, lim_high], [lim_low, lim_high], '--', color='gray', linewidth=3)
    axs[0].set_xlim(lim_low, lim_high)
    axs[0].set_ylim(lim_low, lim_high)

    # # 获取x轴的刻度位置和标签
    # y_ticks = axs[0].get_yticks()
    # y_ticklabels = [tick.get_text() for tick in axs[0].get_yticklabels()]
    # print(y_ticklabels,y_ticks)
    # # 将x轴刻度设置到y轴
    # axs[0].set_yticks(y_ticks)
    # axs[0].set_yticklabels(y_ticklabels)
    # axs[0].set_xticks(y_ticks)
    # axs[0].set_xticklabels(y_ticklabels)

    mae_e = np.mean(np.abs(e_dft - e_sus)) * 1000
    rmse_e = np.sqrt(np.mean((e_dft - e_sus) ** 2)) * 1000
    r2_e = r2_score(e_dft, e_sus)
    if R2:
        axs[0].text(0.05, 0.9,
                    f'MAE: {mae_e:.3f} (meV/atom)\n'
                    f'RMSE: {rmse_e:.3f} (meV/atom)\n'
                    f'R\u00B2: {r2_e:.3f}',
                    transform=axs[0].transAxes, va='top',
                    bbox=dict(facecolor='white', alpha=.8), fontsize=20)
    else:
        axs[0].text(0.05, 0.9,
                    f'MAE: {mae_e:.3f} (meV/atom)\n'
                    f'RMSE: {rmse_e:.3f} (meV/atom)',
                    transform=axs[0].transAxes, va='top',
                    bbox=dict(facecolor='white', alpha=.8), fontsize=20)

    axs[0].set_xlabel('DFT Energy (eV/atom)', fontweight='bold')
    axs[0].set_ylabel('SUS2 Energy (eV/atom)', fontweight='bold')
    axs[0].text(0.03, 0.98, '(a)', transform=axs[0].transAxes,
                fontsize=32, fontweight='bold', va='top')
    axs[0].set_aspect('equal')

    # 子图 (b): 合力
    # 方法2：在绘图后获取范围并统一（推荐）
    sc2 = density_scatter(axs[1], f_dft, f_sus, bins=100)

    # 获取当前范围并统一
    xlim = axs[1].get_xlim()
    ylim = axs[1].get_ylim()
    lim_low = min(xlim[0], ylim[0])
    lim_high = max(xlim[1], ylim[1])

    axs[1].plot([lim_low, lim_high], [lim_low, lim_high], '--', color='gray', linewidth=3)
    axs[1].set_xlim(lim_low, lim_high)
    axs[1].set_ylim(lim_low, lim_high)

    mae_f = np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) / 3) * 1000
    rmse_f = np.sqrt(np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) ** 2) / 3) * 1000
    r2_f = r2_score(f_dft_vec, f_sus_vec)
    if R2:
        axs[1].text(0.05, 0.9,
                    f'MAE: {mae_f:.3f} (meV/Å)\n'
                    f'RMSE: {rmse_f:.3f} (meV/Å)\n'
                    f'R\u00B2: {r2_f:.3f}',
                    transform=axs[1].transAxes, va='top',
                    bbox=dict(facecolor='white', alpha=.8), fontsize=20)
    else:
        axs[1].text(0.05, 0.9,
                    f'MAE: {mae_f:.3f} (meV/Å)\n'
                    f'RMSE: {rmse_f:.3f} (meV/Å)',
                    transform=axs[1].transAxes, va='top',
                    bbox=dict(facecolor='white', alpha=.8), fontsize=20)

    axs[1].set_xlabel('DFT Force (eV/Å)', fontweight='bold')
    axs[1].set_ylabel('SUS2 Force (eV/Å)', fontweight='bold')
    axs[1].text(0.03, 0.98, '(b)', transform=axs[1].transAxes,
                fontsize=32, fontweight='bold', va='top')
    axs[1].set_aspect('equal')

    # 子图 (c): 应力（所有分量合并）
    # ===== 计算应力指标（使用与dft_sus_stress相同的计算方法） =====
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
    mae_s = total_frobenius_mae / N
    rmse_s = np.sqrt(total_frobenius_squared / N)
    r2_s = 1 - (SS_res / SS_tot) if SS_tot != 0 else float('nan')

    print(f'Stress Metrics - MAE: {mae_s:.3f} (GPa), RMSE: {rmse_s:.3f} (GPa), R²: {r2_s:.3f}')

    sc3 = density_scatter(axs[2], dft_all, sus_all, bins=100)

    # 同样使用绘图后统一范围的方法
    xlim = axs[2].get_xlim()
    ylim = axs[2].get_ylim()
    lim_low = min(xlim[0], ylim[0])
    lim_high = max(xlim[1], ylim[1])

    axs[2].plot([lim_low, lim_high], [lim_low, lim_high], '--', color='gray', linewidth=3)
    axs[2].set_xlim(lim_low, lim_high)
    axs[2].set_ylim(lim_low, lim_high)

    if R2:
        axs[2].text(0.05, 0.9,
                    f'MAE: {mae_s:.3f} (GPa)\n'
                    f'RMSE: {rmse_s:.3f} (GPa)\n'
                    f'R\u00B2: {r2_s:.3f}',
                    transform=axs[2].transAxes, va='top',
                    bbox=dict(facecolor='white', alpha=.8), fontsize=20)
    else:
        axs[2].text(0.05, 0.9,
                    f'MAE: {mae_s:.3f} (GPa)\n'
                    f'RMSE: {rmse_s:.3f} (GPa)',
                    transform=axs[2].transAxes, va='top',
                    bbox=dict(facecolor='white', alpha=.8), fontsize=20)

    axs[2].set_xlabel('DFT Stress (GPa)', fontweight='bold')
    axs[2].set_ylabel('SUS2 Stress (GPa)', fontweight='bold')
    axs[2].text(0.03, 0.98, '(c)', transform=axs[2].transAxes,
                fontsize=32, fontweight='bold', va='top')
    axs[2].set_aspect('equal')

    # 添加统一的colorbar
    plt.tight_layout()
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(sc2, cax=cbar_ax, label='density')
    cbar.outline.set_linewidth(2)
    cbar.set_label('Density', fontsize=28, fontweight='bold')
    cbar.ax.tick_params(which='major', labelsize=22, width=2, length=10, direction='in')
    cbar.ax.tick_params(which='minor', labelsize=14, width=2, length=6, direction='in')

    plt.savefig('efs.jpg', dpi=300, bbox_inches='tight')
    plt.close()

# ========= 主程序 =========
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dft_sus2_compare.py dft.xyz sus2.xyz")
        sys.exit(1)

    dft_xyz, sus_xyz = sys.argv[1], sys.argv[2]
    print(f'dft_data: {dft_xyz}  sus_data: {sus_xyz}')
    #
    # t0 = time.time()
    # dft_sus_e(dft_xyz, sus_xyz, num_processes)
    # print(f'Energy done  in {(time.time() - t0) / 60:.2f} min')
    #
    # t0 = time.time()
    # f_dft_vec, f_sus_vec = dft_sus_f(dft_xyz, sus_xyz, num_processes)
    # print(f'|F|   done  in {(time.time() - t0) / 60:.2f} min')
    #
    # t0 = time.time()
    # dft_sus_f_three(dft_xyz, sus_xyz, num_processes)
    # print(f'FxFyFz done  in {(time.time() - t0) / 60:.2f} min')
    #
    # t0 = time.time()
    # dft_sus_stress(dft_xyz, sus_xyz, num_processes)
    # print(f'Stress done  in {(time.time() - t0) / 60:.2f} min')

    """9分量分开应力MAE,RMSE,R2计算方式与合起来MAE,RMSE,R2不同的,以合起来的为标准。合起来的与SUS2程序定义相同"""

    t0 = time.time()
    create_combined_plot(dft_xyz, sus_xyz, num_processes,R2=True)
    print(f'total_efs done  in {(time.time() - t0) / 60:.2f} min')