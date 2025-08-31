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
import time, sys
from sklearn.metrics import r2_score
from joblib import Parallel, delayed

# ========= 开关区 =========
SAVE_DATA     = False       # True → 额外保存 3 个 CSV
num_processes = 10          # 并行进程数
cmap_name     = 'Spectral_r'
scatter_size  = 10          # 每个点半径


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
    order = z.argsort()                # 低密度点先画
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

def select_abs(a, b):      # 保留原作者的小辅助
    return a if abs(a) > abs(b) else b

# ========= 能量 =========
def dft_sus_e(dft_xyz, sus_xyz, num_processes=None):
    t0 = time.time()
    dft = list(iread(dft_xyz)); sus = list(iread(sus_xyz))
    print(f'read_xyz_time: {time.time() - t0:.1f}s  frames={len(dft)}')

    with Pool(processes=num_processes) as pool:
        res = pool.starmap(compute_energy, zip(dft, sus))
    e_dft, e_sus = map(np.array, zip(*res))

    # ===== 绘图 =====
    fig, ax = plt.subplots()
    sc = density_scatter(ax, e_dft, e_sus, bins=100)
    fig.colorbar(sc, ax=ax, label='density')

    lim_low = min(ax.get_xlim()[0], ax.get_ylim()[0])
    lim_high = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lim_low, lim_high], [lim_low, lim_high], '--') #灰阶

    mae  = np.mean(np.abs(e_dft - e_sus)) * 1000
    rmse = np.sqrt(np.mean((e_dft - e_sus) ** 2)) * 1000
    r2   = r2_score(e_dft, e_sus)
    ax.text(0.05, 0.95,
            f'MAE: {mae:.3f} (meV/atom)\n'
            f'RMSE: {rmse:.3f} (meV/atom)\n'
            f'R\u00B2: {r2:.3f}',
            transform=ax.transAxes, va='top',
            bbox=dict(facecolor='white', alpha=.5))

    ax.set_xlabel('DFT Energy (eV/atom)',  fontweight='bold')
    ax.set_ylabel('SUS\u00B2 Energy (eV/atom)', fontweight='bold')
    plt.tight_layout(); plt.savefig('energy.png', dpi=300); plt.close()

    if SAVE_DATA:
        np.savetxt('energy_data.csv',
                   np.column_stack((e_dft, e_sus)),
                   delimiter=',', header='DFT_E,SUS2_E', comments='')

# ========= 总力 =========
def dft_sus_f(dft_xyz, sus_xyz, num_processes=None):
    t0 = time.time()
    dft = list(iread(dft_xyz)); sus = list(iread(sus_xyz))
    print(f'read_xyz_time: {time.time() - t0:.1f}s')

    with Pool(processes=num_processes) as pool:
        res = pool.starmap(compute_forces, zip(dft, sus))

    f_dft_vec = np.vstack([r[0] for r in res])
    f_sus_vec = np.vstack([r[1] for r in res])
    f_dft = np.linalg.norm(f_dft_vec, axis=1)
    f_sus = np.linalg.norm(f_sus_vec, axis=1)

    # ===== 绘图 =====
    fig, ax = plt.subplots()
    sc = density_scatter(ax, f_dft, f_sus, bins=100)
    fig.colorbar(sc, ax=ax, label='density')

    lim_low = min(ax.get_xlim()[0], ax.get_ylim()[0])
    lim_high = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lim_low, lim_high], [lim_low, lim_high], '--')

    mae  = np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) / 3) * 1000
    rmse = np.sqrt(np.mean(np.linalg.norm(f_dft_vec - f_sus_vec, axis=1) ** 2) / 3) * 1000
    r2   = r2_score(f_dft_vec, f_sus_vec)
    ax.text(0.05, 0.95,
            f'MAE: {mae:.3f} (meV/Å)\n'
            f'RMSE: {rmse:.3f} (meV/Å)\n'
            f'R\u00B2: {r2:.3f}',
            transform=ax.transAxes, va='top',
            bbox=dict(facecolor='white', alpha=.5))

    ax.set_xlabel('DFT Force (eV/Å)',  fontweight='bold')
    ax.set_ylabel('SUS\u00B2 Force (eV/Å)', fontweight='bold')
    plt.tight_layout(); plt.savefig('force.png', dpi=300); plt.close()

    if SAVE_DATA:
        np.savetxt('force_data.csv',
                   np.column_stack((f_dft, f_sus)),
                   delimiter=',', header='|F_DFT|,|F_SUS2|', comments='')

    return f_dft_vec, f_sus_vec

# ========= Fx/Fy/Fz =========
def dft_sus_f_three(dft_xyz, sus_xyz, num_processes=None):
    t0 = time.time()
    dft = list(iread(dft_xyz)); sus = list(iread(sus_xyz))
    print(f'read_xyz_time: {time.time() - t0:.1f}s')

    with Pool(processes=num_processes) as pool:
        res = pool.starmap(compute_forces, zip(dft, sus))

    f_dft_vec = np.vstack([r[0] for r in res])
    f_sus_vec = np.vstack([r[1] for r in res])

    names = ['Fx', 'Fy', 'Fz']
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))

    for i, ax in enumerate(axs):
        x = f_dft_vec[:, i]; y = f_sus_vec[:, i]
        sc = density_scatter(ax, x, y, bins=80)
        lim_low = min(ax.get_xlim()[0], ax.get_ylim()[0])
        lim_high = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lim_low, lim_high], [lim_low, lim_high], '--')

        mae  = np.mean(np.abs(x - y)) * 1000
        rmse = np.sqrt(np.mean((x - y) ** 2)) * 1000
        r2   = r2_score(x, y)
        ax.text(0.05, 0.95,
                f'MAE: {mae:.1f}\nRMSE: {rmse:.1f}\nR²: {r2:.3f}',
                transform=ax.transAxes, va='top',
                bbox=dict(facecolor='white', alpha=.5))
        ax.set_xlabel('DFT Force (eV/Å)')
        ax.set_ylabel('SUS\u00B2 Force (eV/Å)')
        ax.set_title(names[i])

        # 只在中间子图加一次 colorbar
        if i == 1:
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label='density')

    plt.tight_layout(); plt.savefig('force_three.png', dpi=300); plt.close()

    if SAVE_DATA:
        np.savetxt('force_components.csv',
                   np.column_stack((f_dft_vec, f_sus_vec)),
                   delimiter=',',
                   header='Fx_DFT,Fy_DFT,Fz_DFT,Fx_SUS2,Fy_SUS2,Fz_SUS2',
                   comments='')

# ========= 主程序 =========
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dft_sus2_compare.py dft.xyz sus2.xyz")
        sys.exit(1)

    dft_xyz, sus_xyz = sys.argv[1], sys.argv[2]
    print(f'dft_data: {dft_xyz}  sus_data: {sus_xyz}')

    t0 = time.time()
    dft_sus_e(dft_xyz, sus_xyz, num_processes)
    print(f'Energy done  in {(time.time()-t0)/60:.2f} min')

    t0 = time.time()
    f_dft_vec, f_sus_vec = dft_sus_f(dft_xyz, sus_xyz, num_processes)
    print(f'|F|   done  in {(time.time()-t0)/60:.2f} min')

    t0 = time.time()
    dft_sus_f_three(dft_xyz, sus_xyz, num_processes)
    print(f'FxFyFz done  in {(time.time()-t0)/60:.2f} min')
