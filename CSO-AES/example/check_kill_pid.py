import os
import psutil
import sys


def get_child_processes(parent_pid):
    """
    获取指定父进程的所有子进程

    参数:
        parent_pid (int): 父进程的进程 ID

    返回:
        list: 包含所有子进程的列表，每个元素是一个 psutil.Process 对象
    """
    # 获取父进程对象
    parent = psutil.Process(parent_pid)
    # 获取所有子进程（包括子孙进程）
    children = parent.children(recursive=True)  # recursive=True 包含所有子孙进程
    return children


def main(main_id, mode):
    """
    主函数，根据模式检查或杀死指定父进程的所有子进程

    参数:
        main_id (int): 父进程的进程 ID
        mode (str): 操作模式，'check' 表示检查子进程，'kill' 表示杀死子进程
    """
    # 获取所有子进程
    child_processes = get_child_processes(main_id)

    if mode == 'check':
        print(f"main PID: {main_id}")
        # 检查模式：打印每个子进程的 PID 和命令
        for child in child_processes:
            print(f"child PID: {child.pid}, command: {child.name()}")
    elif mode == 'kill':
        # 杀死模式：杀死每个子进程
        for child in child_processes:
            os.system(f'kill {child.pid}')
        os.system(f'kill {main_id}')
    else:
        # 如果模式不是 'check' 或 'kill'，抛出错误
        raise ValueError('Usage: python check_kill_pid.py <mode> (mode == check or kill)')


if __name__ == "__main__":
    # 检查是否存在 pid.txt 文件
    if not os.path.exists('pid.txt'):
        raise ValueError('pid.txt file does not exist!')

    # 读取 pid.txt 文件的最后一行，获取父进程的 PID
    with open('pid.txt', 'r') as f:
        lines = f.readlines()
    main_id = int(lines[-1].strip())  # 确保去除行尾的换行符

    # 从命令行参数获取操作模式
    if len(sys.argv) != 2:
        raise ValueError('Usage: python check_kill_pid.py <mode> (mode == check or kill)')
    mode = sys.argv[1]

    # 调用主函数
    main(main_id, mode)
