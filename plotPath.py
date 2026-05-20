import numpy as np
import matplotlib.pyplot as plt

# -------------------------- 文件路径 --------------------------
path_file = "/Users/wenyuyang/Documents/论文写作/bridgeJFR/CanningParameterization-based planning/pillar_path.txt"
cmd_file  = "/Users/wenyuyang/Documents/论文写作/bridgeJFR/CanningParameterization-based planning/pillar_CMD.txt"

# -------------------------- 读取数据 --------------------------
# 路径点 x, y, z
points = np.loadtxt(path_file)
x = points[:, 0]
y = points[:, 1]
z = points[:, 2]

# 姿态 roll, pitch, yaw
rpy = np.loadtxt(cmd_file)
yaw = rpy[:, 2]  # 只用 yaw

# 保证点数对齐（取短的一方）
N = min(len(x), len(yaw))
x, y, z = x[:N], y[:N], z[:N]
yaw = yaw[:N]

# -------------------------- 计算朝向向量 --------------------------
# yaw 是绕 Z 轴的角度，转弧度
yaw_rad = np.deg2rad(yaw)

# 朝向 = 单位前向向量 (cos(yaw), sin(yaw), 0)
dir_x = np.cos(yaw_rad)
dir_y = np.sin(yaw_rad)
dir_z = np.zeros_like(dir_x)

# 箭头长度（按你的场景调整）
arrow_len = 0.3

# -------------------------- 绘图 --------------------------
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 画路径
ax.plot(x, y, z, linewidth=2, color='#2E86AB', label='Path')

# 画每个点的朝向箭头
for i in range(0, N, 2):  # 每隔2个点画一个箭头，避免太密
    ax.quiver(x[i], y[i], z[i],
              dir_x[i]*arrow_len,
              dir_y[i]*arrow_len,
              dir_z[i]*arrow_len,
              color='#E63946', linewidth=1.5, length=1)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Path + Yaw Direction')
ax.legend()

# -------------------------- 等比例坐标轴 --------------------------
max_range = max(np.ptp(x), np.ptp(y), np.ptp(z)) / 2.0
mid_x = (x.max() + x.min()) * 0.5
mid_y = (y.max() + y.min()) * 0.5
mid_z = (z.max() + z.min()) * 0.5

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.tight_layout()
plt.show()