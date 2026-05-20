import numpy as np
import matplotlib.pyplot as plt

# ---------------------- 直接读取你合并好的文件 ----------------------
path_file = "/Users/wenyuyang/Documents/论文写作/bridgeJFR/CanningParameterization-based planning/pillarPathXYZYaw.txt"

# 读取 x y z yaw
data = np.loadtxt(path_file)
x   = data[:, 0]
y   = data[:, 1]
z   = data[:, 2]
yaw = data[:, 3]
yaw = yaw

# 根据 yaw 计算朝向向量（只用到 yaw）
yaw_rad = np.deg2rad(yaw)
# dx = np.cos(yaw_rad)
# dy = np.sin(yaw_rad)
dx = np.sin(yaw_rad)  # 正东
dy = np.cos(yaw_rad)  # 正北
dz = np.zeros_like(dx)

arrow_len = 0.9  # 箭头长度

# ---------------------- 绘图 ----------------------
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 画路径
ax.plot(x, y, z, linewidth=2, color='#2E86AB', label='3D Path')

# 画朝向箭头（每隔2个点画一个，避免太密）
for i in range(0, len(x), 2):
    ax.quiver(
        x[i], y[i], z[i],
        dx[i] * arrow_len,
        dy[i] * arrow_len,
        dz[i] * arrow_len,
        color='#E63946',
        linewidth=1.5
    )

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Path with Yaw Direction (Equal Scale)')
ax.legend()

# ---------------------- 强制等比例显示 ----------------------
max_range = max(np.ptp(x), np.ptp(y), np.ptp(z)) / 2.0
mid_x = (x.max() + x.min()) * 0.5
mid_y = (y.max() + y.min()) * 0.5
mid_z = (z.max() + z.min()) * 0.5

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.tight_layout()
plt.show()