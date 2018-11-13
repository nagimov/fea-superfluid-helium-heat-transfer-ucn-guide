import numpy as np
from matplotlib import pyplot as plt

def crop_mesh(mesh):
    true_points = np.argwhere(mesh)  # argwhere gives coordinates of every non-zero point
    top_left = true_points.min(axis=0)
    bottom_right = true_points.max(axis=0)
    out = mesh[top_left[0]:bottom_right[0]+1,  # plus 1 because slice isn't
               top_left[1]:bottom_right[1]+1]  # inclusive
    out[out==0] = np.nan
    return out


# read the solution
mesh_t = np.loadtxt('mesh_t_post.csv', delimiter=',', dtype=np.float)
mesh_m = np.loadtxt('mesh_m.csv', delimiter=',', dtype=np.float)

# replace zeros with NANs so matplotlib can ignore these cells
mesh_t[mesh_t==0] = np.nan

# rotate mesh_t
mesh_t = mesh_t[::-1,:]
mesh_m = mesh_m[::-1,:]

# plotting...
plt.figure()
# ...whole system
plt.subplot(3, 1, 1)
plt.imshow(mesh_t, cmap='jet', interpolation='nearest', aspect='auto')
plt.colorbar()
plt.tight_layout()
# ...helium-4 only
plt.subplot(3, 1, 2)
mesh_t_he4 = np.copy(mesh_t)
mesh_t_he4[mesh_m!=4] = 0
plt.imshow(crop_mesh(mesh_t_he4), cmap='jet', interpolation='nearest', aspect='auto')
plt.colorbar()
plt.tight_layout()
# ...copper fins only
plt.subplot(3, 1, 3)
mesh_t_cu = np.copy(mesh_t)
mesh_t_cu[mesh_m!=2] = 0
plt.imshow(crop_mesh(mesh_t_cu), cmap='jet', interpolation='nearest', aspect='auto')
plt.colorbar()
plt.tight_layout()

plt.show()
