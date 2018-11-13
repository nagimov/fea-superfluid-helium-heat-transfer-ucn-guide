import numpy as np
import csv
import os
from matplotlib import pyplot as plt
import random

print('# creating mesh...')

# input geometry
# all dimensions are in SI units (a.k.a. "meters")
r_4he = 0.15/2  # helium guide radius (a.k.a. half diameter)
l_hx = 0.5  # length of heat exchanger
l_guide = 2  # length of helium guide
t_cu = 0.001  # thickness of heat exchanger base
h_fin = 0.003  # height of fins
t_fin = 0.001  # thickness of fins
p_fin = 0.002  # pitch of fins
t_3he = 0.001  # thickness of helium-3 layer on top of fins
assert(int(l_hx / p_fin) == l_hx / p_fin)

init_he4_temp = 1.2
init_he3_temp = 0.7

# sizes of mesh elements
# fine mesh
x_size_1 = 0.0005
y_size_1 = 0.001
x_size_2 = x_size_1
y_size_2 = 0.0005
x_size_3 = x_size_2
y_size_3 = 0.0005
x_size_4 = x_size_3
y_size_4 = 0.0001
x_size_5 = 0.005
y_size_5 = y_size_1
x_size_6 = x_size_5
'''
# rough mesh
x_size_1 = 0.001
y_size_1 = 0.01
x_size_2 = x_size_1
y_size_2 = 0.001
x_size_3 = x_size_2
y_size_3 = 0.001
x_size_4 = x_size_3
y_size_4 = 0.0005
x_size_5 = 0.025
y_size_5 = y_size_1
x_size_6 = x_size_5
'''

# meshing:
# first index - length (x), second index - height (y)
# origin - at the left bottom point of HX
l_total = l_hx + l_guide
h_total = t_3he + h_fin + t_cu + r_4he
l_mesh = round(l_hx / x_size_1) + round(l_guide / x_size_5)  # length
h_mesh = round(r_4he / y_size_1) + round(t_cu / y_size_2) + round(h_fin / y_size_3) + round(t_3he / y_size_4)  # height
# initializing empty mesh arrays
mesh_m = np.ones([l_mesh, h_mesh], dtype=np.int)  # element material (1 - 'empty', 2 - 'cu', 3 - '3he', 4 - '4he')
mesh_l = np.zeros([l_mesh, h_mesh])  # element length
mesh_h = np.zeros([l_mesh, h_mesh])  # element height
mesh_t = np.zeros([l_mesh, h_mesh], dtype=np.float64)  # element temperature
mesh_t[:, :] = np.nan
mesh_r = np.zeros([l_mesh, h_mesh])  # radius from symmetry axis to the center of the element
mesh_p = np.zeros([l_mesh, h_mesh], dtype=np.bool)  # element is processed (false - not processed yet, true - element processed)
mesh_v = np.zeros([l_mesh, h_mesh])  # volume of the cell
mesh_q = np.zeros([l_mesh, h_mesh])  # external heat


# meshing stuff (mesh is spaced unevenly)
# helium guide under hx
m_block = 4
t_block = init_he4_temp
h_block = r_4he
l_block = l_hx
l_mesh_block = round(l_block / x_size_1)
h_mesh_block = round(h_block / y_size_1)
x_block_a = 0
x_block_b = x_block_a + l_mesh_block
y_block_a = 0
y_block_b = y_block_a + h_mesh_block
r_block_lo = y_size_1 / 2
r_block_hi = r_block_lo + h_block
r_block = np.linspace(r_block_lo, r_block_hi, num=h_mesh_block, endpoint=True)
mesh_m[x_block_a:x_block_b, y_block_a:y_block_b] = m_block
mesh_l[x_block_a:x_block_b, y_block_a:y_block_b] = l_block / l_mesh_block
mesh_h[x_block_a:x_block_b, y_block_a:y_block_b] = h_block / h_mesh_block
mesh_t[x_block_a:x_block_b, y_block_a:y_block_b] = t_block
mesh_r[x_block_a:x_block_b, y_block_a:y_block_b] = r_block
# copper base of heat exchanger
m_block = 2
t_block = init_he4_temp
h_block = t_cu
l_block = l_hx
l_mesh_block = round(l_block / x_size_2)
h_mesh_block = round(h_block / y_size_2)
x_block_a = 0
x_block_b = x_block_a + l_mesh_block
y_block_a = y_block_b
y_block_b = y_block_a + h_mesh_block
r_block_lo = r_block_hi
r_block_hi = r_block_lo + h_block
r_block = np.linspace(r_block_lo, r_block_hi, num=h_mesh_block, endpoint=True)
mesh_m[x_block_a:x_block_b, y_block_a:y_block_b] = m_block
mesh_l[x_block_a:x_block_b, y_block_a:y_block_b] = l_block / l_mesh_block
mesh_h[x_block_a:x_block_b, y_block_a:y_block_b] = h_block / h_mesh_block
mesh_t[x_block_a:x_block_b, y_block_a:y_block_b] = t_block
mesh_r[x_block_a:x_block_b, y_block_a:y_block_b] = r_block
# heat exchanger fins
m_block_1 = 2
m_block_2 = 3
t_block = init_he4_temp
t_block_pumping = init_he3_temp
h_block = h_fin
l_block = l_hx
l_block_1 = t_fin
l_block_2 = p_fin - t_fin
h_mesh_block = round(h_block / y_size_3)
l_mesh_block = round(l_block / x_size_3)
l_mesh_block_1 = l_block_1 / x_size_3
l_mesh_block_2 = l_block_2 / x_size_3
x_block_a = 0
x_block_b = x_block_a + l_mesh_block
y_block_a = y_block_b
y_block_b = y_block_a + h_mesh_block
r_block_lo = r_block_hi
r_block_hi = r_block_lo + h_block
r_block = np.linspace(r_block_lo, r_block_hi, num=h_mesh_block, endpoint=True)
hx_periods = int(l_hx / p_fin)
l_mesh_period = l_block / x_size_3 / hx_periods
mesh_t[x_block_a:x_block_b, y_block_a:y_block_b] = t_block_pumping
for i in range(0, hx_periods):
    x_block_a_1 = round(x_block_a + i * l_mesh_period)
    x_block_b_1 = round(x_block_a + i * l_mesh_period + l_mesh_block_1)
    mesh_m[x_block_a_1:x_block_b_1, y_block_a:y_block_b] = m_block_1
    x_block_a_2 = round(x_block_a + i * l_mesh_period + l_mesh_block_1)
    x_block_b_2 = round(x_block_a + (i + 1) * l_mesh_period)
    mesh_m[x_block_a_2:x_block_b_2, y_block_a:y_block_b] = m_block_2
    mesh_t_gradient = np.linspace(t_block, t_block_pumping, num=h_mesh_block, endpoint=True)
    mesh_t[x_block_a_1:x_block_b_1, y_block_a:y_block_b] = mesh_t_gradient
mesh_l[x_block_a:x_block_b, y_block_a:y_block_b] = l_block / l_mesh_block
mesh_h[x_block_a:x_block_b, y_block_a:y_block_b] = h_block / h_mesh_block
mesh_r[x_block_a:x_block_b, y_block_a:y_block_b] = r_block
# helium-3 layer
m_block = 3
t_block_pumping = init_he3_temp
h_block = t_3he
l_block = l_hx
l_mesh_block = round(l_block / x_size_4)
h_mesh_block = round(h_block / y_size_4)
x_block_a = 0
x_block_b = x_block_a + l_mesh_block
y_block_a = y_block_b
y_block_b = y_block_a + h_mesh_block
r_block_lo = r_block_hi
r_block_hi = r_block_lo + h_block
r_block = np.linspace(r_block_lo, r_block_hi, num=h_mesh_block, endpoint=True)
mesh_m[x_block_a:x_block_b, y_block_a:y_block_b] = m_block
mesh_l[x_block_a:x_block_b, y_block_a:y_block_b] = l_block / l_mesh_block
mesh_h[x_block_a:x_block_b, y_block_a:y_block_b] = h_block / h_mesh_block
mesh_t[x_block_a:x_block_b, y_block_a:y_block_b] = t_block_pumping
mesh_r[x_block_a:x_block_b, y_block_a:y_block_b] = r_block
# helium-4 guide between UCN bottle and HX
m_block = 4
t_block = init_he4_temp
h_block = r_4he
l_block = l_guide
l_mesh_block = round(l_block / x_size_5)
h_mesh_block = round(h_block / y_size_5)
x_block_a = round(l_hx / x_size_1)
x_block_b = x_block_a + l_mesh_block
y_block_a = 0
y_block_b = y_block_a + h_mesh_block
r_block_lo = y_size_5 / 2
r_block_hi = r_block_lo + h_block
r_block = np.linspace(r_block_lo, r_block_hi, num=h_mesh_block, endpoint=True)
mesh_m[x_block_a:x_block_b, y_block_a:y_block_b] = m_block
mesh_l[x_block_a:x_block_b, y_block_a:y_block_b] = l_block / l_mesh_block
mesh_h[x_block_a:x_block_b, y_block_a:y_block_b] = h_block / h_mesh_block
mesh_t[x_block_a:x_block_b, y_block_a:y_block_b] = t_block
mesh_r[x_block_a:x_block_b, y_block_a:y_block_b] = r_block

print('# mesh {} x {} created'.format(l_mesh, h_mesh))


print('# analyzing mesh...')

'''
make lists of for finite elements operations
--------------------------------------------
Operation is presented as a tuple ((x1, y1), (x2, y2), dV1, dV2, dL, dA, op_type)
    * (x1, y1) is a coordinate of the first cell
    * (x2, y2) is a corrdinate of the second cell
    * dV1 and dV2 are the volumes of cells
    * dL is the distance between centers of two cells
    * dA is the surface area of interface between two cells
    * op_type is marked as two-digit number indicating type of both cells
        - order of digits in cell operation is cell-centered, e.g. for helium-4
        cell neighboring with copper cell, cell operation will be marked as '42'
        - this coincidence has nothing to do with the answer to the ultimate
        question of life
e.g. operation in between copper cell (2, 2) and helium-4 cell (2, 3) can be
written as ((2, 2), (2, 3), 1.0, 2.0, 3.0, 4.0, '24') or ((2, 3), (2, 2), 2.0, 1.0, 3.0, 4.0, '42')
There are absolutely no optimizations here, since this analysis runs only once!
Code is written very explicitly on purpose, please do not fold this into a
function called for each pair of cells.
'''

def ext_point(np_array, in_tuple):
    return np_array[in_tuple[0], in_tuple[1]]

# I feel uncomfortable without this, so just live with it
dfi = 1.0  # this is the only case where unit doesn't matter
           # (well, it kinda affects iteration timestep)

# calculate volume of each cell
for i in range(l_mesh):  # x
    for j in range(h_mesh):  # y
        cell = (i, j)
        if ext_point(mesh_m, cell) != 1:  # cell is not empty
            mesh_v[i, j] = ext_point(mesh_l, cell) * ext_point(mesh_h, cell) * dfi * ext_point(mesh_r, cell)
        else:
            mesh_v[i, j] = 0

mesh_ops = []
for i in range(l_mesh):  # x
    for j in range(h_mesh):  # y
        if not mesh_p[i, j]:  # cell not processed yet
            mesh_p[i, j] = True  # mark cell as processed
            if mesh_m[i, j] != 1:  # cell is not empty
                cell_current = (i, j)
                neighbors = []
                dA = []
                dL = []
                # left neighbor
                if (i != 0):  # not on the left boundary of the mesh
                    neighbors.append((i-1, j))
                # right neighbor
                if (i != l_mesh - 1):  # not on the right boundary of the mesh
                    neighbors.append((i+1, j))
                # upper neighbor
                if (j != 0):  # not on the upper boundary of the mesh
                    neighbors.append((i, j-1))
                # lower neighbor
                if (j != h_mesh - 1):  # not on the lower boundary of the mesh
                    neighbors.append((i, j+1))
                for cell_proc in neighbors:
                    if not ext_point(mesh_p, cell_proc):  # neighbor cell not processed yet
                        if ext_point(mesh_m, cell_proc) != 1:  # neighbor cell is not empty
                            op_type = str(ext_point(mesh_m, cell_current)) + \
                                      str(ext_point(mesh_m, cell_proc))
                            if cell_current[1] == cell_proc[1]:  # neighbors are horizontal
                                dL = 0.5 * (ext_point(mesh_l, cell_current) + ext_point(mesh_l, cell_proc))
                                dA = ext_point(mesh_h, cell_current) * dfi * ext_point(mesh_r, cell_current)
                            else:
                                if cell_current[0] == cell_proc[0]:  # neighbors are vertical
                                    dL = 0.5 * (ext_point(mesh_h, cell_current) + ext_point(mesh_h, cell_proc))
                                    dA = ext_point(mesh_l, cell_current) * dfi * ext_point(mesh_r, cell_current)
                                else: print('mesh is screwed'); 1/0
                            mesh_ops.append((cell_current, cell_proc, dL, dA, op_type))  # mesh_ops data structure
print('# mesh is analyzed, {} operable cell pairs are found'.format(len(mesh_ops)))


# cells where to apply some heat
Q = 10  # W, total heat deposition
l_q = 0.5  # m, length of region to apply heat
r_q = 0.05  # m, radius of region to apply heat
l_q_cells = round(l_q/x_size_5)
r_q_cells = round(r_q/y_size_5)
V_q = l_q * r_q**2 * np.pi
for i_heat in range(l_q_cells):
    for j_heat in range(r_q_cells):
        i = l_mesh - 1 - i_heat
        j = j_heat
        mesh_q[i][j] = Q / V_q * mesh_v[i][j]

print('# saving mesh info...')
with open('mesh_ops.csv', 'w', newline='\n') as f:
    writer = csv.writer(f)
    for mesh_op in mesh_ops:
        line = [mesh_op[0][0], mesh_op[0][1], mesh_op[1][0], mesh_op[1][1]] + list(mesh_op[2:])
        writer.writerow(line)
np.savetxt('mesh_t.csv', np.transpose(mesh_t), fmt='%.6f', delimiter=',')
np.savetxt('mesh_m.csv', np.transpose(mesh_m), fmt='%.6f', delimiter=',')
np.savetxt('mesh_v.csv', np.transpose(mesh_v), fmt='%.14f', delimiter=',')
np.savetxt('mesh_q.csv', np.transpose(mesh_q), fmt='%.14f', delimiter=',')
filesize = 0
for file in ['mesh_ops.csv', 'mesh_t.csv', 'mesh_m.csv', 'mesh_v.csv', 'mesh_q.csv']:
    filesize += os.path.getsize(file)
print('# mesh info saved, {:.2f}MB total'.format(filesize / 1024.0 / 1024.0))

# pass some info to C solver via header
# so static arrays can be defined
with open('params.h', 'w', newline='\n') as f:
    f.write('#define L_MESH {}\n'.format(l_mesh))
    f.write('#define H_MESH {}\n'.format(h_mesh))
    f.write('#define OPS_MESH {}\n'.format(len(mesh_ops)))
