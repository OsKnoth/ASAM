from netCDF4 import Dataset
import numpy as np
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from shapely.geometry import Polygon
from descartes import PolygonPatch
mpl.use('Agg')


nfiles = 6
file_base = "output"

z_slice = 1e-20


solutions = []
gradients_u = []
gradients_v = []
gradients_w = []
partitions = []
x2_arrs = []
y2_arrs = []
z2_arrs = []


dataset = Dataset("CanyonHof.nc")
x2_glob = dataset['x2glob'][:]
y2_glob = dataset['y2glob'][:]
z2_glob = dataset['z2glob'][:]
dataset.close()

x_glob = 0.5 * (x2_glob[1:] + x2_glob[:-1])
y_glob = 0.5 * (y2_glob[1:] + y2_glob[:-1])
z_glob = 0.5 * (z2_glob[1:] + z2_glob[:-1])

nz = z_glob.size
ny = y_glob.size
nx = x_glob.size

solution = np.empty([nz, ny, nx])
block_bnds = []
block_ids = []
block_ranks = []
gradient_u = np.empty([nz, ny, nx + 1])
gradient_v = np.empty([nz, ny + 1, nx])
gradient_w = np.empty([nz + 1, ny, nx])

for i in range(nfiles):
    dataset = Dataset('{}{:03d}.nc'.format(file_base, i), 'r')
    nblocks = dataset.dimensions['block'].size
    blockids = dataset['blockids'][:]

    for j in range(nblocks):
        z_arr = dataset['zb' + str(blockids[j])][:]
        z2_arr = dataset['z2b' + str(blockids[j])][:]
        y2_arr = dataset['y2b' + str(blockids[j])][:]
        x2_arr = dataset['x2b' + str(blockids[j])][:]

        zind = np.argmin(np.absolute(z_arr - z_slice))
        sol = dataset['rho_b' + str(blockids[j])][:, :, :]
        sol[sol < 1] = np.nan
        nz_sub, ny_sub, nx_sub = sol.shape[:]
        grad_u = dataset['rhovel_u_b' + str(blockids[j])][:]
        grad_v = dataset['rhovel_v_b' + str(blockids[j])][:]
        grad_w = dataset['rhovel_w_b' + str(blockids[j])][:]

        if z2_arr[0] < z_slice and z2_arr[-1] >= z_slice:
            solutions.append(sol[zind])
            gradients_u.append(grad_u[zind])
            gradients_v.append(grad_v[zind])
            gradients_w.append(grad_w[zind])
            partitions.append(np.full_like(solutions[-1], i))
            block_bnds.append(Polygon([(x2_arr[0], y2_arr[0]),(x2_arr[-1], y2_arr[0]),(x2_arr[-1], y2_arr[-1]),(x2_arr[0], y2_arr[-1]),(x2_arr[0], y2_arr[0])]))
            block_ids.append(blockids[j])
            block_ranks.append(i)
            z2_arrs.append(z2_arr)
            x2_arrs.append(x2_arr)
            y2_arrs.append(y2_arr)
        
        iz = np.argmin(np.absolute(z2_glob - z2_arr[0]))
        iy = np.argmin(np.absolute(y2_glob - y2_arr[0]))
        ix = np.argmin(np.absolute(x2_glob - x2_arr[0]))

        try:
            solution[iz:iz + nz_sub, iy:iy + ny_sub, ix:ix + nx_sub] = sol
            gradient_u[iz:iz + nz_sub, iy:iy + ny_sub, ix:ix + nx_sub + 1] = grad_u
            gradient_v[iz:iz + nz_sub, iy:iy + ny_sub + 1, ix:ix + nx_sub] = grad_v
            gradient_w[iz:iz + nz_sub + 1, iy:iy + ny_sub, ix:ix + nx_sub] = grad_w
        except:
            pass

    dataset.close()

globmin = np.nanmin([np.nanmin(sol) for sol in solutions])
globmax = np.nanmax([np.nanmax(sol) for sol in solutions])

#globmin = np.min([np.nanmin(sol) for sol in partitions])
#globmax = np.max([np.nanmax(sol) for sol in partitions])

fig = plt.figure(figsize=[10, 10])
ax = fig.add_subplot()
ax.set_aspect(1.0)

#for i, sol in enumerate(partitions):
for i, sol in enumerate(solutions):
    mp = ax.pcolormesh(x2_arrs[i], y2_arrs[i], sol, vmin=globmin, vmax=globmax, cmap=cm.jet)#, rasterized=True)
    ax.add_patch(PolygonPatch(block_bnds[i], fc=[0,0,0,0], ec="black"))
    ax.text(block_bnds[i].centroid.coords[0][0], block_bnds[i].centroid.coords[0][1], "{:d}, {:d}".format(block_ids[i], block_ranks[i]), fontsize=5)

cax = plt.axes([0.92, 0.45, 0.03, 0.2])
plt.colorbar(mappable=mp, cax=cax)

plt.savefig("solution_z{:d}.eps".format(int(z_slice)))


dataset = Dataset("output.nc", 'w')
dataset.createDimension('x', x_glob.size)
dataset.createDimension('x2', x2_glob.size)
dataset.createDimension('y', y_glob.size)
dataset.createDimension('y2', y2_glob.size)
dataset.createDimension('z', z_glob.size)
dataset.createDimension('z2', z2_glob.size)

dataset.createVariable('x2', np.float64, 'x2')
dataset.createVariable('y2', np.float64, 'y2')
dataset.createVariable('z2', np.float64, 'z2')
dataset.createVariable('x', np.float64, 'x')
dataset.createVariable('y', np.float64, 'y')
dataset.createVariable('z', np.float64, 'z')
dataset.createVariable('solution', np.float64, ('z', 'y', 'x'))
dataset.createVariable('gradient_u', np.float64, ('z', 'y', 'x2'))
dataset.createVariable('gradient_v', np.float64, ('z', 'y2', 'x'))
dataset.createVariable('gradient_w', np.float64, ('z2', 'y', 'x'))

dataset['x2'][:] = x2_glob
dataset['y2'][:] = y2_glob
dataset['z2'][:] = z2_glob

dataset['x'][:] = x_glob
dataset['y'][:] = y_glob
dataset['z'][:] = z_glob

dataset['solution'][:] = solution
dataset['gradient_u'][:] = gradient_u
dataset['gradient_v'][:] = gradient_v
dataset['gradient_w'][:] = gradient_w

dataset.close()
