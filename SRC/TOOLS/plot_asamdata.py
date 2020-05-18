#**********#
# 1. Header
#**********#
import sci_meshgrid
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import os
import types

import Scientific.IO.NetCDF as netcdf
# import Numeric
# import NetCDF

import re

#----- Import my own package
import ptpack.func as ppf
import ptpack.parameter as ppp
import ptpack.datapack as ppd


# varmeteo = ppd.VarInfoMeteo('u', 'u', 'ms', 'ms', 'uuu')
varlib = ppd.VarLib()
print varlib.get_disname('diffmom')
print varlib
# sys.exit()
#**********#
# 2. Read data
#**********#
# asam_dir = '/media/My Passport/Data/ASAM_P1/CanopyTest/'
asam_dir = '/home/pzzhou/Scripts/Fortran/ASAM_SYNC/ParExample/Canopy_rans/'
# case_name = 'case_20140528_1-1-42_c-c+ce-s+_08-72'
# case_name = 'case_20140530_1-1-42_c-c-ce-s+tkedis_08-72'
# case_name = 'case_20140522_144-96-42_nosoil_noslip'
# case_name = 'case_20140605_1-1-42_00-48_rans-lng0-meteo'
# case_name = 'case_20140606_1-1-42_00-48_rans-lng0-radprof-meteo'
# case_name = 'case_20140606_1-1-42_08-48_0-rp-met'
# case_name = 'case_20140610_1-1-42_08-72_rp-met'
# case_name = 'case_20140611_1-1-66_08-72_rp-met'
# case_name = 'case_20140611_1-1-66_l4-rp-met'
# case_name = 'case_20140611_1-1-100_08-72_rp-met_l4'
case_name = 'case_20140613_1-1-66_rp-met_wrf-damp400-diff1'
# asam_case_dir = 'case_20140522_144-96-42_nosoil_noslip/'
# asam_case_dir = 'case_20140526_1-1-42_nochem_nocnpy_nosoil_72/'
asam_case_dir = case_name + '/'
asam_dir = asam_dir + asam_case_dir
asam_ncfile = 'Canopy_rans_1-1-66_rp-met.nc'
asam_gdfile = 'Canopy_rans_1-1-66_rp-met.grid'
asam_dp = ppd.AsamDataPack(asam_dir, asam_ncfile, asam_gdfile)

save_dir = './' + case_name + '/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
##### variables' list, in the order of:
##### var_name: (0-name in asam, 1-display name, 2-unit, 3-min, 4-max)
##### ASAM data format: data[time, lev, yc, xc]
##### var_prop: name in nc file, display name, display unit, min, max
var_dict = {
  'u': ('U', 'u', 'm/s', -4.0, 6.0),
  'v': ('V', 'v', 'm/s', -2.0, 4.0),
  'w': ('W', 'w', 'm/s', -1.0, 1.0),
  'theta': ('THETA', r'$\theta$', 'K', 290.0, 315.0),
  # 'thetarho': ('ThDens', r'$\theta_\rho$', 'K', 290.0, 315.0),
  'rho': ('RHO', r'$\rho$', r'kg m$^{-3}$', 0.8, 1.3),
  'tke': ('TKE', 'tke', r'm$^2$ s$^{-2}$', 0.0, 0.15),
  'qv': ('QV', r'$q_v$', r'kg m$^{-3}$', 0.0, 0.03),
  'pre': ('PRE', 'pressure', 'hPa', 70000.0, 105000.0),
  # 'diffpot': ('DIFFPOT', 'diffpot', r'm$^2$ s$^{-1}$', 0.0, 10.0),
  # 'diffmom': ('DIFFMOM', 'diffmom', r'm$^2$ s$^{-1}$', 0.0, 5.0),
  'abstemp': ('AbsTemp', 'T', 'K', 270.0, 312.0),
  # 'iso': ('ISO', 'ISO', r'molec cm$^{-3}$', 0.0, 5.0e9),
  # 'api': ('API', r'$\alpha$-pinene', r'molec cm$^{-3}$', 0.0, 1.0e10),
  # 'bpi': ('BPI', r'$\beta$-pinene', r'molec cm$^{-3}$', 0.0, 6.0e9),
  # 'limn': ('LIMN', 'LIMN', r'molec cm$^{-3}$', 0.0, 2.0e10),
  # 'car3': ('CAR3', 'CAR3', r'molec cm$^{-3}$', 0.0, 2.0e8),
  # 'oh': ('HO', 'OH', r'molec cm$^{-3}$', 0.0, 2.0e7),
  # 'ho2': ('HO2', r'HO$_2$', r'molec cm$^{-3}$', 0.0, 1.0e9),
  # 'o3': ('O3', r'O$_3$', r'molec cm$^{-3}$', 4.0e11, 5.5e11),
  # 'so2': ('SO2', r'SO$_2$', r'molec cm$^{-3}$', 0.8e10, 1.30e10),
  # 'no': ('NO', 'NO', r'molec cm$^{-3}$', 0.0, 1.0e9),
  # 'no2': ('NO2', r'NO$_2$', r'molec cm$^{-3}$', 0.0, 6.0e9),
  # 'no3': ('NO3', r'NO$_3$', r'molec cm$^{-3}$', 0.0, 1.0e6),
  # 'co': ('CO', 'CO', r'molec cm$^{-3}$', 2.12e12, 2.50e12),
  # 'bcar': ('BCAR', 'BCAR', r'molec cm$^{-3}$', 0.0, 3.0e8),
  # 'farn': ('FARN', 'FARN', r'molec cm$^{-3}$', 0.0, 5.0e7),
  }

xx = asam_dp.gd.xx
yy = asam_dp.gd.yy
zz = asam_dp.gd.zz
# print zz
#**********#
# 3. Plot figures
#**********#
plt.rc('font', size=18)
#===== Time array for the data =====#
time_array = np.linspace(asam_dp.gd.t0, asam_dp.gd.t1, (asam_dp.gd.t1 - asam_dp.gd.t0)/asam_dp.gd.odt + 1) # asam_dp.get_data('time')

# print asam_dp.gd.odt
# print time_array
# print zz
# sys.exit()

#===== Plot profiles in one axis =====#
def plot_profiles_in_one(time_interv):
  '''
  The profiles at all time points are in one axis
  '''
  ind_array = np.arange(0, time_array.size)
  ind_plot = ind_array[20:70:8] # time point every two hours
  # ind_plot = ind_array[4:48:8] # time point every two hours
  num_plot = ind_plot.size
  
  cm = plt.get_cmap('gist_rainbow')
  # linecolor = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f'] # lighter colors
  # linecolor = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d'] # darker colors
  # linecolor = ['#b2182b', '#d6604d', '#f4a582', '#fddbc7', '#92c5de', '#4393c3', '#2166ac']
  linecolor = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#999999', '#a65628']
  linestyle = ['-', '--', '-.', ':']
  num_linestyle = len(linestyle)
  # print cm(0)
  # sys.exit()
  for var, var_prop in var_dict.iteritems():
    print var
  
    fg = plt.figure(figsize=(12,10), dpi=72)
    ax = fg.add_subplot(1, 1, 1)
    ax.set_color_cycle([cm(1.*i/num_plot) for i in range(num_plot)])
  
    data = asam_dp.get_data(var_prop[0])
    data_h_mean = np.mean(np.mean(data, axis=3), axis=2)
    for i in range(num_plot):
      color = cm(1.*i/num_plot)
      ax.plot(data_h_mean[ind_plot[i], :], zz,
        linestyle[i % num_linestyle],
        color = linecolor[i],
        linewidth = 2,
        label='{0:02.0f}:00'.format(time_array[ind_plot[i]]/3600. % 24))
    
    # ax.set_xlim(var_prop[3], var_prop[4])
    ax.set_title(var_prop[1])
    ax.set_xlabel(var_prop[1] + ' [' + var_prop[2] + ']')
    ax.set_ylabel('Height [m]')
    ax.legend()
  
    fg.savefig(save_dir + 'profile_' + var + '.png')

#===== Plot profiles as contour =====#
def plot_profiles_as_contour():
  for var, var_prop in var_dict.iteritems():
    print var

    fg = plt.figure(figsize=(12,10), dpi=72)
    ax = fg.add_subplot(1,1,1)
    
    data = asam_dp.get_data(var_prop[0])
    data_h_mean = np.mean(np.mean(data, axis=3), axis=2) # horizontal mean
    # print data_h_mean.shape
    # sys.exit()
    
    # levels = np.linspace(var_prop[3], var_prop[4], 10)
    cax = ax.contourf(time_array/3600./24, zz, data_h_mean.T)
    # cax.set_clim(vmin=var_prop[3], vmax=var_prop[4])
    # cax.cmap.set_under('yellow')
    # cax.cmap.set_over('cyan')

    # ax.set_xlim()
    ax.set_ylim(0, 2000)
    ax.set_title(var_prop[1])
    fg.colorbar(cax)

    #===== chemistry, canopy, soil,
    fg.savefig(save_dir + 'time_height_' + var + '.png')

plot_profiles_as_contour()
plot_profiles_in_one(1)
# plt.show()
'''
fg = plt.figure(figsize=(12,10), dpi=72)
ax = fg.add_subplot(1,1,1)

u = asam_dp.get_data('U')
v = asam_dp.get_data('V')
wind = np.sqrt(u*u+v*v)

wind_h_mean = np.mean(np.mean(wind, axis=3), axis=2) # horizontal mean
# print data_h_mean.shape
# sys.exit()

# levels = np.linspace(var_prop[3], var_prop[4], 10)
cax = ax.contourf(time_array/3600., zz, wind_h_mean.T)
# cax.set_clim(vmin=var_prop[3], vmax=var_prop[4])
# cax.cmap.set_under('yellow')
# cax.cmap.set_over('cyan')

# ax.set_xlim()
ax.set_ylim(0, 2000)
ax.set_title('wind')
fg.colorbar(cax)

fg.savefig(save_dir + 'asam_rans_time_height_' + 'wind' + '.png')
'''
