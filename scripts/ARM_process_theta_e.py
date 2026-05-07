'arm process'
# constants
rd = 287.04 # gas constant of dry air
rv = 461.50 # gas constant of water vapor
cvd = 719. # heat capacity at constant volume for dry air
cpd = 1005.7 # heat capacity at constant pressure for dry air
cvv = 1410. # heat capacity at constant volume of water vapor
cpv = 1870. # heat capacity at constant pressure of water vapor
cl = 4190. # heat capacity of liquid water (above freezing)
cpi = 2093.
cpvmcl = cl-cpv # cpvmcl seems to be a common notation for this value

# Units: dimensionless
epsilon = rd/rv

# Units: J/Kg
lv0 = 2.501E6 # latent heat of vaporization at 0-deg-C
ls = 2.834E6 # latent heat of sublimation (-100<=T<= 0-deg-C)
lf = 0.3337E6 # latent heat of fusion at 0-deg-C (lv0-ls)

# Units: m/s**2
g = 9.80665 # standard gravity
# import packages
import os
import numpy as np
import netCDF4 as nc
import datetime as dt
import sys
from glob import glob
import matplotlib.pyplot as plt
import calendar
import numpy.ma as ma
import string
import matplotlib.cm as cmap
from scipy import integrate
from scipy import interpolate
from scipy import stats
from matplotlib.ticker import NullFormatter
import scipy
from scipy import signal
import xarray as xr
from tqdm import tqdm

from beam_model import beam_model_calc,eval_es,eval_theta_e,eval_q,eval_qs




f_dir = sorted(glob('/neelin2020/ARM/TWPC1/*.cdf'))+ sorted(glob('/neelin2020/ARM/TWPC2' + '/*.cdf'))
f_lwc_dir = sorted(glob('/neelin2020/ARM/TWPC1/ice/*.nc'))+ sorted(glob('/neelin2020/ARM/TWPC2' + '/ice/*.nc'))
site = 'Manus+ Nauru'
lev = nc.Dataset(f_dir[0],'r').variables['p'][:]*100



lft_t = np.argmin(lev>600e2)

uft_ind = np.argmin(lev>=600e2)
uft_t_ind = np.argmin(lev>=200e2)

lft_b = np.argmin(lev>900e2)
bl_t = np.argmin(lev>900e2)


# Calculate theta 
t_dim = 8760*17*2
t_count = 0
theta_es_lft = np.zeros((t_dim))*np.nan
theta_e_lft= np.zeros((t_dim))*np.nan
ta_lft= np.zeros((t_dim))*np.nan

time = np.zeros((t_dim),dtype=object)

theta_e_bl= np.zeros((t_dim))*np.nan
q_bl= np.zeros((t_dim))*np.nan
qs_bl= np.zeros((t_dim))*np.nan

q_lft= np.zeros((t_dim))*np.nan
qs_lft= np.zeros((t_dim))*np.nan

pr= np.zeros((t_dim))*np.nan
rh_full = np.zeros((t_dim,len(lev)))*np.nan
theta_e = np.zeros((t_dim,len(lev)))*np.nan
theta_es = np.zeros((t_dim,len(lev)))*np.nan
lwc = np.zeros((t_dim,len(lev)))*np.nan
lwc_mwrp = np.zeros((t_dim,len(lev)))*np.nan

ice = np.zeros((t_dim,len(lev)))*np.nan

omega = np.zeros((t_dim,len(lev)))*np.nan

cl_t= np.zeros((t_dim))*np.nan
cl_b= np.zeros((t_dim))*np.nan
cl_type= np.zeros((t_dim,8))*np.nan

q = np.zeros((t_dim,len(lev)))*np.nan
qs_full  = np.zeros((t_dim,len(lev)))*np.nan

ta_full  = np.zeros((t_dim,len(lev)))*np.nan

es_full  = np.zeros((t_dim,len(lev)))*np.nan

t_count= 0
for a,m in enumerate(f_dir):
    # if site_num<=2:
    if site_num==3:
        f = nc.Dataset(m,'r')
        ta = f.variables['temperature_p'][:].data
        rh = f.variables['relative_humidity_p'][:].data
        t = len(ta)
        pr[t_count:t_count + t] = f.variables['precip_rate_sfc'][:].data    
        lwc_pre = nc.Dataset(f_lwc_dir[a],'r').variables['lwc'][:]
        lwc_mwrp_pre = nc.Dataset(f_lwc_dir[a],'r').variables['lwc_mwrp'][:]

        ice_pre = nc.Dataset(f_lwc_dir[a],'r').variables['iwc'][:]

        cl_type_pre = nc.Dataset(f_cl_dir[a],'r').variables['cl_type'][:]
        cl_b_pre = nc.Dataset(f_cl_dir[a],'r').variables['cl_b'][:]
        cl_t_pre = nc.Dataset(f_cl_dir[a],'r').variables['cl_top'][:]

        omega_pre_nan = f.variables['omega_nwp_p'][:]
        masked_arr = np.ma.masked_array(omega_pre_nan, mask=omega_pre_nan.mask)
        
        # Fill masked values with NaN
        omega_pre = masked_arr.filled(np.nan)
        omega[t_count:t_count + t] = omega_pre
        lwc[t_count:t_count + t] = lwc_pre
        lwc_mwrp[t_count:t_count + t] = lwc_mwrp_pre

        ice[t_count:t_count + t] = ice_pre

        cl_b[t_count:t_count + t] = cl_b_pre
        cl_t[t_count:t_count + t] = cl_t_pre
        cl_type[t_count:t_count + t] = cl_type_pre

        time_var = f.variables['time']
    
        time_values = time_var[:]
        units = time_var.units
        calendar = time_var.calendar if hasattr(time_var, 'calendar') else 'standard'
        
        # Convert numerical time values to datetime objects
        datetime_values = nc.num2date(time_values, units=units, calendar=calendar)
        time[t_count:t_count+t] = datetime_values

    else: 
        f = nc.Dataset(m,'r')
        ta = f.variables['T_p'][:].data
        rh = f.variables['rh_p'][:].data
        t = len(ta)
        pr[t_count:t_count + t] = f.variables['prec_sfc'][:].data
        lwc_pre = nc.Dataset(f_lwc_dir[a],'r').variables['lwc'][:]
        # lwc_mwrp_pre = nc.Dataset(f_lwc_dir[a],'r').variables['lwc_mwrp'][:]

        ice_pre = nc.Dataset(f_lwc_dir[a],'r').variables['iwc'][:]
        time_var = f.variables['time']
        # tv = ta*(1+0.61*
        # rho = 
        lwc[t_count:t_count + t] = lwc_pre

        ice[t_count:t_count + t] = ice_pre
        time_values = time_var[:]
        units = time_var.units
        calendar = time_var.calendar if hasattr(time_var, 'calendar') else 'standard'
        
        # Convert numerical time values to datetime objects
        datetime_values = nc.num2date(time_values, units=units, calendar=calendar)
        time[t_count:t_count+t] = datetime_values

            
    for i in np.arange(t):
        if any(np.logical_or(ta[i] < 0, rh[i] < 0)):
            continue
        rh_full[t_count + i] = rh[i]
        last_work = i
        Es = eval_es(ta[i])
        # Es = es_calc(ta[i])
        e = Es*rh[i]/100
        es_full[t_count + i] = Es
        # w = e*epsilon/(lev-e)#epsilon*Es*rh[i]/100/lev
        q_temp = eval_q(e,lev)

        # q_temp = w/(w+1)
        q_temp[q_temp<0] = np.nan
        q[t_count+i,:] = q_temp
        qs = eval_qs(ta[i,uft_ind],lev[uft_ind])
        qs_col = eval_qs(ta[i,bl_t:uft_ind],lev[bl_t:uft_ind])
        
        theta_e_temp = eval_theta_e(ta[i],q_temp,0,0,lev)
        theta_e_bl[t_count + i] = np.nanmean(theta_e_temp[:bl_t+1])
        

        theta_e[t_count+i] = theta_e_temp
          
        theta_e_bl[t_count + i] = np.nanmean(theta_e_temp[:bl_t+1])
        theta_e_lft[t_count + i] = np.nanmean(theta_e_temp[lft_b:lft_t+1])#np.nanmean(theta_e_calc2d(ta[i,bl_t:uft_ind],q_temp[bl_t:uft_ind],lev[bl_t:uft_ind]))
        ta_lft[t_count + i] = np.nanmean(ta[i,lft_b:lft_t+1])
        q_bl[t_count + i] = np.nanmean(q_temp[:bl_t+1])
        q_lft[t_count + i] = np.nanmean(q_temp[lft_b:lft_t+1])
    #     qs_uft = qs_calc2d(ta[i,uft_t_ind:uft_ind,:,:],lev[uft_t_ind:uft_ind])
        qs_temp  = eval_qs(ta[i],lev)
        qs_bl[t_count + i] = np.nanmean(qs_temp[:bl_t + 1])
        qs_lft_full  = qs_temp[lft_b:lft_t+1]
    #     qs_bl = qs_calc2d(ta[i,bl_t:,:,:],lev[bl_t:])
        qs_full[t_count + i] = qs_temp
        ta_full[t_count + i] = ta[i]
    
    #     qs_lft_avg[i] = np.nanmean(qs_lft,axis=0)
#         theta_es_uft_avg[i] = np.nanmean(theta_e_calc2d(ta[i,uft_t_ind:uft_ind],qs_uft,lev[uft_t_ind:uft_ind]),axis=0)
        theta_es_temp = eval_theta_e(ta[i],qs_temp,0,0,lev)
        theta_es[t_count + i] = theta_es_temp
        theta_es_lft[t_count + i] = np.nanmean(theta_es_temp[lft_b:lft_t+1])
        qs_lft[t_count+i] = np.nanmean(qs_lft_full)#eval_qs(ta_lft[t_count+i],750e2)#np.nanmean(qs_lft_full)

        
        # Extract numerical time values, units, and calendar attributes

    #     theta_es_bl_avg[i] = np.nanmean(theta_e_calc2d(ta[i,bl_t:,:,:],qs_bl,lev[bl_t:]),axis=0)
    #     w_lft[i] = np.nanmean(w[i,uft_ind:bl_t,:,:],axis=0)
    #     w_uft[i] = np.nanmean(w[i,uft_t_ind:uft_ind,:,:],axis=0)
  
    t_count = t_count + t

theta_es_lft = theta_es_lft[:t_count]

theta_e_lft= theta_e_lft[:t_count] 
ta_lft= ta_lft[:t_count] 

theta_e_bl= theta_e_bl[:t_count] 
q_bl= q_bl[:t_count] 
q_lft= q_lft[:t_count] 
qs_lft= qs_lft[:t_count] 

pr= pr[:t_count] 
rh_full = rh_full[:t_count,:] 
theta_e = theta_e[:t_count,:] 
theta_es = theta_es[:t_count,:]
lwc = lwc[:t_count,:]
lwc_mwrp = lwc_mwrp[:t_count,:]

ice = ice[:t_count,:]

omega = omega[:t_count,:] 
cl_type = cl_type[:t_count]
cl_b = cl_b[:t_count]
cl_t = cl_t[:t_count]


q = q[:t_count,:] 
qs_full = qs_full[:t_count,:] 
ta_full = ta_full[:t_count,:] 
es_full = es_full[:t_count,:] 

ta_full[ta_full<0] = np.nan

time = time[:t_count]




tv = ta_full*(1+0.61*q)
rho = lev/rd/tv 
lwc_gkg = (lwc/rho).data
iwc_gkg = (ice/rho).data

qc_env = iwc_gkg + lwc_gkg






pr[pr<0] = 0
theta_e_bl[theta_e_bl <= 300] = np.nan
theta_e_lft[theta_e_lft <= 273] = np.nan
# theta_es_lft[theta_es_lft <= 273] = np.nan
ta_lft[ta_lft<0] = np.nan
# theta_e[theta_e <= 273] = np.nan
# theta_es[theta_es <= 273] = np.nan
q_lft[q_lft<0] = np.nan
q_bl[q_bl<0] = np.nan

bad_inds = np.nanmean(theta_e[:,:bl_t+1],axis=1)<300
theta_e[bad_inds] = np.ones((np.sum(bad_inds),len(lev)))*np.nan
theta_es[bad_inds] = np.ones((np.sum(bad_inds),len(lev)))*np.nan
q[bad_inds]= np.ones((np.sum(bad_inds),len(lev)))*np.nan
qs_full[bad_inds]= np.ones((np.sum(bad_inds),len(lev)))*np.nan



