# IMERG directory
imergf = '/neelin2020/IMERG_V06/'

# COSMIC directory to match
cosmicf = nc.Dataset('/home/twemmen/COSMIC2_interp25hpa.nc','r')


cosmic_lat = cosmicf.variables['lat'][:]
cosmic_lon = cosmicf.variables['lon'][:]
t = cosmicf.variables['time']
cosmic_t = nc.num2date(t,units = t.units) 
cosmic_yrs = np.asarray([cosmic_t[x].year for x in np.arange(len(cosmic_t))])
cosmic_mons = np.asarray([cosmic_t[x].month for x in np.arange(len(cosmic_t))])
yrs = np.unique(cosmic_yrs)
mons  = np.unique(cosmic_mons)
pre = 6
post = 6
res = 2

import cftime
import re
date_format = '%Y%m%d-S%H%M%S'

count = 0
imlist = []
for a,y in enumerate(yrs):
    for b,m in enumerate(mons):
        mstr = str(m)
        if m < 10:
            mstr = '0' + str(m)
        imlist = imlist + glob(imergf + str(y) + '/' + mstr + '/*.HDF5')

filebeg = re.search('3B',imlist[0]).start()
dd = []
for c,i in enumerate(imlist):
    imdays = re.search('3IMERG.',i[filebeg:]).end()
    datestr = i[filebeg:][imdays:imdays+16]
    dd.append(datetime.strptime(datestr, date_format))        

    
f = h5py.File(imlist[0],'r')['Grid']
# qc = f['precipitationQualityIndex'][:]
lon = f['lon'][:].flatten()-.05
np.append(lon,f['lon'][-1]+.05)
lat = f['lat'][:].flatten()-.05
np.append(lat,f['lat'][-1]+.05)
t_data_ts = [cftime.date2num(t,'seconds since 1970-01-01') for t in cosmic_t]
t_edges_ts = [cftime.date2num(t,'seconds since 1970-01-01') for t in dd]
#### because there is strange gap
t_edges_ts = sorted(t_edges_ts)
###############
cosmic_t_arr = np.zeros((len(cosmic_t),(pre+post)*res))
cosmic_lat_arr = np.zeros((len(cosmic_t),(pre+post)*res))
cosmic_lon_arr = np.zeros((len(cosmic_t),(pre+post)*res))

from datetime import timedelta
for i in np.arange(-pre*res,post*res+1):
    cosmic_t_arr[:,i] = np.asarray(t_data_ts) + np.sign(i)*timedelta(minutes=60/res*i).seconds
    cosmic_lat_arr[:,i] = np.asarray(cosmic_lat) 
    cosmic_lon_arr[:,i] = np.asarray(cosmic_lon) 
 
t_inds = np.digitize(cosmic_t_arr,bins = t_edges_ts)-1
lon_inds= np.digitize(cosmic_lon_arr,bins = lon)-1
lat_inds = np.digitize(cosmic_lat_arr,bins = lat )-1


gap_inds = np.reshape(np.in1d(t_inds.flatten(),np.where([(np.diff(t_edges_ts)>1800)])[1]),t_inds.shape)
t_inds[gap_inds] = -1
t_inds[t_inds>=(len(t_edges_ts)-1)] = -1

pr_arr = np.zeros((len(cosmic_yrs),(pre + post)*res))*np.nan
start = t_inds[t_inds>0].min()
end = t_inds.max() 
for i in np.arange(t_inds[t_inds>0].min(),t_inds.max()):
    pr_imerg = h5py.File(imlist[i],'r')['Grid']['precipitationCal'][:].squeeze()
    inds = t_inds==i
    pr_arr[inds] = pr_imerg[lon_inds[inds],lat_inds[inds]]
    sys.stdout.write('\r Progress # {}/{}'.format(i-start,end-start))
    sys.stdout.flush()



t_vec = np.arange(-pre,post,.5)
cfmip_netcdf = nc.Dataset('/home/twemmen/COSMIC2_IMERG_match_021325.nc','w',format='NETCDF4')
n = len(t_inds)
cfmip_netcdf.createDimension('profiles',size=n)
t_val = 'seconds since 1970-01-01'

time_val = cfmip_netcdf.createVariable('date',np.float64,('profiles',))

time_val.units = t.units
time_val[:] =  t_data_ts


cfmip_netcdf.createDimension('time_rel',size=len(t_vec))
time_rel = cfmip_netcdf.createVariable('hour',np.float64,('time_rel',))
time_rel.units = 'hours'
time_rel[:] = t_vec

lon_val= cfmip_netcdf.createVariable('lon',np.float32,('profiles','time_rel'))
lon_val.units="degree"
lon_val[:]=lon[lon_inds]

lat_val= cfmip_netcdf.createVariable('lat',np.float32,('profiles','time_rel'))
lat_val.units="degree"
lat_val[:]= lat[lat_inds]

pr_series = cfmip_netcdf.createVariable('pr_prof',np.float32,('profiles','time_rel'))
pr_series.units="mm/hr"
pr_series[:]= pr_arr

cfmip_netcdf.close()  
