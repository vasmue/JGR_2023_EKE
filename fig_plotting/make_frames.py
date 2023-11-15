import pyfesom2 as pf
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.tri as mtri
from matplotlib.ticker import LogFormatter

import cartopy.crs as ccrs

import cmocean
import sys

# input year
yy = int(sys.argv[1])
# input depth 
depth = int(sys.argv[2])

# chunk sizes
horizontal_split = 5000000
vertical_split = 1
time_split = 1

#paths 
data_path = '/p/scratch/chhb19/mueller29/AO_40_nopc/'
mesh_path = '/p/project/chhb19/mueller29/meshes/AO_40/'
plot_path  = '/p/project/chhb19/mueller29/animation_frames/'

#the right depth level
z1 = xr.open_dataset(data_path+'unod.fesom.2015.nc').coords['nz1'] 
iz1 = np.argmin(np.abs((z1-depth).values))

#define euler angles 
alpha=-90
beta=90 
gamma=90

#region
left = -180
right = 180
south = 70
north = 90

#load mesh and rotate to equator
mesh = pf.load_mesh(mesh_path)
elements=mesh.elem.astype('int32')
lons_rot,lats_rot = pf.ut.scalar_g2r(alpha,beta,gamma,mesh.x2,mesh.y2)
d = lons_rot[elements].max(axis=1) - lons_rot[elements].min(axis=1)
no_cyclic_elem = np.argwhere(d < 100).ravel()


#regular grid for interpolation at equator
dx=0.01
dy=0.01
nx2=north-south+dx #just to make sure the grid isn't square (makes it easier to keep track of dimensions)
ny2=north-south
lon_reg = np.arange(-nx2, nx2, dx) 
lat_reg = np.arange(-ny2, ny2, dy)

xx_eq, yy_eq = np.meshgrid(lon_reg, lat_reg)
xx_eq=xx_eq.T;
yy_eq=yy_eq.T;
xx_pol,yy_pol = pf.ut.scalar_r2g(alpha,beta,gamma,xx_eq,yy_eq) #regular grid from equator rotated to North (for loading mask)
#linear interpolator
triang = mtri.Triangulation(lons_rot, lats_rot, elements[no_cyclic_elem])
tri = triang.get_trifinder()


#for ice interpolate without rotation with lower resolution
#lons = mesh.x2
#lats = mesh.y2
#d = lons[elements].max(axis=1) - lons[elements].min(axis=1)
#no_cyclic_elem = np.argwhere(d < 100).ravel()
#triang_ice = mtri.Triangulation(lons, lats, elements[no_cyclic_elem])
#tri_ice = triang_ice.get_trifinder()



#load all u/v data
fn_u = data_path+'unod.fesom.'+str(yy)+'.nc'
fn_v = data_path+'vnod.fesom.'+str(yy)+'.nc'

ds_u = xr.open_dataset(fn_u, chunks={'time': time_split,'nz1':vertical_split})['unod'].astype('float32')[:,iz1,:]
ds_v = xr.open_dataset(fn_v, chunks={'time': time_split,'nz1':vertical_split})['vnod'].astype('float32')[:,iz1,:]

plt.rcParams['figure.dpi']=288
for dd in range(len(ds_u.time.values)):
    #calculate speed
    speed = np.sqrt(ds_u[dd,:]**2+ds_v[dd,:]**2).drop_vars('nz1')
    #interpolate for plot
    speed_grid = mtri.LinearTriInterpolator(triang,speed,trifinder=tri)(xx_eq, yy_eq)
    speed_grid[np.where(speed_grid==0)] = np.nan
    

    frame_n=plot_path+str(depth)+'/'+'speed_'+str(depth).zfill(3)+'m'+str(yy)+'_'+str(dd+1).zfill(3)+'.png'
    figsize = (10,10) #1/4 page
    fig, ax = plt.subplots(
            1,
            1,
            subplot_kw=dict(projection=ccrs.RotatedPole(pole_longitude=240.0, pole_latitude=45.0, central_rotated_longitude=0.0)),
            constrained_layout=True,
            figsize=figsize,
        )

    ax.set_extent([-50,170,77,90], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m')

    image = ax.pcolormesh(xx_pol, yy_pol, speed_grid, cmap = cmocean.cm.ice, norm=colors.LogNorm(vmin=1e-2, vmax=1), shading='auto', transform=ccrs.PlateCarree())
    ax.text(0.01,0.93,str(yy)+'_'+str(dd+1).zfill(3), size=15,transform=ax.transAxes,zorder=100)

    
    formatter = LogFormatter(10, labelOnlyBase=False) 
    cb = fig.colorbar(image, ticks=[0.01,0.1,1], format=formatter, orientation='horizontal', ax=ax, pad=0.02, shrink=0.85, extend='both')
    cb.ax.tick_params(labelsize=15)
    cb.ax.set_xticklabels(['0.01', '0.1', '1'])
    cb.set_label('speed [m/s]', size=15)
    fig.set_size_inches((10,4.25), forward=True)
    fig.savefig(frame_n,format='png')
    plt.close(fig)
    
    
    
print(str(yy)+' all frames done')