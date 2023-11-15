overwrite=True


import xarray as xr
import os
import sys

#get top and bottom from input
top=int(sys.argv[1])
bottom=int(sys.argv[2])

#paths 
data_path = '/p/scratch/chhb19/mueller29/AO_40_nopc/'
mesh_path = '/p/project/chhb19/mueller29/meshes/AO_40/'
out_path  = '/p/scratch/chhb19/mueller29/maps/'

# chunk sizes
horizontal_split = 5000000
vertical_split = 1
time_split = 1

# depth indices for averaging
z = xr.open_dataset(data_path+'w.fesom.2014.nc').coords['nz'] 
g = 9.81 #[m/s2]
rho0 = 1027 #[kg/m3]
iz = (z<=bottom)&(z>=top)

# stratification 3D monthly 
var='N2'
for yy in range(2014,2021):
    ds = xr.open_mfdataset(data_path+var+'.fesom.' + str(yy) + '.nc', chunks={'time': time_split,'nod2':horizontal_split}, concat_dim='time', combine='nested',data_vars='minimal', coords='minimal', compat='override', parallel=True)[var].astype('float32')[:,iz,:].weighted(z[iz]).mean(dim="nz") 


    #make xarray dataset
    ex_stri='ds_out = xr.Dataset(data_vars=dict('+var+'=ds,),attrs=dict(description="average '+var+'"),)'
    exec(ex_stri)

    outfile = (out_path+var+'.fesom.'+str(top)+'_'+str(bottom)+'_mean.'+str(yy)+'.nc')

    #write dataset to disk
    if (os.path.exists(outfile))&(overwrite==True):
        os.remove(outfile)
        ds_out.to_netcdf(outfile, compute=True)
        print(var+' overwritten')
    elif (os.path.exists(outfile))&(overwrite==False):
        print(var+' skipped')
    else:   
        ds_out.to_netcdf(outfile, compute=True)
        print(var+' done')
    
    
# temperature/salinity 3D monthly 
vars = ('temp','salt')
for var in vars:
    for yy in range(2014,2021):
        ds = xr.open_mfdataset(data_path+var+'.fesom.' + str(yy) + '.nc', chunks={'time': time_split,'nod2':horizontal_split}, concat_dim='time', combine='nested',data_vars='minimal', coords='minimal', compat='override', parallel=True)[var].astype('float32')[:,iz,:].weighted(z[iz]).mean(dim="nz") 
        #make xarray dataset
        ex_stri='ds_out = xr.Dataset(data_vars=dict('+var+'=ds,),attrs=dict(description="average '+var+'"),)'
        exec(ex_stri)
        #write dataset to disk
        outfile = (out_path+var+'.fesom.'+str(top)+'_'+str(bottom)+'_mean.'+str(yy)+'.nc')
        if (os.path.exists(outfile))&(overwrite==True):
            os.remove(outfile)
            ds_out.to_netcdf(outfile, compute=True)
            print(var+' overwritten')
        elif (os.path.exists(outfile))&(overwrite==False):
            print(var+' skipped')
        else:   
            ds_out.to_netcdf(outfile, compute=True)
            print(var+' done')
        
        
