import xarray as xr
import numpy as np
import dask
from dask.distributed import Client
dask.config.config.get('distributed').get('dashboard').update({'link':'{JUPYTERHUB_SERVICE_PREFIX}proxy/{port}/status'})
import shutil
import os
import sys
import random
import string

#saving
def save_file(ds_out,outfile,overwrite):
    if (os.path.exists(outfile))&(overwrite==True):
        os.remove(outfile)
        ds_out.to_netcdf(outfile, compute=True)
        print(outfile+' overwritten')
    elif (os.path.exists(outfile))&(overwrite==False):
        print(outfile+' skipped')
    else:   
        ds_out.to_netcdf(outfile, compute=True)
        print(outfile+' done')    

def generate_random_string(length):
    # Define the characters to choose from
    characters = string.ascii_letters + string.digits  # You can add more characters if needed
    # Use random.choices to generate a list of random characters
    random_characters = random.choices(characters, k=length)
    # Join the characters into a string
    random_string = ''.join(random_characters)
    return random_string
        
def main():

    overwrite=True
    #get top and bottom from input
    top=int(sys.argv[1])
    bottom=int(sys.argv[2])

    n_cores = 5
    mem_lim = str(int(np.floor(95/n_cores)))+'GB'
    dask_dir = '/p/scratch/chhb19/mueller29/dask_dir/'+generate_random_string(10)
    if os.path.exists(dask_dir):
        shutil.rmtree(dask_dir)
    if 'client' in locals() or 'client' in globals():
        client.close()
    client = Client(local_directory=dask_dir,n_workers=n_cores, threads_per_worker=1,memory_limit=mem_lim)
    client.amm.start()

    
    #paths 
    data_path = '/p/scratch/chhb19/mueller29/AO_40_nopc/'
    out_path  = '/p/scratch/chhb19/mueller29/maps/'

    # chunk sizes
    horizontal_split = 4000000
    time_split = 1

    # depth indices for averaging

    z1 = xr.open_dataset(data_path+'unod.fesom.2015.nc').coords['nz1'] 
    if top==bottom:
        iz1 = xr.zeros_like(z1, dtype='bool')
        iz1[np.argmin(np.abs(z1.values-top))] = True
        iz1 = xr.DataArray (iz1)
    else:
        iz1 = (z1<=bottom)&(z1>=top)
    
    
    
    for yy in range(2015,2021):
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            u=xr.open_dataset((data_path+'unod_monthly.fesom.' + str(yy) + '.nc'),chunks={"time":time_split,"nod2":horizontal_split})['unod'].astype('float32')
            v=xr.open_dataset((data_path+'vnod_monthly.fesom.' + str(yy) + '.nc'),chunks={"time":time_split,"nod2":horizontal_split})['vnod'].astype('float32')
            #take average here, since uu/vv are already squared
            uu=xr.open_dataset((data_path+'uu.fesom.' + str(yy) + '.nc'),chunks={"time":time_split,"nod2":horizontal_split})['uu'].astype('float32')[:,iz1,:].weighted(z1[iz1]).mean(dim="nz1") 
            vv=xr.open_dataset((data_path+'vv.fesom.' + str(yy) + '.nc'),chunks={"time":time_split,"nod2":horizontal_split})['vv'].astype('float32')[:,iz1,:].weighted(z1[iz1]).mean(dim="nz1")
            #make same time coordinates
            u.coords['time']=uu.coords['time'].values
            v.coords['time']=vv.coords['time'].values

            MKE = 0.5*(u**2+v**2)[:,iz1,:].weighted(z1[iz1]).mean(dim="nz1") #take average AFTER squaring!!!!!
            TKE = 0.5*(uu+vv)
            EKE = TKE - MKE

            #make xarray dataset
            ds = xr.Dataset(
                data_vars=dict(
                    EKE=EKE,
                    TKE=TKE,
                    MKE=MKE,
                ),
                attrs=dict(description="KE components wrt monthly averages"),
            )

        #write dataset to disk

        outfile = (out_path+'KE.fesom.'+str(top)+'_'+str(bottom)+'_mean.'+str(yy)+'.nc')
        save_file(ds,outfile,overwrite)

        

    client.close()
    shutil.rmtree(dask_dir)

        
if __name__ == '__main__':
    main()