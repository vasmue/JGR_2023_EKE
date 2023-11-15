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
    z = xr.open_dataset(data_path+'w.fesom.2015.nc').coords['nz'] 
    g = 9.81 #[m/s2]
    rho0 = 1027 #[kg/m3]
    iz = (z<=bottom)&(z>=top)

    for yy in range(2015,2021):
        with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            w=xr.open_dataset((data_path+'w.fesom.' + str(yy) + '.nc'),chunks={"time":time_split,"nod2":horizontal_split})['w'].astype('float32')
            rhof=xr.open_dataset((data_path+'rhof.fesom.' + str(yy) + '.nc'),chunks={"time":time_split,"nod2":horizontal_split})['rhof'].astype('float32')
            #take average here, since wrhof is already product
            wrhof=xr.open_dataset((data_path+'wrhof.fesom.' + str(yy) + '.nc'),chunks={"time":time_split,"nod2":horizontal_split})['wrhof'].astype('float32')[:,iz,:].weighted(z[iz]).mean(dim="nz")



            #calculate MCR/TCR/ECR (CR = Conversion Rate)
            MCR = (-g/rho0)*w*rhof
            MCR = MCR[:,iz,:].weighted(z[iz]).mean(dim="nz") #take average AFTER product!!!!!

            TCR = (-g/rho0)*wrhof
            ECR = TCR - MCR


            #make xarray dataset
            ds = xr.Dataset(
                data_vars=dict(
                    ECR=ECR,
                    TCR=TCR,
                    MCR=MCR,
                ),
                attrs=dict(description="Energy conversion rates wrt monthly averages"),
            )


        #write dataset to disk

        outfile = (out_path+'wpbp.fesom.'+str(top)+'_'+str(bottom)+'_mean.'+str(yy)+'.nc')
        save_file(ds,outfile,overwrite)


    client.close()
    shutil.rmtree(dask_dir)

if __name__ == '__main__':
    main()