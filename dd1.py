import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime

def era5_nc_to_csv(nc_file, csv_file, lat_point=22.0, lon_point=81.0):
    """Convert ERA5 NetCDF to CSV format for specific location"""
    
    # Open NetCDF file
    ds = xr.open_dataset(nc_file)
    
    print("Dataset info:")
    print(ds)  # Show dataset summary
    print("\nVariables:", list(ds.data_vars))
    print("Coordinates:", list(ds.coords))
    
    # Select specific location (nearest point)
    if 'latitude' in ds.coords and 'longitude' in ds.coords:
        ds_point = ds.sel(latitude=lat_point, longitude=lon_point, method='nearest')
    elif 'lat' in ds.coords and 'lon' in ds.coords:
        ds_point = ds.sel(lat=lat_point, lon=lon_point, method='nearest')
    else:
        print("Could not find lat/lon coordinates")
        ds_point = ds
    
    # Convert to DataFrame
    df = ds_point.to_dataframe().reset_index()
    
    # Clean up time column if present
    if 'time' in df.columns:
        df['Date'] = pd.to_datetime(df['time']).dt.date
        df['Time'] = pd.to_datetime(df['time']).dt.time
    
    # Rename common ERA5 variables
    column_mapping = {
        't': 'Temperature',
        'z': 'Geopotential',
        'q': 'Specific_humidity'
    }
    
    df = df.rename(columns=column_mapping)
    
    # Save to CSV
    df.to_csv(csv_file, index=False)
    
    print(f"Converted {nc_file} to {csv_file}")
    print(f"Final shape: {df.shape}")
    print(f"Final columns: {list(df.columns)}")
    
    return df

if __name__ == "__main__":
    df = era5_nc_to_csv('qvalue.nc', 'qvalue.csv')
