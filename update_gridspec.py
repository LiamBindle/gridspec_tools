import os.path

import xarray as xr

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generates a stretched grid initial restart file for GCHP.')
    parser.add_argument('input_file',
                        type=str,
                        help='path to the input file')

    parser.add_argument('output_file',
                        type=str,
                        help='output filename')
    args = parser.parse_args()
    ds = xr.open_dataset(args.input_file)
    # Get prefix
    prefix = os.path.commonprefix(ds.tile_files.astype(str).values.tolist())
    if not os.path.isdir(prefix):
        prefix = os.path.dirname(prefix)
    ds = ds.merge(xr.open_dataset(f'{prefix}/{ds["contact_files"].astype(str).item()}'))
    ds = ds.rename({'tile_files': 'gridfiles'})
    ds['gridtiles'] = xr.DataArray([f'tile{i+1}' for i in range(len(ds.ntiles))], dims=['ntiles'])
    ds['gridlocation'] = prefix
    ds['gridfiles'] = xr.DataArray([os.path.relpath(tilefile, prefix) for tilefile in ds.gridfiles.astype(str).values.tolist()], dims=['ntiles'])

    string_vars = ['gridfiles', 'gridtiles', 'gridlocation']
    ds.to_netcdf(args.output_file, encoding={sv: {'dtype':'S255'} for sv in string_vars})
    print('Done')


