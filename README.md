### Build instruction

```
mkdir build
cd build
cmake /path/to/gridspec_tools
make -j
```

### Instructions for generating C48 gridpec file
```
PATH=$PATH:/path/to/build/dir/bin
make_hgrid --grid_type gnomonic_ed --nlon 48
make_mosaic --num_tiles 6 --tile_file $(pwd)/horizontal_grid --generate_contact  # horizontal_grid is prefix of tile file names
python update_gridspec.py mosaic.nc mosaic2.nc # mosaic2 is the fixed one
```
