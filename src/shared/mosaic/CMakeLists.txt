add_library(shared STATIC
   create_xgrid.c
   gradient_c2l.c
   interp.c
   mosaic_util.c
   read_mosaic.c
)

target_link_libraries(shared PUBLIC NetCDF-C)
target_compile_definitions(shared PUBLIC use_netCDF)
target_include_directories(shared INTERFACE ${CMAKE_CURRENT_LIST_DIR})