add_library(mpp STATIC
   mpp.c
   mpp_domain.c
   mpp_io.c
   tool_util.c
)

target_link_libraries(mpp PUBLIC shared)
target_include_directories(mpp INTERFACE ${CMAKE_CURRENT_LIST_DIR})