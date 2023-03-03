  cmake -D Trilinos_ENABLE_OPTIONAL_PACKAGES:BOOL=ON \
  -D CMAKE_CXX_FLAGS:STRING="-O3" \
  -D CMAKE_C_FLAGS:STRING="-O3" \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D CMAKE_INSTALL_PREFIX:PATH=/usr/local/ \
  -D deal.II_DIR=/usr/local/dealii-9.4.0/lib/cmake \
  -D deal2lkit_DIR=/usr/local/deal2lkit_master_20230124 \
  -D CMAKE_CXX_STANDARD:="17" ../
