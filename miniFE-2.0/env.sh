find `pwd` -name "*.hpp" > sstmac_headers
find `pwd` -name "*.h" >> sstmac_headers
export SSTMAC_HEADERS=`pwd`/sstmac_headers
export SSTMAC_SKELETONIZE=1
