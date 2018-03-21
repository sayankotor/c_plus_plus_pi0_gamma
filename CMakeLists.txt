################################################################################
# Package: CaloTools
################################################################################
gaudi_subdir(CaloTools v6r13)

gaudi_depends_on_subdirs(Calo/CaloInterfaces
                         Calo/CaloUtils
                         Det/CaloDet
                         Event/LinkerEvent
                         Event/RecEvent
                         Event/TrackEvent
                         GaudiAlg
                         Kernel/LHCbKernel
                         Kernel/LHCbMath
                         Kernel/Relations
                         Tr/TrackInterfaces)

find_package(ROOT)
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
include_directories(/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_88/xgboost/0.60/x86_64-slc6-gcc49-opt/)
include_directories(/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_88/xgboost/0.60/x86_64-slc6-gcc49-opt/include)

gaudi_add_module(CaloTools
                 src/*.cpp
                 INCLUDE_DIRS Tr/TrackInterfaces
                 LINK_LIBRARIES CaloUtils CaloDetLib LinkerEvent RecEvent TrackEvent GaudiAlgLib LHCbKernel LHCbMathLib RelationsLib)

message("hello world!")

set(CMAKE_PREFIX_PATH /cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_88/xgboost/0.60/x86_64-slc6-gcc49-opt/lib)
include_directories(/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_88/xgboost/0.60/x86_64-slc6-gcc49-opt/rabit/include/)
message(${CMAKE_PREFIX_PATH})
message(${CMAKE_INCLUDE_PATH})
find_library(XGB_LIB22 xgboost)
message("xgb find here:")
message(${XGB_LIB22})
#message("include dirs:")
#message(${LIBXGBOOST_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(CaloTools ${XGB_LIB22})


