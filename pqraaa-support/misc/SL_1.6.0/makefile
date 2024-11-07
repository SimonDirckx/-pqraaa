CXXFLAGS =      -fno-exceptions -pedantic -fprefetch-loop-arrays -fopenmp  -I/Users/karlm/soft/SL_1.6.0/
LDFLAGS =   
LIBS =  -L/Users/karlm/soft/SL_1.6.0 -lsl_1_6 -lstdc++

all:
	ar -r libSL_1_6.a   /Users/karlm/soft/glas/SL_1.6.0/CuHyb.o /Users/karlm/soft/glas/SL_1.6.0/FileToVT.o /Users/karlm/soft/glas/SL_1.6.0/MachineInfo.o /Users/karlm/soft/glas/SL_1.6.0/Matrix2HilbertCoordinates.o /Users/karlm/soft/glas/SL_1.6.0/McDMV.o /Users/karlm/soft/glas/SL_1.6.0/McDirect.o /Users/karlm/soft/glas/SL_1.6.0/McShared.o /Users/karlm/soft/glas/SL_1.6.0/SBDTree.o /Users/karlm/soft/glas/SL_1.6.0/driver.o /Users/karlm/soft/glas/SL_1.6.0/micbench.o /Users/karlm/soft/glas/SL_1.6.0/mm2htrp.o /Users/karlm/soft/glas/SL_1.6.0/mmio.o /Users/karlm/soft/glas/SL_1.6.0/mtx2crstrp.o /Users/karlm/soft/glas/SL_1.6.0/sbd2trp-main.o /Users/karlm/soft/glas/SL_1.6.0/sbdmv.o /Users/karlm/soft/glas/SL_1.6.0/util.o

CUHYB_DEPS = \
	./CuHyb.hpp\
	./SparseMatrix.hpp\
	./ULIdef.hpp\
	./Matrix.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./Triplet.hpp\
	./alignment.hpp

CuHyb.o : CuHyb.cpp $(CUHYB_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o CuHyb.o CuHyb.cpp



FILETOVT_DEPS = \
	./FileToVT.hpp\
	./mmio.h\
	./Triplet.hpp\
	./ULIdef.hpp

FileToVT.o : FileToVT.cpp $(FILETOVT_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o FileToVT.o FileToVT.cpp



MACHINEINFO_DEPS = \
	./MachineInfo.hpp

MachineInfo.o : MachineInfo.cpp $(MACHINEINFO_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o MachineInfo.o MachineInfo.cpp



MATRIX2HILBERTCOORDINATES_DEPS = \
	./Matrix2HilbertCoordinates.hpp

Matrix2HilbertCoordinates.o : Matrix2HilbertCoordinates.cpp $(MATRIX2HILBERTCOORDINATES_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o Matrix2HilbertCoordinates.o Matrix2HilbertCoordinates.cpp



MCDMV_DEPS = \
	./Matrix.hpp\
	./Triplet.hpp\
	./ULIdef.hpp\
	./Upscaler.hpp\
	./SBDTree.hpp\
	./MinCRS.hpp\
	./BlockOrderer.hpp\
	./CBICRS.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./BICRS.hpp\
	./SparseMatrix.hpp\
	./alignment.hpp\
	./CCSWrapper.hpp\
	./HBICRS.hpp\
	./TS.hpp\
	./CRS.hpp\
	./ICRS.hpp\
	./ZZ_CRS.hpp\
	./ZZ_ICRS.hpp\
	./SVM.hpp\
	./HTS.hpp\
	./HilbertTriplet.hpp\
	./Matrix2HilbertCoordinates.hpp\
	./HilbertTripletCompare.hpp\
	./BetaHilbert.hpp\
	./FBICRS.hpp\
	./MachineInfo.hpp\
	./RDBHilbert.hpp\
	./BigInt.hpp

McDMV.o : McDMV.cpp $(MCDMV_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o McDMV.o McDMV.cpp



MCDIRECT_DEPS = \
	./Triplet.hpp\
	./ULIdef.hpp\
	./ICRS.hpp\
	./SparseMatrix.hpp\
	./Matrix.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./alignment.hpp\
	./CCSWrapper.hpp

McDirect.o : McDirect.cpp $(MCDIRECT_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o McDirect.o McDirect.cpp



MCSHARED_DEPS = \
	./Triplet.hpp\
	./ULIdef.hpp\
	./CRS.hpp\
	./SparseMatrix.hpp\
	./Matrix.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./alignment.hpp\
	./ICRS.hpp\
	./Hilbert.hpp\
	./HilbertTriplet.hpp\
	./Matrix2HilbertCoordinates.hpp\
	./HilbertTripletCompare.hpp\
	./CBICRS.hpp\
	./BICRS.hpp\
	./CCSWrapper.hpp

McShared.o : McShared.cpp $(MCSHARED_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o McShared.o McShared.cpp



SBDTREE_DEPS = \
	./SBDTree.hpp

SBDTree.o : SBDTree.cpp $(SBDTREE_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o SBDTree.o SBDTree.cpp



DRIVER_DEPS = \
	./SparseMatrix.hpp\
	./ULIdef.hpp\
	./Matrix.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./Triplet.hpp\
	./alignment.hpp\
	./TS.hpp\
	./CRS.hpp\
	./ICRS.hpp\
	./vecBICRS.hpp\
	./ZZ_CRS.hpp\
	./ZZ_ICRS.hpp\
	./SVM.hpp\
	./HTS.hpp\
	./HilbertTriplet.hpp\
	./Matrix2HilbertCoordinates.hpp\
	./HilbertTripletCompare.hpp\
	./BICRS.hpp\
	./CBICRS.hpp\
	./CCSWrapper.hpp\
	./Hilbert.hpp\
	./BlockHilbert.hpp\
	./HBICRS.hpp\
	./BetaHilbert.hpp\
	./FBICRS.hpp\
	./MachineInfo.hpp\
	./RDBHilbert.hpp\
	./BigInt.hpp\
	./BisectionHilbert.hpp\
	./RDScheme.hpp\
	./McCRS.hpp\
	./CompressedHilbert.hpp\
	./custom_quicksort.hpp\
	./HilbertArray.hpp\
	./MKLCRS.hpp\
	./CuHyb.hpp\
	./RDCSB.hpp

driver.o : driver.cpp $(DRIVER_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o driver.o driver.cpp



MICBENCH_DEPS = \
	./Matrix.hpp\
	./vecBICRS.hpp\
	./alignment.hpp\
	./Triplet.hpp\
	./ULIdef.hpp\
	./SparseMatrix.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./FBICRS.hpp\
	./ICRS.hpp\
	./CRS.hpp\
	./RDBHilbert.hpp\
	./Matrix2HilbertCoordinates.hpp\
	./MachineInfo.hpp\
	./BigInt.hpp\
	./TS.hpp\
	./RDCSB.hpp

micbench.o : micbench.cpp $(MICBENCH_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o micbench.o micbench.cpp



MM2HTRP_DEPS = \
	./HTS.hpp\
	./SparseMatrix.hpp\
	./ULIdef.hpp\
	./Matrix.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./Triplet.hpp\
	./alignment.hpp\
	./HilbertTriplet.hpp\
	./Matrix2HilbertCoordinates.hpp\
	./HilbertTripletCompare.hpp

mm2htrp.o : mm2htrp.cpp $(MM2HTRP_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o mm2htrp.o mm2htrp.cpp



MMIO_DEPS = \
	./mmio.h

mmio.o : mmio.cpp $(MMIO_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o mmio.o mmio.cpp



MTX2CRSTRP_DEPS = \
	./SparseMatrix.hpp\
	./ULIdef.hpp\
	./Matrix.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./Triplet.hpp\
	./alignment.hpp\
	./TS.hpp

mtx2crstrp.o : mtx2crstrp.cpp $(MTX2CRSTRP_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o mtx2crstrp.o mtx2crstrp.cpp



SBD2TRP-MAIN_DEPS = \
	./util.hpp\
	./MinCCS.hpp\
	./BlockOrderer.hpp\
	./SBDTree.hpp\
	./Triplet.hpp\
	./ULIdef.hpp\
	./BlockCRS.hpp\
	./U2.hpp\
	./Duck.hpp\
	./SepLast.hpp\
	./MinCRS.hpp\
	./FileToVT.hpp\
	./mmio.h

sbd2trp-main.o : sbd2trp-main.cpp $(SBD2TRP-MAIN_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o sbd2trp-main.o sbd2trp-main.cpp



SBDMV_DEPS = \
	./Triplet.hpp\
	./ULIdef.hpp\
	./FileToVT.hpp\
	./mmio.h\
	./HBICRS.hpp\
	./Matrix.hpp\
	./SparseMatrix.hpp\
	./alignment.hpp\
	./TS.hpp\
	./CRS.hpp\
	./ICRS.hpp\
	./ZZ_CRS.hpp\
	./ZZ_ICRS.hpp\
	./SVM.hpp\
	./HTS.hpp\
	./HilbertTriplet.hpp\
	./Matrix2HilbertCoordinates.hpp\
	./HilbertTripletCompare.hpp\
	./BICRS.hpp\
	./CCSWrapper.hpp\
	./MinCCS.hpp\
	./BlockOrderer.hpp\
	./SBDTree.hpp\
	./MinCRS.hpp\
	./BlockCRS.hpp\
	./U2.hpp\
	./Duck.hpp\
	./SepLast.hpp\
	./util.hpp

sbdmv.o : sbdmv.cpp $(SBDMV_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o sbdmv.o sbdmv.cpp



UTIL_DEPS = \
	./util.hpp

util.o : util.cpp $(UTIL_DEPS)
	@clang++ -c -std=c++11 -Wall -fopenmp   -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -DNDEBUG $(CXXFLAGS) -o util.o util.cpp


