CC=gcc-4.5
CPP=g++-4.5

OPT=-O3
DEBUG=-g -G -O0
CUDA_CC_FLAGS=-w -arch=sm_20  -Xptxas -dlcm=ca -Xcompiler -rdynamic -Xcompiler -fopenmp -maxrregcount 42#-Xptxas -v #
LD_LIBS =
LIB_INCLUDE=-I include -I src/include -I /usr/local/cuda/include -I ~/NVIDIA_GPU_Computing_SDK/C/common/inc -I metis-4.0.3/Lib

#default to debug builds for now
OPT_FLAGS=$(DEBUG) 

all: lib/levelset.a levelset/levelset_d

profile: 
	$(MAKE) DEFINES=-DPROFILE

verbose: 
	$(MAKE) DEFINES=-DDEBUG

opt:
	$(MAKE) OPT_FLAGS="$(OPT)"

debug:
	$(MAKE) OPT_FLAGS="$(DEBUG)"

SRCS=$(wildcard src/*.cu) $(wildcard src/*/*.cu) $(wildcard src/*/*/*.cu)
OBJS=$(SRCS:.cu=.o)
INCLUDE_FILES=$(wildcard include/*.h) $(wildcard include/*/*.h) $(wildcard include/*/*/*.h) $(wildcard src/include/*.h)

%.o : %.cu $(INCLUDE_FILES)
	    nvcc ${LIB_INCLUDE} ${OPT_FLAGS} ${CXX_FLAGS} ${CUDA_CC_FLAGS} -c $< -o $@ 

lib/levelset.a: lib $(OBJS)
	ar ruv lib/levelset.a $(OBJS)

lib:
	mkdir lib


.PHONY clean:
				rm -f $(OBJS) \
				rm -f lib/levelset.a \
				rm -f levelset/levelset_d

levelset/levelset_d: levelset/levelset.cu lib/levelset.a metis-4.0.3/libmetis.a
	nvcc -o levelset/levelset_d ${LIB_INCLUDE} ${OPT_FLAGS} ${LD_LIBS} ${CUDA_CC_FLAGS} levelset/levelset.cu lib/levelset.a --linker-options metis-4.0.3/libmetis.a 

#examples/cusp_sa: examples/CUSP_smoothed_aggregation.cu lib/amg.a
#	nvcc -o examples/cusp_sa -I include ${OPT_FLAGS} ${CUDA_CC_FLAGS} lib/amg.a examples/CUSP_smoothed_aggregation.cu 

#examples/BoomerAMG: examples/BoomerAMG.cu lib/amg.a
#	nvcc -o examples/BoomerAMG -I include -I /home/sci/zhisong/hypre-2.8.0b/src/hypre/include -I /usr/lib64/mpi/gcc/openmpi/include ${OPT_FLAGS} ${CUDA_CC_FLAGS} lib/amg.a examples/BoomerAMG.cu \
#				--linker-options /home/sci/zhisong/hypre-2.8.0b/src/hypre/lib64/libHYPRE.a \
#				--linker-options /usr/lib64/mpi/gcc/openmpi/lib64/libmpi.so
