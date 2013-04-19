
#flags to enable cuda
ifndef CUDA_PATH
	CUDA_PATH = /usr/local/cuda
endif
CUDA_FLG = -L$(CUDA_PATH)/lib -lcuda -lcudart -lm




#if 64bit os
ifeq ($(shell file `which pwd` | grep -c -o 64), 1)
	NVCC = nvcc -m64
else
	NVCC = nvcc
endif

.SUFFIXES: .cu .cc
	
.cc.cu:
	cat $< > $@


ecmp : ecmp.cu
	$(NVCC) -o $@ $< $(CUDA_FLG)

emu : ecmp.cc
	$(CXX) -o addvec addvec.cc $(CUDA_FLG)


clean :
	rm -rf *.o ecmp


