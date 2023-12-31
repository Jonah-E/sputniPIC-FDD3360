# From https://x.momo86.net/?p=29

CXX=g++
CXXFLAGS=-std=c++11 -I./include -O3 -g -Xcompiler -Wall

NVCC=nvcc
ARCH=sm_70
NVCCFLAGS= -I./include -arch=$(ARCH) -std=c++11 -O3 -g -Xcompiler -Wall --compiler-bindir=$(CXX)

SRCDIR=src
SRCS=$(shell find $(SRCDIR) -name '*.cu' -o -name '*.cpp')

OBJDIR=src
OBJS=$(subst $(SRCDIR),$(OBJDIR), $(SRCS))
OBJS:=$(subst .cpp,_cpp.o,$(OBJS))
OBJS:=$(subst .cu,_cu.o,$(OBJS))

DATA_DIR := ./data
$(shell mkdir -p $(DATA_DIR))

BIN := ./bin
TARGET=sputniPIC.out
TARGET_GPU=sputniPIC_gpu.out
TARGET_COMP=sputniPIC_comp.out

ifdef USE_DOUBLE
	CXXFLAGS += -DUSE_DOUBLE
	NVCCFLAGS += -DUSE_DOUBLE
endif

cpu: dir $(BIN)/$(TARGET)

gpu: CXXFLAGS += -DUSE_GPU
gpu: dir $(BIN)/$(TARGET_GPU)

compare: CXXFLAGS += -DCOMPARE
compare: dir $(BIN)/$(TARGET_COMP)

dir: ${BIN}

${BIN}:
	mkdir -p $(BIN)

$(BIN)/$(TARGET_COMP): $(OBJS)
	$(NVCC) $(NVCCFLAGS) $+ -o $@

$(BIN)/$(TARGET_GPU): $(OBJS)
	$(NVCC) $(NVCCFLAGS) $+ -o $@

$(BIN)/$(TARGET): $(OBJS)
	$(NVCC) $(NVCCFLAGS) $+ -o $@

$(SRCDIR)/%_cu.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

$(OBJDIR)/%_cpp.o: $(SRCDIR)/%.cpp
	[ -d $(OBJDIR) ] || mkdir $(OBJDIR)
	$(NVCC) $(CXXFLAGS) $< -c -o $@

.PHONY: clean clean-data

clean:
	rm -rf $(OBJS)
	rm -rf $(TARGET)

clean-data:
	rm -rf $(DATA_DIR)/*
