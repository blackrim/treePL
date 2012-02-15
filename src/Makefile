
#CC := /opt/intel/composer_xe_2011_sp1/bin/icc
#CP := /opt/intel/composer_xe_2011_sp1/bin/icc
#probably need something like this if using intel
#export LD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.6.233/compiler/lib/intel64/:$LD_LIBRARY_PATH
#get some errors with intel with really big trees
CC := gcc
CP := g++

OPENMP := -fopenmp
#OPENMP := -openmp -parallel

LIBS := -lm $(OPENMP) -ladolc -L/usr/lib64 -lnlopt_cxx

RM := rm -rf

SRC_DIR = ./

prefix = /usr/bin/

CPP_SRCS += \
./main.cpp \
./node.cpp \
./pl_calc_parallel.cpp \
./tree.cpp \
./tree_reader.cpp \
./tree_utils.cpp \
./optimize_tnc.cpp \
./optimize_nlopt.cpp \
./myradops.cpp \
./siman_calc_par.cpp \
./optim_options.cpp \
./utils.cpp

C_SRCS += \
./tnc.c 

OBJS += \
./main.o \
./node.o \
./pl_calc_parallel.o \
./tnc.o \
./tree.o \
./tree_reader.o \
./tree_utils.o \
./optimize_tnc.o \
./optimize_nlopt.o \
./myradops.o \
./siman_calc_par.o \
./optim_options.o \
./utils.o

C_DEPS += \
./tnc.d 

CPP_DEPS += \
./main.d \
./node.d \
./pl_calc_parallel.d \
./tree.d \
./tree_reader.d \
./tree_utils.d \
./optimize_tnc.d \
./optimize_nlopt.d \
./myradops.d \
./siman_calc_par.d \
./optim_options.d \
./utils.d

%.o: %.cpp
	$(CP) -O3 -g3 -c -fmessage-length=0 $(OPENMP) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o $@ $(SRC_DIR)$<

%.o: %.c
	$(CC) -O3 -g3 -c -fmessage-length=0 $(OPENMP) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o $@ $<

all: treePL

treePL: $(OBJS)
	$(CP) -O3 -g3 -o treePL $(OBJS) $(LIBS)

clean:
	-$(RM) *.d *.o

install:
	install -m 0755 treePL $(prefix)
	
uninstall:
	-rm $(prefix)/treePL

