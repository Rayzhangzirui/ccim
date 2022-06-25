HOSTNAME = $(shell echo $$HOSTNAME)

CC=g++ -std=c++11 
# CC=mpic++ -std=c++11 

HYPRE_DIR = /Users/zzirui/hypre-2.11.2/src/hypre
# HYPRE_DIR = /Users/Ray/project/hypre/src/hypre
HYPRE_LIBS     = -L$(HYPRE_DIR)/lib -lHYPRE -lm

HYPRE_INC = -I$(HYPRE_DIR)/include -DHAVE_CONFIG_H -DHYPRE_TIMING

debug=1
$(info    debug is $(debug))
$(info    CC is $(CC))

ifeq ($(debug), 1)
	OPT = -g 
else
	OPT = -Ofast
endif


SOLVER = amg.o gmres.o solver.o hypresolver.o

METHOD = cim12.o iim.o icim.o cim345cond.o ccim.o

COMMON = global.o matrix.o input.o sparse.o finitediff.o numerics.o interface.o pb.o helper.o

MOTION = storage.o march.o getvn.o

cim: cim.o $(COMMON) tryvn.o  advance.o $(SOLVER) $(METHOD) $(MOTION)
	$(CC) $(OPT) $(HYPRE_LIBS) -o $@ $^ 

epde: epde.o $(COMMON) storage.o hypresolver.o ccim.o solver.o icim.o
	$(CC) $(OPT) $(HYPRE_LIBS) -o $@ $^ 

motion: motion.o $(COMMON) hypresolver.o solver.o $(MOTION) $(METHOD)
	$(CC) $(OPT) $(HYPRE_LIBS) -o $@ $^ 

hypresolver.o: hypresolver.cpp hypresolver.h
	$(CC) $(OPT) $(HYPRE_INC) -c $<

# %.o: %.cpp %.h
# 	$(CC) $(OPT) -c $< 


# makefile header dependency
# answer by Sophie https://stackoverflow.com/questions/2394609/makefile-header-dependencies?noredirect=1&lq=1
# https://codereview.stackexchange.com/questions/2547/makefile-dependency-generation/11109#11109
srcs = $(wildcard *.cpp)
objs = $(srcs:.cpp=.o)
deps = $(srcs:.cpp=.d)

%.o: %.cpp
	$(CC) $(OPT) -MMD -MP -c $< -o $@


.PHONY: clean
clean:
	rm  *.o

-include $(deps)
