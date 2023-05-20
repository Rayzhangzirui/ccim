HOSTNAME = $(shell echo $$HOSTNAME)

CC=g++ -std=c++11 
# CC=mpic++ -std=c++11 
# CC=mpicc -std=c++11 


debug=0
$(info    debug is $(debug))
$(info    CC is $(CC))

ifeq ($(debug), 1)
	OPT = -g 
else
	OPT = -Ofast
endif

usehypre=1
# ifeq ($(usehypre), 1)
# 	HYPRE_DIR =  /Users/Ray/project/hypre/src/hypre
# 	HYPRE_LIBS = -L$(HYPRE_DIR)/lib -lHYPRE -lm
# 	HYPRE_INC = -I$(HYPRE_DIR)/include -DHAVE_CONFIG_H -DHYPRE_TIMING
# 	OPT+=-DUSEHYPRE
# endif

SOLVER = amg.o gmres.o solver.o hypresolver.o

METHOD = cim12.o iim.o icim.o ccim.o

COMMON = global.o matrix.o input.o sparse.o finitediff.o numerics.o interface.o pb.o helper.o storage.o

MOTION = storage.o march.o getvn.o

# elliptic PDE solver
epde: epde.o $(COMMON) $(SOLVER) $(METHOD)
	$(CC) $(OPT) $(HYPRE_LIBS) -o $@ $^ 

# Motion under jump of normal derivative example
motion: motion.o $(COMMON)  $(SOLVER) $(MOTION) $(METHOD)
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
