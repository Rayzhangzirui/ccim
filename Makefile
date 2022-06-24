HOSTNAME = $(shell echo $$HOSTNAME)

CC=g++ -std=c++11 
# CC=mpic++ -std=c++11 

HYPRE_DIR = /Users/zzirui/hypre-2.11.2/src/hypre
# HYPRE_DIR = /Users/Ray/project/hypre/src/hypre
HYPRE_LIBS     = -L$(HYPRE_DIR)/lib -lHYPRE -lm

CINCLUDES = -I$(HYPRE_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS    = $(OPTS) $(CINCLUDES) $(CDEFS) 

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

cim.o: cim.cpp
	$(CC) $(OPT) -c $<

motion.o: motion.cpp
	$(CC) $(OPT) -c $<
epde.o: epde.cpp
	$(CC) $(OPT) -c $<

motion.o: motion.cpp
	$(CC) $(OPT) -c $<

hypresolver.o: hypresolver.cpp hypresolver.h
	$(CC) $(OPT) $(CFLAGS) -c $<

%.o: %.cpp %.h
	$(CC) $(OPT) -c $< 


.PHONY: clean
clean:
	rm  *.o