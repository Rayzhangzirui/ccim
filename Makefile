HOSTNAME = $(shell echo $$HOSTNAME)

CC=g++ -std=c++11 

HYPRE_DIR = /Users/zzirui/hypre-2.11.2/src/hypre
HYPRE_LIBS     = -L$(HYPRE_DIR)/lib -lHYPRE -lm

CINCLUDES = -I$(HYPRE_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS    = $(OPTS) $(CINCLUDES) $(CDEFS)

debug=0
$(info    debug is $(debug))
$(info    CC is $(CC))
ifeq ($(debug), 1)
	OPT = -g 
else
	OPT = -Ofast -lomp
endif

SOLVER = amg.o gmres.o solver.o hypresolver.o

METHOD = cim12.o iim.o icim.o cim345cond.o ccim.o

cim: cim.o global.o tryvn.o helper.o  matrix.o input.o sparse.o advance.o storage.o march.o numerics.o interface.o pb.o $(SOLVER) $(METHOD)
	$(CC) $(OPT) $(HYPRE_LIBS) -o $@ $^ 


epde: epde.o matrix.o input.o sparse.o storage.o numerics.o interface.o pb.o hypresolver.o ccim.o global.o helper.o
	$(CC) $(OPT) $(HYPRE_LIBS) -o $@ $^ 

epde.o: epde.cpp
	$(CC) $(OPT) -c $<


solvers: $(SOLVER)
	$(CC) $(OPT) -o $@ $^ 

methods: $(METHOD)
	$(CC) $(OPT) -o $@ $^ 

cim.o: cim.cpp
	$(CC) $(OPT) -c $<

test: test.o tests_main.o tryvn.o helper.o icim.o
	$(CC) $(OPT) -o $@ $^ 

test.o: test.cpp
	$(CC) $(OPT) -c $< 

tests_main.o : tests_main.cpp
	$(CC) $(OPT) -c $< 

# icim 
icim.o: icim.cpp icim.h
	$(CC) $(OPT) -c $<	

test_icim.o: test_icim.cpp
	$(CC) $(OPT) -c $<	

test_icim: test_icim.o tests_main.o tryvn.o helper.o icim.o
	$(CC) $(OPT) -o $@ $^ 

icim: icimmain.o tryvn.o helper.o icim.o
	$(CC) $(OPT) -o $@ $^ 

hypresolver.o: hypresolver.cpp
	$(CC) $(OPT) $(CFLAGS) -c $<


%.o: %.cpp %.h
	$(CC) $(OPT) -c $< 


.PHONY: clean
clean:
	rm  *.o