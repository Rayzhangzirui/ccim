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



cim: cim.o tryvn.o helper.o icim.o amgsolver.o
	$(CC) $(OPT) $(HYPRE_LIBS) -o $@ $^ 

cim.o: cim.cpp
	$(CC) $(OPT) -c $<

tryvn.o: tryvn.C extratest.h
	$(CC) $(OPT) -c $<

test: test.o tests_main.o tryvn.o helper.o icim.o
	$(CC) $(OPT) -o $@ $^ 

test.o: test.cpp
	$(CC) $(OPT) -c $< 

helper.o: helper.cpp helper.h
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

amgsolver.o: amgsolver.cpp
	$(CC) $(OPT) $(CFLAGS) -c $<


.PHONY: clean
clean:
	rm  *.o