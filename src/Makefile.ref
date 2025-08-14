#****************************************************************************#
#SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001
#$Id$
#****************************************************************************#
#this file can be started only from parent directory


.PHONY:	all tests objs clean tgz tests

all: ccvars gibbslfm.mak
	@if [ ! -d ./obj ]; then mkdir ./obj && echo "Creating ./obj...";fi
	$(MAKE) -f gibbslfm.mak

clean:
	$(MAKE) -f gibbslfm.mak clean

tgz:
	$(MAKE) -f gibbslfm.mak tgz

tests:
	$(MAKE) -f routine-tests/AtgcTest.mak
	$(MAKE) -f routine-tests/SequencesPileTest.mak
	$(MAKE) -f routine-tests/SymbolsCounterTest.mak
	$(MAKE) -f routine-tests/MarkovChainStateTest.mak
	$(MAKE) -f routine-tests/MotifsSetTest.mak
