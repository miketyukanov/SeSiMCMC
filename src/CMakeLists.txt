cmake_minimum_required(VERSION 3.23)

project(SeSiMCMC VERSION 1.0)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

add_executable(SeSiMCMC)

set(common_SeSiSources confread.c Random.c MCMC.cpp ResultsSet.cpp
    Sequences.cpp SymbolsCounter.cpp)
target_sources(SeSiMCMC PRIVATE ${common_SeSiSources} gibbslfm.cpp)
target_sources(SeSiMCMC PRIVATE FILE_SET HEADERS FILES confread.h Random.h Atgc.hpp
    Diagnostics.hpp Exception.hpp KullbakCounter.hpp Logger.hpp MarkovChainState.hpp
    MCMC.hpp RandomMappings.hpp ResultsSet.hpp Sequences.hpp SymbolsCounter.hpp)



#Favorov's tests (those that compile)

add_executable(AtgcTest)
target_sources(AtgcTest PRIVATE ${common_SeSiSources} routine-tests/AtgcTest.cpp)

add_executable(DiagnosticsTest)
target_sources(DiagnosticsTest PRIVATE ${common_SeSiSources} routine-tests/DiagnosticsTest.cpp)

add_executable(MarkovChainStateTest)
target_sources(MarkovChainStateTest PRIVATE ${common_SeSiSources} routine-tests/MarkovChainStateTest.cpp)


#to do: make this test compileable
#add_executable(ResultSetTest)
#target_sources(ResultSetTest PRIVATE ${common_SeSiSources} routine-tests/ResultSetTest.cpp)

add_executable(SequencesPileTest)
target_sources(SequencesPileTest PRIVATE ${common_SeSiSources} routine-tests/SequencesPileTest.cpp)

add_executable(StupidTest)
target_sources(StupidTest PRIVATE ${common_SeSiSources} routine-tests/StupidTest.cpp)

add_executable(SymbolsCounterTest)
target_sources(SymbolsCounterTest PRIVATE ${common_SeSiSources} routine-tests/SymbolsCounterTest.cpp)
