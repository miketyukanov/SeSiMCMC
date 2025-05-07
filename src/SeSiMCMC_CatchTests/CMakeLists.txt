cmake_minimum_required(VERSION 3.16)

project(SMTests LANGUAGES CXX)

enable_testing()

SET(CMAKE_CXX_STANDARD 11)

include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.8.1
FIND_PACKAGE_ARGS
)

FetchContent_MakeAvailable(Catch2)

if ($ENV{CATCH2_INSTALL_DIR})
    get_filename_component(_CATCH2_BASE $ENV{CATCH2_INSTALL_DIR} REALPATH)
elseif (EXISTS "") # set by QC
    get_filename_component(_CATCH2_BASE "" REALPATH) # set by QC
endif()

if (NOT Catch2_DIR)
    if (_CATCH2_BASE)
        if (EXISTS "${_CATCH2_BASE}/lib/cmake/Catch2")
            set(_CATCH2_BASE "${_CATCH2_BASE}/lib/cmake/Catch2")
        endif()
        message("Setting Catch2_DIR to ${_CATCH2_BASE}")
        set(Catch2_DIR ${_CATCH2_BASE})
    else()
        message("No Catch2_DIR specified - using system defaults.")
    endif()
endif()

find_package(Catch2 3 REQUIRED)

add_executable(SMTests tst_smtests.cpp)
target_link_libraries(SMTests PRIVATE Catch2::Catch2WithMain PRIVATE SeSiMCMCLib)
target_include_directories(SMTests PRIVATE "$<TARGET_PROPERTY:SeSiMCMCLib,INCLUDE_DIRECTORIES>")

add_test(NAME SMTests COMMAND SMTests)
