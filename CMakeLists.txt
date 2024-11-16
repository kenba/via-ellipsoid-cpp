# Copyright (c) 2019-2024 Ken Barker
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required (VERSION 3.30)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0167 OLD) # For MSVC
project (via-ellipsoid LANGUAGES CXX)

option(INSTALL_PYTHON "Install Python Library." ON)
option(CPP_UNIT_TESTS "Build C++ Unit Tests." OFF)
option(CODE_COVERAGE "Add gcc code coverage options." OFF)

add_library(${PROJECT_NAME} INTERFACE)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER ellipsoid.hpp)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_23)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

if (WIN32)
  # Add the C++ Guidelines Support Library include files
  set(CPP_GSL_DIR $ENV{CPP_GSL_DIR})
  include_directories(${CPP_GSL_DIR}/include)

  # Add the Eigen Library include files
  set(Eigen_DIR $ENV{Eigen_DIR})
  include_directories(${Eigen_DIR})

  # Add the Via Angle Library include files
  set(ViaAngle_DIR $ENV{ViaAngle_DIR})
  include_directories(${ViaAngle_DIR}/include)

  # Add the Via Sphere Library include files
  set(ViaSphere_DIR $ENV{ViaSphere_DIR})
  include_directories(${ViaSphere_DIR}/include)
endif()

if (INSTALL_PYTHON)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 CONFIG REQUIRED)

  pybind11_add_module(via_ellipsoid src/via_ellipsoid_python_bindings.cpp)
  target_compile_features(via_ellipsoid PRIVATE cxx_std_23)
  target_include_directories(via_ellipsoid PRIVATE
    $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>)
  install(TARGETS via_ellipsoid DESTINATION .)
endif(INSTALL_PYTHON)

if (CPP_UNIT_TESTS)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)
    if(Boost_FOUND)
    add_executable(${PROJECT_NAME}_test "")

    file(GLOB_RECURSE INCLUDE_FILES include/via/*.hpp)
    target_sources(${PROJECT_NAME}_test
        PRIVATE
        ${INCLUDE_FILES}
    )

    target_sources(${PROJECT_NAME}_test
        PRIVATE
        tests/test_main.cpp

        tests/test_coefficients.cpp
        tests/test_ellipsoid_functions_double.cpp
        tests/test_Ellipsoid.cpp
        tests/test_geodesic_functions_double.cpp
        tests/test_Geodesic_double.cpp
        tests/test_ellipsoid_double.cpp
    )

    target_compile_definitions(${PROJECT_NAME}_test PRIVATE BOOST_TEST_DYN_LINK)
    target_include_directories(${PROJECT_NAME}_test PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME}_test
                            PRIVATE
                            ${PROJECT_NAME}
                            Boost::unit_test_framework)

    if (MSVC)
        target_compile_options(${PROJECT_NAME}_test PRIVATE /W4)
    else()
        target_compile_options(${PROJECT_NAME}_test PRIVATE -Wall -Wextra -Wpedantic)

        if (CODE_COVERAGE)
        target_compile_options(${PROJECT_NAME}_test PRIVATE --coverage)
        target_link_libraries(${PROJECT_NAME}_test PRIVATE --coverage)
        endif()
    endif()

    enable_testing()
    add_test(NAME via-ellipsoid_test COMMAND  ${PROJECT_NAME}_test)
    endif(Boost_FOUND)
endif(CPP_UNIT_TESTS)

# Install headers:
include(GNUInstallDirs)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/via" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
