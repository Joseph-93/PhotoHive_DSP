cmake_minimum_required(VERSION 3.20)
project(PhotoHive-DSP C)

option(BUILD_TESTS "Build the test suite" OFF)
option(DEBUGGING "Build for debugging" ON)

set(CMAKE_C_STANDARD 11)
if(DEBUGGING)
    set(CMAKE_BUILD_TYPE Debug)
else()
    set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast")
endif()

set(FFTW3_DIR "/usr/local/lib/cmake/fftw3")
find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})

set(SOURCE_FILES
    image_processing.c
    fft_processing.c
    color_quantization.c
    blur_profile.c
    utilities.c
    filtering.c
    debug.c
)

find_library(FFTW3_THREADS_LIBRARY NAMES fftw3_threads PATHS /usr/local/lib)

if(BUILD_TESTS)
    set(TEST_FILES
        test/test.c
        interface.c
        ${SOURCE_FILES}
    )
    add_executable(test_suite ${TEST_FILES})
    target_link_libraries(test_suite ${FFTW3_LIBRARIES} ${FFTW3_THREADS_LIBRARY} m)
    target_compile_definitions(test_suite PRIVATE TESTING=1)
else()
    add_library(report_data SHARED ${SOURCE_FILES} interface.c)
    target_link_libraries(report_data ${FFTW3_LIBRARIES} ${FFTW3_THREADS_LIBRARY} m)
    set_target_properties(report_data PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../PhotoHive_DSP_lib"
    )
endif()
