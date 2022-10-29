set(TIRA_SOURCE main.cpp tira.cpp tools.cpp)

add_executable(tira ${TIRA_SOURCE})

target_include_directories(tira PUBLIC ${CMAKE_SOURCE_DIR}/src)

find_library(GINAC_LIBRARY ginac)
if(NOT GINAC_LIBRARY)
    message(FATAL_ERROR "GiNaC library not found.")
endif()

target_link_libraries(tira ${GINAC_LIBRARY})