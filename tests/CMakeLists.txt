if(DEFINED ENV{GTEST_ROOT})
  set(GTEST_ROOT $ENV{GTEST_ROOT})
endif()

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Créer l'exécutable de test
add_executable(tests 
    tests.cpp
    ${PROJECT_SOURCE_DIR}/src/Matrix.cpp
    ${PROJECT_SOURCE_DIR}/src/solver.cpp
    ${PROJECT_SOURCE_DIR}/src/functions.cpp
)

# Inclure les headers du projet
target_include_directories(tests PRIVATE 
    ${PROJECT_SOURCE_DIR}/src
)

# Lier les bibliothèques
target_link_libraries(tests
    ${MPI_LIBRARIES}
    GTest::GTest
    GTest::Main
    pthread
)