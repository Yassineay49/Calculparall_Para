find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Créer l'exécutable de test
add_executable(tests 
    tests.cpp
)

# Lier les bibliothèques
target_link_libraries(tests
    solver_lib
    ${MPI_LIBRARIES}
    ${GTEST_LIBRARIES}
    pthread
)

# Copier les fichiers de données nécessaires pour les tests
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Créer le dossier results
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)