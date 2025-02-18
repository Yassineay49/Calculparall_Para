# Fichier CMakeLists.txt pour le dossier tests/
cmake_policy(SET CMP0079 NEW)

# Configuration pour les tests
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Création de l'exécutable de test avec un nom différent
add_executable(run_tests tests.cpp)
target_link_libraries(run_tests 
    solver_lib
    ${GTEST_LIBRARIES}
    ${GTEST_MAIN_LIBRARIES}
    pthread
    MPI::MPI_CXX
)

# S'assurer que les dossiers de résultats existent
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)

# Donner les permissions d'écriture au dossier results
execute_process(COMMAND chmod -R 777 ${CMAKE_CURRENT_BINARY_DIR}/results)

# Créer une cible de nettoyage pour les tests uniquement
add_custom_target(clean_tests
    COMMAND ${CMAKE_COMMAND} -E remove run_tests
    COMMAND ${CMAKE_COMMAND} -E remove_directory results
    COMMAND ${CMAKE_COMMAND} -E echo "Tests cleaned"
)