cmake_minimum_required(VERSION 3.10)
project(RunTests)

# Configurer quelques politiques pour éviter les avertissements
cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0028 NEW)
cmake_policy(SET CMP0003 NEW)

# Trouver les paquets requis
find_package(MPI REQUIRED)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

# Inclure les répertoires d'en-têtes nécessaires
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${MPI_CXX_INCLUDE_PATH})

# Création de l'exécutable de test
add_executable(run_tests tests.cpp)

# Lier les bibliothèques nécessaires en utilisant les cibles modernes
target_link_libraries(run_tests 
    solver_lib           # Bibliothèque personnalisée, supposée accessible dans le chemin de linkage
    MPI::MPI_CXX         # Cible MPI fournie par find_package(MPI)
    GTest::GTest         # Cible GTest
    GTest::Main          # Cible GTest Main
    Threads::Threads     # Cible Threads
)

# S'assurer que le dossier results existe
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)

# Donner les permissions d'écriture au dossier results
execute_process(COMMAND chmod -R 777 ${CMAKE_CURRENT_BINARY_DIR}/results)

# Créer une cible de nettoyage pour les tests uniquement
add_custom_target(clean_tests
    COMMAND ${CMAKE_COMMAND} -E remove run_tests
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/results
    COMMAND ${CMAKE_COMMAND} -E echo "Tests cleaned"
)

