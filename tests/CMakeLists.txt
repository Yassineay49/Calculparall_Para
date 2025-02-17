# Version minimale de CMake requise
cmake_minimum_required(VERSION 3.10)

# Définir le projet
project(TestSolver)

# Trouver MPI, GTest et Threads
find_package(MPI REQUIRED)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

# Comme libsolver_lib.a se trouve dans le dossier parent, on définit une variable pour le chemin
set(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Ajoute le répertoire parent aux répertoires de recherche de l'éditeur de liens
link_directories(${PARENT_DIR})

# Créer l'exécutable des tests
add_executable(tests tests.cpp)

# Lier les bibliothèques nécessaires
target_link_libraries(tests
    solver_lib       # Cible trouvée grâce à link_directories
    MPI::MPI_CXX     # Assurez-vous que MPI est bien installé
    GTest::GTest
    GTest::Main
    Threads::Threads
)

# Activer les tests
enable_testing()
add_test(NAME solver_tests COMMAND tests)

# Copier les fichiers de données nécessaires pour les tests (depuis le dossier parent)
file(COPY ${PARENT_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Créer le dossier results
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)

