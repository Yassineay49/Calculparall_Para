cmake_minimum_required(VERSION 3.0)
project(CHP)

# Configuration C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Chemins des bibliothèques
set(LIBRARY_PATH_EIGEN "../../../libraries/eigen/")

# Options de compilation
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_BUILD_TYPE "Release")

# Trouver MPI
find_package(MPI REQUIRED)

# Définir les chemins des sources
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)

# Liste des fichiers sources
set(PROJECT_SOURCES
    ${SRC_DIR}/settings.cpp
    ${SRC_DIR}/charge.cpp
    ${SRC_DIR}/functions.cpp
    ${SRC_DIR}/Matrix.cpp
    ${SRC_DIR}/solver.cpp
    ${PROJECT_SOURCE_DIR}/main.cpp
)

# Créer la bibliothèque pour les tests
add_library(solver_lib
    ${SRC_DIR}/settings.cpp
    ${SRC_DIR}/charge.cpp
    ${SRC_DIR}/functions.cpp
    ${SRC_DIR}/Matrix.cpp
    ${SRC_DIR}/solver.cpp
)

# Configurer la bibliothèque
target_include_directories(solver_lib PUBLIC 
    ${SRC_DIR}
    ${LIBRARY_PATH_EIGEN}
    "src/libraries"
)
target_link_libraries(solver_lib PUBLIC MPI::MPI_CXX)

# Créer l'exécutable principal
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})

# Configurer l'exécutable
target_include_directories(${PROJECT_NAME} PUBLIC 
    ${SRC_DIR}
    ${LIBRARY_PATH_EIGEN}
    "src/libraries"
)
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)

# Créer les dossiers pour les résultats et données
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data)

# Ajouter le sous-répertoire des tests
add_subdirectory(tests)