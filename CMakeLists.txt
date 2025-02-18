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

# Créer les dossiers pour les résultats et données - avec permissions write
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data)

# S'assurer que les dossiers ont les bonnes permissions
execute_process(COMMAND chmod -R 777 ${PROJECT_BINARY_DIR}/results)
execute_process(COMMAND chmod -R 777 ${PROJECT_BINARY_DIR}/data)

# Copier tous les fichiers data.txt du répertoire source vers le répertoire build
file(GLOB DATA_FILES ${PROJECT_SOURCE_DIR}/*.txt ${PROJECT_SOURCE_DIR}/data/*.txt)
foreach(DATA_FILE ${DATA_FILES})
    get_filename_component(FILENAME ${DATA_FILE} NAME)
    configure_file(${DATA_FILE} ${PROJECT_BINARY_DIR}/${FILENAME} COPYONLY)
    # Si le fichier est dans un sous-dossier data, le copier aussi dans data/
    if(${DATA_FILE} MATCHES ".*data/.*\\.txt$")
        get_filename_component(FILENAME ${DATA_FILE} NAME)
        configure_file(${DATA_FILE} ${PROJECT_BINARY_DIR}/data/${FILENAME} COPYONLY)
    endif()
endforeach()

# Configuration pour GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Target pour nettoyer les exécutables et les fichiers générés
add_custom_target(cleanall
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/results
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_NAME}
    COMMAND ${CMAKE_COMMAND} -E echo "Clean completed"
)

# Ajouter le sous-répertoire des tests
add_subdirectory(tests)