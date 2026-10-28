# CMake 3.28 introduced supported dependency scanning for C++20 modules.
# Keep this requirement aligned with the CXX_MODULES file set below.
cmake_minimum_required(VERSION 4.4)
project(ray1 VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#set(CMAKE_CXX_FLAGS "/W4")
#set(CMAKE_CXX_FLAGS_DEBUG "/g")
# These are MSVC switches.  Passing them to Emscripten's Clang makes its
# dependency scanner interpret them as input files.
if (MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /wd5244")
endif()

if (NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    set(CMAKE_PREFIX_PATH "C:/Qt/6.12.0/msvc2022_64" CACHE PATH "Qt prefix for native builds")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON) 
set(CMAKE_CXX_EXTENSIONS ON)

# `FlipFluid` is imported by normal C++ sources (indirectly through fluid.h).
# Scan every C++ source so CMake supplies MSVC with the generated BMI/IFC.
set(CMAKE_CXX_SCAN_FOR_MODULES ON)

#set(CMAKE_EXPERIMENTAL_CXX_MODULE_DYNDEP 1)
# turn on the experimental API
#set(CMAKE_EXPERIMENTAL_CXX_MODULE_CMAKE_API  "2182bf5c-ef0d-489a-91da-49dbc3090d2a")


# Set up AUTOMOC and some sensible defaults for runtime execution
# When using Qt 6.3, you can replace the code block below with
# qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
include(GNUInstallDirs)
set(CMAKE_AUTOUIC ON)

find_package(QT NAMES Qt6 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Concurrent Gui Widgets OpenGLWidgets)


qt_add_executable(ray1 WIN32
    Ray.h
    calc.cpp calc.h
    camera.h
    hitPosition.cpp hitPosition.h
    main.cpp
    FluidBenchmark.cpp FluidBenchmark.h
    rayview.cpp rayview.h rayview.ui
    fluid.cpp fluid.h fluid.ui
    shapes.h 
    FlipFluid.cxx
)

# --- Configuration pour l'analyse de performance sous MSVC (Visual Studio) ---
if(MSVC)
    # 1. Force la génération des symboles (.pdb) en configuration Release
    set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Release>:ProgramDatabase>")
    
    # 2. Configure l'Éditeur de liens (Linker) pour générer les symboles et optimiser le binaire
    target_link_options(ray1 PRIVATE 
        "$<$<CONFIG:Release>:/DEBUG>"          # Génère le fichier PDB
        "$<$<CONFIG:Release>:/OPT:REF>"        # Élimine les fonctions jamais appelées
        "$<$<CONFIG:Release>:/OPT:ICF>"        # Fusionne les fonctions identiques (COMDAT)
        "$<$<CONFIG:Release>:/INCREMENTAL:NO>" # Désactive le lien incrémental (qui fausse les mesures)
    )
endif()

target_sources(ray1 
    PRIVATE
        FILE_SET cxx_modules TYPE CXX_MODULES FILES
        FlipFluid.ixx 
        mdspan.ixx
        fluid/BaseStructures.ixx
        fluid/CellCalculations.ixx 
        fluid/CellColor.ixx
        fluid/Constants.ixx
        fluid/CellVelocities.ixx
        fluid/Obstacle.ixx
        fluid/ParticleCollision.ixx
        fluid/ParticleColor.ixx
        fluid/ParticleDensity.ixx 
        fluid/ParticleIntegration.ixx 
        fluid/ParticleIncompressibility.ixx
        fluid/ParticleVelocities.ixx
        fluid/PushParticles.ixx
)

target_link_libraries(ray1 PRIVATE
    Qt::Concurrent
    Qt::Core
    Qt::Gui
    Qt::Widgets
    Qt::OpenGLWidgets
)

install(TARGETS ray1
    BUNDLE DESTINATION .
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Consider using qt_generate_deploy_app_script() for app deployment if
# the project can use Qt 6.3. In that case rerun qmake2cmake with
# --min-qt-version=6.3.
