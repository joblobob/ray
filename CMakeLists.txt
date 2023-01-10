cmake_minimum_required(VERSION 3.20)
project(ray1 VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#set(CMAKE_CXX_FLAGS "/W4")
#set(CMAKE_CXX_FLAGS_DEBUG "/g")
set(CMAKE_CXX_FLAGS_RELEASE "/O2")

set(CMAKE_PREFIX_PATH C:/Qt/5.15.2/msvc2019_64)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set up AUTOMOC and some sensible defaults for runtime execution
# When using Qt 6.3, you can replace the code block below with
# qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
include(GNUInstallDirs)
set(CMAKE_AUTOUIC ON)

find_package(QT NAMES Qt5 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Concurrent Gui Widgets)

add_executable(ray1 WIN32
    Ray.h
    calc.cpp calc.h
    camera.h
    hitPosition.cpp hitPosition.h
    main.cpp
    rayview.cpp rayview.h rayview.ui
    fluid.h fluid.cpp fluid.ui
    shapes.h 
    Module1.ixx)

target_link_libraries(ray1 PRIVATE
    Qt::Concurrent
    Qt::Core
    Qt::Gui
    Qt::Widgets
)

install(TARGETS ray1
    BUNDLE DESTINATION .
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Consider using qt_generate_deploy_app_script() for app deployment if
# the project can use Qt 6.3. In that case rerun qmake2cmake with
# --min-qt-version=6.3.
