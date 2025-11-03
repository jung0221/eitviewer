cmake_minimum_required(VERSION 3.16)
# Enable C as well because some dependencies (VTK) may require MPI C targets
project(eitviewer LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# This project requires Qt6 exclusively. If VTK or any other imported target
# already brought Qt5 (for example VTK was built against Qt5), we fail early
# with a clear message. Mixing Qt5 and Qt6 in one target is unsupported.

find_package(MPI REQUIRED)

# VTK: require the Qt GUI support and rendering
find_package(VTK REQUIRED COMPONENTS 
  CommonCore
  FiltersCore
  FiltersGeometry
  FiltersGeneral
  FiltersExtraction
  IOXML
  IOGeometry
  GUISupportQt
  RenderingCore
  RenderingOpenGL2
  RenderingAnnotation
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "VTK not found. Install VTK with Qt support (GUISupportQt).")
endif()

# Require Qt6 and fail if any Qt5 target is present (commonly caused when
# VTK or other libs were built against Qt5).
find_package(Qt6 REQUIRED COMPONENTS Widgets)

if (TARGET Qt5::Core)
  message(FATAL_ERROR "A Qt5 target (Qt5::Core) was detected in the CMake graph. This project requires Qt6 exclusively.\nPossible causes: VTK or other packages were built against Qt5.\nOptions: rebuild VTK and other dependencies with Qt6, install matching VTK/Qt6 packages, or adjust CMAKE_PREFIX_PATH/Qt6_DIR to point to a Qt6 installation that matches your VTK.\n")
endif()

message(STATUS "Selected Qt major version: 6 (Qt6 found and required)")

add_executable(${PROJECT_NAME}
  src/main.cpp
  src/ConductivityViewer.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE src)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
    Qt6::Widgets
    ${VTK_LIBRARIES}
)

# For UIC/MOC if needed later
set_target_properties(${PROJECT_NAME} PROPERTIES
  AUTOMOC ON
  AUTOUIC ON
  AUTORCC ON
)

# Helpful packaging / compile definitions from VTK
if (VTK_USE_FILE)
  include(${VTK_USE_FILE})
endif()

message(STATUS "Configuração CMake concluída. Use: mkdir build && cd build && cmake .. && cmake --build .")
