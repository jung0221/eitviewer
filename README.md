eitviewer
===========

Simple Qt6 + VTK viewer for visualizing conductivity fields on 3D meshes.

Features (minimal demo):
- Qt6 window with an embedded VTK widget (QVTKOpenGLNativeWidget)
- Open legacy VTK / .vtu/.vtk meshes
- Display point scalar array as colors
- Slider to move a slicing plane and visualize a 2D slice via vtkCutter

Dependencies
- CMake >= 3.16
- Qt6 (Widgets)
- VTK built with Qt support (GUISupportQt) and RenderingOpenGL2

Quick build

```bash
mkdir build
cd build
cmake ..
cmake --build .
# run
./condviz_gui
```

Notes
- This is a minimal starter. On your machine you may need to install VTK and Qt6 from your distro or build VTK with Qt support. On Ubuntu, packages or vcpkg/conan can help.
- The UI is intentionally small: open a mesh file from the menu, select a scalar array (if multiple), move the slider to change the slice position.

Next steps (suggested)
- Add support for Gmsh (.msh) files via vtkGmshReader if available
- Add colormap presets and automatic scalar range estimation
- Add isosurface extraction and opacity controls
- Add multiple views (3D + orthographic slices)
