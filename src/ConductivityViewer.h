#pragma once

#include <QMainWindow>
#include <QSlider>
#include <QComboBox>
#include <QString>

#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>

class vtkRenderer;
class vtkActor;
class vtkDataSet;
class vtkDataSetSurfaceFilter;
class vtkPolyDataMapper;
class vtkLookupTable;
class vtkScalarBarActor;
class vtkPlane;
class vtkCutter;
class vtkProbeFilter;

class ConductivityViewer : public QMainWindow {
  Q_OBJECT
public:
  ConductivityViewer(QWidget* parent = nullptr, const QString& initialMesh = QString());
  ~ConductivityViewer() override = default;

private slots:
  void openMesh();
  void onSliceSliderChanged(int value);
  void onRenderModeChanged(int index);
  void onOpacityChanged(int value);
  void resetCamera();

private:
  void setupUi();
  bool loadMeshFile(const QString& path);

  // UI / widgets
  QVTKOpenGLNativeWidget* vtkWidget_ = nullptr;
  QSlider* sliceSlider_ = nullptr;
  QComboBox* renderModeCombo_ = nullptr;
  QSlider* opacitySlider_ = nullptr;

  // VTK pipeline objects
  vtkSmartPointer<vtkRenderer> renderer_;
  vtkSmartPointer<vtkActor> surfaceActor_;
  vtkSmartPointer<vtkPolyDataMapper> surfaceMapper_;
  vtkSmartPointer<vtkLookupTable> lookupTable_;
  vtkSmartPointer<vtkScalarBarActor> scalarBar_;
  // slice pipeline
  vtkSmartPointer<vtkPlane> slicePlane_;
  vtkSmartPointer<vtkCutter> sliceCutter_;
  vtkSmartPointer<vtkActor> sliceActor_;
  vtkSmartPointer<vtkPolyDataMapper> sliceMapper_;
  vtkSmartPointer<vtkProbeFilter> sliceProbe_;
  // slice axis (0=x,1=y,2=z)
  int sliceAxis_ = 2;
  double dataBounds_[6] = {0,0,0,0,0,0};

  QString currentMeshPath_;
};