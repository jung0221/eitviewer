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
  ConductivityViewer(QWidget* parent = nullptr, const QString& initialLeftMesh = QString(), const QString& initialRightMesh = QString());
  ~ConductivityViewer() override = default;

private slots:
  void openMesh();
  void onSliceSliderChanged(int value);
  void onRenderModeChanged(int index);
  void onOpacityChanged(int value);
  void resetCamera();

private:
  void setupUi();
  bool loadMeshFile(const QString& path, int side = 0);

  // UI / widgets
  QVTKOpenGLNativeWidget* vtkWidget_ = nullptr;
  QVTKOpenGLNativeWidget* vtkWidgetRight_ = nullptr;
  QSlider* sliceSlider_ = nullptr;
  QComboBox* renderModeCombo_ = nullptr;
  QSlider* opacitySlider_ = nullptr;

  // VTK pipeline objects
  vtkSmartPointer<vtkRenderer> renderer_;
  vtkSmartPointer<vtkRenderer> rendererRight_;
  vtkSmartPointer<vtkActor> surfaceActor_;
  vtkSmartPointer<vtkActor> surfaceActorRight_;
  vtkSmartPointer<vtkPolyDataMapper> surfaceMapper_;
  vtkSmartPointer<vtkPolyDataMapper> surfaceMapperRight_;
  vtkSmartPointer<vtkLookupTable> lookupTable_;
  vtkSmartPointer<vtkLookupTable> lookupTableRight_;
  vtkSmartPointer<vtkScalarBarActor> scalarBar_;
  vtkSmartPointer<vtkScalarBarActor> scalarBarRight_;
  // slice pipeline
  vtkSmartPointer<vtkPlane> slicePlane_;
  vtkSmartPointer<vtkPlane> slicePlaneRight_;
  vtkSmartPointer<vtkCutter> sliceCutter_;
  vtkSmartPointer<vtkCutter> sliceCutterRight_;
  vtkSmartPointer<vtkActor> sliceActor_;
  vtkSmartPointer<vtkActor> sliceActorRight_;
  vtkSmartPointer<vtkPolyDataMapper> sliceMapper_;
  vtkSmartPointer<vtkPolyDataMapper> sliceMapperRight_;
  vtkSmartPointer<vtkProbeFilter> sliceProbe_;
  vtkSmartPointer<vtkProbeFilter> sliceProbeRight_;
  // slice axis (0=x,1=y,2=z)
  int sliceAxis_ = 2;
  int sliceAxisRight_ = 2;
  double dataBounds_[6] = {0,0,0,0,0,0};
  double dataBoundsRight_[6] = {0,0,0,0,0,0};

  QString currentMeshPath_;
  QString currentMeshPathRight_;

  // slots to open left/right meshes
private slots:
  void openPredictedMesh();
  void openTargetMesh();
};