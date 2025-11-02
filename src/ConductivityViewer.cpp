#include "ConductivityViewer.h"

#include <QAction>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QMenuBar>
#include <QSlider>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QApplication>
#include <QWidget>

#include <vtkActor.h>
#include <vtkCutter.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkLookupTable.h>
#include <vtkNew.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkClipDataSet.h>
#include <vtkPlaneSource.h>
#include <vtkProbeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>

#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellData.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkPointLocator.h>

#include <unordered_map>
#include <array>

#include <QVTKOpenGLNativeWidget.h>
#include <QProcess>
#include <QDir>
#include <QFileInfo>
#include <QTemporaryFile>
#include <QDebug>
#include <QDateTime>
#include <QLabel>

/**
 * @brief Load a mesh file and build the visualization pipelines.
 *
 * Attempts to convert .msh files via meshio or gmsh CLI when necessary,
 * reads the resulting dataset into VTK, builds surface and probe pipelines
 * and sets up mappers and lookup table for scalar coloring.
 *
 * @param path Path to mesh file (.msh, .vtu, .vtk, ...)
 * @return true on success, false on failure
 */
bool ConductivityViewer::loadMeshFile(const QString &path)
{
  // try several readers based on extension
  std::string ext = QFileInfo(path).suffix().toLower().toStdString();

  vtkSmartPointer<vtkDataSet> data;

  // Handle Gmsh .msh: first try converting with gmsh CLI to a .vtk file (no Python required).
  QString localPath = path;
  if (ext == "msh")
  {
    QString tmpdir = QDir::tempPath();
    QString base = QFileInfo(path).baseName();
    // creates a temporary mesh with vtu format (converted from .msh) 
    QString out = tmpdir + "/condviz_" + base + "_" + QString::number(QDateTime::currentSecsSinceEpoch()) + ".vtu";

    // Prefer using Python + meshio so we avoid launching gmsh GUI. Try several python executables
    QStringList pythonCandidates;
    // explicit env override
    QString pyEnv = qEnvironmentVariable("EIT_PYTHON");
    if (!pyEnv.isEmpty())
      pythonCandidates << pyEnv;
    QString condapath = QDir::homePath() + "/anaconda3/envs/eit/bin/python";
    if (QFileInfo::exists(condapath))
      pythonCandidates << condapath;
    pythonCandidates << "python3" << "python";

    bool convertedOk = false;
    QString pythonCmd = QString("import meshio,sys; mesh=meshio.read(sys.argv[1]); meshio.write(sys.argv[2], mesh)");
    for (const QString &pythonExec : pythonCandidates)
    {
      QStringList args;
      args << "-c" << pythonCmd << path << out;
      QProcess proc;
      proc.start(pythonExec, args);
      if (!proc.waitForFinished(20000))
      { // 20s timeout for python conversion
        proc.kill();
        proc.waitForFinished(2000);
      }
      if (proc.exitStatus() == QProcess::NormalExit && proc.exitCode() == 0)
      {
        convertedOk = true;
        break;
      }
      else
      {
        qDebug() << pythonExec << "meshio conversion failed (stderr):" << proc.readAllStandardError();
        qDebug() << pythonExec << "meshio conversion failed (stdout):" << proc.readAllStandardOutput();
      }
    }

    if (convertedOk)
    {
      localPath = out;
      ext = QFileInfo(localPath).suffix().toLower().toStdString();
    }
    else
    {
      // As last resort, try gmsh CLI conversion (may open GUI on some systems)
      QStringList gmshArgs;
      gmshArgs << "-noenv" << "-nopopup" << "-n" << "-format" << "vtk" << "-o" << out << path;
      QProcess gmshProc;
      gmshProc.start("gmsh", gmshArgs);
      bool started = gmshProc.waitForStarted(3000);
      if (started)
      {
        if (!gmshProc.waitForFinished(15000))
        {
          gmshProc.terminate();
          gmshProc.waitForFinished(2000);
          if (gmshProc.state() != QProcess::NotRunning)
            gmshProc.kill();
          gmshProc.waitForFinished(2000);
        }
      }
      if (gmshProc.exitStatus() == QProcess::NormalExit && gmshProc.exitCode() == 0 && QFileInfo::exists(out))
      {
        localPath = out;
        ext = QFileInfo(localPath).suffix().toLower().toStdString();
      }
      else
      {
        qDebug() << "gmsh conversion stderr:" << gmshProc.readAllStandardError();
        qDebug() << "gmsh conversion stdout:" << gmshProc.readAllStandardOutput();
        statusBar()->showMessage(tr("Failed to convert .msh to .vtk. Install meshio or set EIT_PYTHON to a Python with meshio."));
        return false;
      }
    }
  }

  if (ext == "vtu")
  {
    vtkNew<vtkXMLUnstructuredGridReader> r;
    r->SetFileName(localPath.toStdString().c_str());
    r->Update();
    data = r->GetOutput();
  }
  else
  {
    // try generic legacy reader
    vtkNew<vtkGenericDataObjectReader> r;
    r->SetFileName(localPath.toStdString().c_str());
    r->Update();
    if (r->IsFileUnstructuredGrid())
    {
      data = vtkUnstructuredGrid::SafeDownCast(r->GetOutput());
    }
    else if (r->IsFilePolyData())
    {
      data = vtkPolyData::SafeDownCast(r->GetOutput());
    }
    else
    {
      statusBar()->showMessage(tr("Unsupported file type or empty output"));
      return false;
    }
  }

  if (!data)
  {
    statusBar()->showMessage(tr("Failed to read mesh file"));
    return false;
  }
  // Ensure point data exists (convert cell data to point data if necessary)
  vtkSmartPointer<vtkDataSet> ds = data;
  vtkSmartPointer<vtkDataSet> pointDataSource = ds;

  // If there is cell data but no point data, convert
  if (ds->GetPointData()->GetNumberOfArrays() == 0 && ds->GetCellData()->GetNumberOfArrays() > 0)
  {
    vtkNew<vtkCellDataToPointData> c2p;
    c2p->SetInputData(ds);
    c2p->PassCellDataOn();
    c2p->Update();
    pointDataSource = vtkSmartPointer<vtkDataSet>::Take(c2p->GetOutput());
  }

  // Find a scalar array in point data (prefer the active scalars)
  vtkPointData* pd = pointDataSource->GetPointData();
  vtkDataArray* scalars = pd->GetScalars();
  if (!scalars && pd->GetNumberOfArrays() > 0)
    scalars = pd->GetArray(0);

  if (!scalars)
  {
    statusBar()->showMessage(tr("No point data arrays found (expecting conductivity in $NodeData)"));
    return false;
  }

  double range[2];
  scalars->GetRange(range);

  // Build lookup table
  lookupTable_ = vtkSmartPointer<vtkLookupTable>::New();
  lookupTable_->SetNumberOfTableValues(256);
  lookupTable_->SetRange(range);
  lookupTable_->Build();

  // Create surface polydata for rendering
  vtkSmartPointer<vtkPolyData> surfacePoly;
  if (vtkUnstructuredGrid::SafeDownCast(pointDataSource))
  {
    vtkNew<vtkDataSetSurfaceFilter> surf;
    surf->SetInputData(pointDataSource);
    surf->Update();
    surfacePoly = surf->GetOutput();
  }
  else if (vtkPolyData::SafeDownCast(pointDataSource))
  {
    surfacePoly = vtkPolyData::SafeDownCast(pointDataSource);
  }
  else
  {
    statusBar()->showMessage(tr("Unsupported dataset type for surface extraction"));
    return false;
  }

  // Setup mapper & actor
  surfaceMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
  surfaceMapper_->SetInputData(surfacePoly);
  surfaceMapper_->SetLookupTable(lookupTable_);
  surfaceMapper_->SetScalarModeToUsePointData();
  surfaceMapper_->ScalarVisibilityOn();
  surfaceMapper_->SetColorModeToMapScalars();
  surfaceMapper_->SetScalarRange(range);

  surfaceActor_ = vtkSmartPointer<vtkActor>::New();
  surfaceActor_->SetMapper(surfaceMapper_);

  // Add to renderer
  if (!renderer_)
    renderer_ = vtkSmartPointer<vtkRenderer>::New();
  renderer_->RemoveAllViewProps();
  renderer_->AddActor(surfaceActor_);

  // Scalar bar
  scalarBar_ = vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar_->SetLookupTable(lookupTable_);
  scalarBar_->SetTitle(scalars->GetName() ? scalars->GetName() : "Scalar");
  vtkTextProperty* tp = scalarBar_->GetLabelTextProperty();
  tp->SetFontSize(12);
  renderer_->AddActor2D(scalarBar_);

  // Store bounds and choose largest axis for slicing
  pointDataSource->GetBounds(dataBounds_);
  double dx = dataBounds_[1] - dataBounds_[0];
  double dy = dataBounds_[3] - dataBounds_[2];
  double dz = dataBounds_[5] - dataBounds_[4];
  if (dx >= dy && dx >= dz) sliceAxis_ = 0;
  else if (dy >= dx && dy >= dz) sliceAxis_ = 1;
  else sliceAxis_ = 2;

  // Create slice plane and cutter
  slicePlane_ = vtkSmartPointer<vtkPlane>::New();
  double normal[3] = {0.0, 0.0, 0.0};
  normal[sliceAxis_] = 1.0;
  slicePlane_->SetNormal(normal);
  // initial position at slider value
  double pos0 = dataBounds_[2*sliceAxis_] + 0.5*(dataBounds_[2*sliceAxis_+1] - dataBounds_[2*sliceAxis_]);
  double origin0[3] = { (dataBounds_[0]+dataBounds_[1])/2.0, (dataBounds_[2]+dataBounds_[3])/2.0, (dataBounds_[4]+dataBounds_[5])/2.0 };
  origin0[sliceAxis_] = pos0;
  slicePlane_->SetOrigin(origin0);

  sliceCutter_ = vtkSmartPointer<vtkCutter>::New();
  sliceCutter_->SetCutFunction(slicePlane_.Get());
  sliceCutter_->SetInputData(pointDataSource);
  sliceCutter_->Update();

  sliceMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
  sliceMapper_->SetInputData(sliceCutter_->GetOutput());
  sliceMapper_->SetLookupTable(lookupTable_);
  sliceMapper_->SetScalarModeToUsePointData();
  sliceMapper_->ScalarVisibilityOn();
  sliceMapper_->SetScalarRange(range);

  sliceActor_ = vtkSmartPointer<vtkActor>::New();
  sliceActor_->SetMapper(sliceMapper_);
  sliceActor_->GetProperty()->SetLineWidth(2.0);
  renderer_->AddActor(sliceActor_);

  // Attach renderer to widget's render window
  if (vtkWidget_ && vtkWidget_->renderWindow())
  {
    // Remove previous renderers to avoid duplicates
    vtkRenderWindow* rw = vtkWidget_->renderWindow();
    // ensure renderer is present (remove any previous instance then add)
    rw->RemoveRenderer(renderer_.Get());
    rw->AddRenderer(renderer_.Get());
    rw->Render();
  }

  statusBar()->showMessage(tr("Loaded mesh: %1").arg(path));
  currentMeshPath_ = path;

  return true;
}

ConductivityViewer::ConductivityViewer(QWidget* parent, const QString& initialMesh)
  : QMainWindow(parent)
{
  setupUi();
  if (!initialMesh.isEmpty())
    loadMeshFile(initialMesh);
}

void ConductivityViewer::setupUi()
{
  // central widget and layout
  QWidget* central = new QWidget(this);
  setCentralWidget(central);
  QHBoxLayout* mainLayout = new QHBoxLayout(central);

  // left controls
  QWidget* left = new QWidget(central);
  QVBoxLayout* leftLayout = new QVBoxLayout(left);
  left->setLayout(leftLayout);

  renderModeCombo_ = new QComboBox(left);
  renderModeCombo_->addItem(tr("Surface"));
  renderModeCombo_->addItem(tr("Slice"));
  renderModeCombo_->addItem(tr("Clipped"));
  leftLayout->addWidget(new QLabel(tr("Render mode:"), left));
  leftLayout->addWidget(renderModeCombo_);

  sliceSlider_ = new QSlider(Qt::Horizontal, left);
  sliceSlider_->setRange(0, 100);
  sliceSlider_->setValue(50);
  leftLayout->addWidget(new QLabel(tr("Slice:"), left));
  leftLayout->addWidget(sliceSlider_);

  opacitySlider_ = new QSlider(Qt::Horizontal, left);
  opacitySlider_->setRange(0, 100);
  opacitySlider_->setValue(100);
  leftLayout->addWidget(new QLabel(tr("Opacity:"), left));
  leftLayout->addWidget(opacitySlider_);

  leftLayout->addStretch(1);

  // VTK widget
  vtkWidget_ = new QVTKOpenGLNativeWidget(central);

  // Make window reasonably large by default
  resize(1200, 800);
  setMinimumSize(800, 600);

  mainLayout->addWidget(left);
  mainLayout->addWidget(vtkWidget_, 1);

  // Menu
  QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
  QAction* openAct = new QAction(tr("Open..."), this);
  fileMenu->addAction(openAct);
  connect(openAct, &QAction::triggered, this, &ConductivityViewer::openMesh);

  // Connect controls
  connect(sliceSlider_, &QSlider::valueChanged, this, &ConductivityViewer::onSliceSliderChanged);
  connect(renderModeCombo_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &ConductivityViewer::onRenderModeChanged);
  connect(opacitySlider_, &QSlider::valueChanged, this, &ConductivityViewer::onOpacityChanged);

  // status bar
  statusBar()->showMessage(tr("Ready"));
}

void ConductivityViewer::openMesh()
{
  QString path = QFileDialog::getOpenFileName(this, tr("Open mesh"), QString(), tr("Meshes (*.msh *.vtu *.vtk);;All files (*)"));
  if (!path.isEmpty())
    loadMeshFile(path);
}

void ConductivityViewer::onSliceSliderChanged(int value)
{
  if (!slicePlane_ || !sliceCutter_)
  {
    statusBar()->showMessage(tr("No slice available"));
    return;
  }

  double t = value / 100.0;
  double minv = dataBounds_[2*sliceAxis_];
  double maxv = dataBounds_[2*sliceAxis_+1];
  double pos = minv + t * (maxv - minv);

  double origin[3] = { (dataBounds_[0]+dataBounds_[1])/2.0, (dataBounds_[2]+dataBounds_[3])/2.0, (dataBounds_[4]+dataBounds_[5])/2.0 };
  origin[sliceAxis_] = pos;
  slicePlane_->SetOrigin(origin);
  sliceCutter_->Update();

  if (vtkWidget_ && vtkWidget_->renderWindow())
    vtkWidget_->renderWindow()->Render();

  statusBar()->showMessage(tr("Slice moved to %1%").arg(value));
}

void ConductivityViewer::onRenderModeChanged(int index)
{
  Q_UNUSED(index);
  // Placeholder: mode switching (surface/clip/probe) can be implemented
  statusBar()->showMessage(tr("Render mode changed"));
}

void ConductivityViewer::onOpacityChanged(int value)
{
  double op = value / 100.0;
  if (surfaceActor_)
    surfaceActor_->GetProperty()->SetOpacity(op);
  if (vtkWidget_ && vtkWidget_->renderWindow())
    vtkWidget_->renderWindow()->Render();
}

void ConductivityViewer::resetCamera()
{
  if (renderer_)
  {
    renderer_->ResetCamera();
    if (vtkWidget_ && vtkWidget_->renderWindow())
      vtkWidget_->renderWindow()->Render();
  }
}