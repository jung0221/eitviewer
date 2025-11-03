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
#include <iostream>

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

#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>

#include <QVTKOpenGLNativeWidget.h>
#include <vtkIntArray.h>
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

  // Inline parser for Gmsh .msh files: reads $Nodes, $Elements and $NodeData and
  // constructs a vtkUnstructuredGrid with a point array named "conductivity".
  auto parseMshFile = [](const QString &mshPath) -> vtkSmartPointer<vtkUnstructuredGrid>
  {
    std::ifstream ifs(mshPath.toStdString());
    if (!ifs)
    {
      qDebug() << "Failed to open" << mshPath;
      return nullptr;
    }

    // helper: trim whitespace (including CR/LF) from both ends of a line
    auto trim = [&](std::string &s)
    {
      while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))
        s.pop_back();
      while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())))
        s.erase(0, 1);
    };

    vtkNew<vtkPoints> points;
    vtkNew<vtkUnstructuredGrid> ug;
    vtkNew<vtkIntArray> idArr;
    idArr->SetName("gmsh_node_id");
    idArr->SetNumberOfComponents(1);

    std::unordered_map<int, vtkIdType> nodeIdToIndex;
    std::vector<int> nodeOrder; // preserves the order nodes were defined

    struct ElemRecord
    {
      int type;
      std::vector<int> nodeIds;
    };
    std::vector<ElemRecord> elements;

    std::string line;
    while (std::getline(ifs, line))
    {
      trim(line);
      if (line.size() == 0)
        continue;
      if (line == "$Nodes")
      {
        // read number of nodes
        if (!std::getline(ifs, line))
          break;
        trim(line);
        int n = std::stoi(line);
        for (int i = 0; i < n; ++i)
        {
          if (!std::getline(ifs, line))
            break;
          trim(line);
          if (line.size() == 0)
          {
            --i;
            continue;
          }
          std::istringstream iss(line);
          int id;
          double x, y, z;
          iss >> id >> x >> y >> z;
          vtkIdType idx = points->InsertNextPoint(x, y, z);
          nodeIdToIndex[id] = idx;
          nodeOrder.push_back(id);
          idArr->InsertNextValue(id);
        }
      }
      else if (line == "$Elements")
      {
        if (!std::getline(ifs, line))
          break;
        trim(line);
        int ne = std::stoi(line);
        for (int i = 0; i < ne; ++i)
        {
          if (!std::getline(ifs, line))
            break;
          trim(line);
          if (line.size() == 0)
          {
            --i;
            continue;
          }
          std::istringstream iss(line);
          int elemId, elemType, numTags;
          iss >> elemId >> elemType >> numTags;
          for (int t = 0; t < numTags; ++t)
          {
            int tmp;
            iss >> tmp;
          }
          std::vector<int> nids;
          int nid;
          while (iss >> nid)
            nids.push_back(nid);
          elements.push_back({elemType, nids});
        }
      }
      else if (line == "$NodeData")
      {
        // Parse tags according to Gmsh MSH specification
        auto readNonEmpty = [&](std::string &out) -> bool
        {
          while (std::getline(ifs, out))
          {
            trim(out);
            if (!out.empty())
              return true;
          }
          return false;
        };

        std::string s;
        if (!readNonEmpty(s))
          break;
        int numStringTags = std::stoi(s);
        for (int i = 0; i < numStringTags; ++i)
        {
          std::getline(ifs, s);
        }
        if (!readNonEmpty(s))
          break;
        int numRealTags = std::stoi(s);
        for (int i = 0; i < numRealTags; ++i)
        {
          std::getline(ifs, s);
        }
        if (!readNonEmpty(s))
          break;
        int numIntTags = std::stoi(s);
        std::vector<int> intTags;
        for (int i = 0; i < numIntTags; ++i)
        {
          if (!readNonEmpty(s))
            break;
          intTags.push_back(std::stoi(s));
        }

        // Now, read data lines until $EndNodeData
        std::vector<double> orderedValues;
        std::vector<std::pair<int, double>> pairs;
        while (std::getline(ifs, s))
        {
          if (s == "$EndNodeData")
            break;
          if (s.empty())
            continue;
          std::istringstream iss(s);
          std::vector<std::string> toks;
          std::string tok;
          while (iss >> tok)
            toks.push_back(tok);
          if (toks.size() >= 2)
          {
            // try nodeID value pair
            try
            {
              int nid = std::stoi(toks[0]);
              double val = std::stod(toks[1]);
              pairs.emplace_back(nid, val);
            }
            catch (...)
            {
              try
              {
                orderedValues.push_back(std::stod(toks[0]));
              }
              catch (...)
              {
              }
            }
          }
          else if (toks.size() == 1)
          {
            try
            {
              orderedValues.push_back(std::stod(toks[0]));
            }
            catch (...)
            {
            }
          }
        }

        // Build conductivity array
        vtkNew<vtkDoubleArray> cond;
        cond->SetName("conductivity");
        cond->SetNumberOfComponents(1);
        vtkIdType npoints = points->GetNumberOfPoints();
        cond->SetNumberOfTuples(npoints);
        double nanv = std::numeric_limits<double>::quiet_NaN();
        for (vtkIdType i = 0; i < npoints; ++i)
          cond->SetTuple1(i, nanv);

        if (!pairs.empty())
        {
          for (auto &pr : pairs)
          {
            int nid = pr.first;
            double val = pr.second;
            auto it = nodeIdToIndex.find(nid);
            if (it != nodeIdToIndex.end())
              cond->SetTuple1(it->second, val);
            else
              qDebug() << "Node id" << nid << "in $NodeData not present in $Nodes";
          }
        }
        else if (!orderedValues.empty() && (vtkIdType)orderedValues.size() == npoints)
        {
          for (vtkIdType i = 0; i < npoints; ++i)
            cond->SetTuple1(i, orderedValues[(size_t)i]);
        }
        else if (!orderedValues.empty())
        {
          qDebug() << "NodeData values count" << (int)orderedValues.size() << "does not match number of nodes" << (int)npoints;
        }

        ug->GetPointData()->AddArray(cond);
        // also expose original gmsh node ids if available
        if (idArr->GetNumberOfTuples() == points->GetNumberOfPoints())
          ug->GetPointData()->AddArray(idArr);
        ug->GetPointData()->SetScalars(cond);
      }
    }

    // Insert elements as cells into unstructured grid
    for (const auto &er : elements)
    {
      int gmshType = er.type;
      int vtkType = -1;
      switch (gmshType)
      {
      case 1:
        vtkType = VTK_LINE;
        break;
      case 2:
        vtkType = VTK_TRIANGLE;
        break;
      case 3:
        vtkType = VTK_QUAD;
        break;
      case 4:
        vtkType = VTK_TETRA;
        break;
      case 5:
        vtkType = VTK_HEXAHEDRON;
        break;
      case 6:
        vtkType = VTK_WEDGE;
        break;
      case 7:
        vtkType = VTK_PYRAMID;
        break;
      default:
        vtkType = -1;
        break;
      }
      if (vtkType < 0)
        continue; // skip unsupported types

      vtkNew<vtkIdList> idList;
      bool ok = true;
      for (int nid : er.nodeIds)
      {
        auto it = nodeIdToIndex.find(nid);
        if (it == nodeIdToIndex.end())
        {
          ok = false;
          break;
        }
        idList->InsertNextId(it->second);
      }
      if (!ok)
        continue;
      ug->InsertNextCell(vtkType, idList);
    }

    ug->SetPoints(points);
    return ug;
  };
  // Handle Gmsh .msh: first try converting with gmsh CLI to a .vtk file (no Python required).
  QString localPath = path;
  if (ext == "msh")
  {
    QString tmpdir = QDir::tempPath();
    QString base = QFileInfo(path).baseName();
    // creates a temporary mesh with vtu format (converted from .msh)
    QString out = tmpdir + "/eitviewer_" + base + "_" + QString::number(QDateTime::currentSecsSinceEpoch()) + ".vtu";
    // Try parsing the .msh directly (pure C++). If successful, we will skip
    // external conversion and use the in-memory unstructured grid.
    vtkSmartPointer<vtkUnstructuredGrid> parsedGrid = parseMshFile(path);

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

  if (!data && ext == "vtu")
  {
    vtkNew<vtkXMLUnstructuredGridReader> r;
    r->SetFileName(localPath.toStdString().c_str());
    r->Update();
    data = r->GetOutput();
  }
  else
  {
    // try generic legacy reader
    if (!data)
    {
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
  vtkPointData *pd = pointDataSource->GetPointData();
  vtkDataArray *scalars = pd->GetScalars();
  if (!scalars && pd->GetNumberOfArrays() > 0)
    scalars = pd->GetArray(0);

  if (!scalars)
  {
    statusBar()->showMessage(tr("No point data arrays found (expecting conductivity in $NodeData)"));
    return false;
  }

  // Ensure array has a name and make it active on the source dataset (and on surface after extraction)
  const char *arrNameC = scalars->GetName();
  std::string arrNameStd;
  if (!arrNameC)
  {
    arrNameStd = "conductivity";
    scalars->SetName(arrNameStd.c_str());
    arrNameC = scalars->GetName();
  }
  else
  {
    arrNameStd = arrNameC;
  }
  qDebug() << "Using scalar array:" << arrNameC;
  pd->SetActiveScalars(arrNameC);

  double range[2];
  scalars->GetRange(range);
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
    // Ensure point scalars are present on the generated surface polydata.
    // Use vtkProbeFilter to sample the original pointDataSource onto the surface
    // geometry — this attaches the exact point-data arrays (by interpolation or
    // direct mapping) to the surface points so the mapper can use them.
    vtkNew<vtkProbeFilter> surfaceProbe;
    surfaceProbe->SetInputData(surfacePoly);
    surfaceProbe->SetSourceData(pointDataSource);
    surfaceProbe->Update();
    surfacePoly = vtkPolyData::SafeDownCast(surfaceProbe->GetOutput());
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
  // Use point field data by name to ensure correct array is mapped
  surfaceMapper_->SetLookupTable(lookupTable_);
  surfaceMapper_->SetScalarModeToUsePointFieldData();
  surfaceMapper_->SelectColorArray(arrNameStd.c_str());
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
  vtkTextProperty *tp = scalarBar_->GetLabelTextProperty();
  tp->SetFontSize(12);
  renderer_->AddActor2D(scalarBar_);

  // Store bounds and choose largest axis for slicing
  pointDataSource->GetBounds(dataBounds_);
  double dx = dataBounds_[1] - dataBounds_[0];
  double dy = dataBounds_[3] - dataBounds_[2];
  double dz = dataBounds_[5] - dataBounds_[4];
  if (dx >= dy && dx >= dz)
    sliceAxis_ = 0;
  else if (dy >= dx && dy >= dz)
    sliceAxis_ = 1;
  else
    sliceAxis_ = 2;

  // Create slice plane and cutter
  slicePlane_ = vtkSmartPointer<vtkPlane>::New();
  double normal[3] = {0.0, 0.0, 0.0};
  normal[sliceAxis_] = 1.0;
  slicePlane_->SetNormal(normal);
  // initial position at slider value
  double pos0 = dataBounds_[2 * sliceAxis_] + 0.5 * (dataBounds_[2 * sliceAxis_ + 1] - dataBounds_[2 * sliceAxis_]);
  double origin0[3] = {(dataBounds_[0] + dataBounds_[1]) / 2.0, (dataBounds_[2] + dataBounds_[3]) / 2.0, (dataBounds_[4] + dataBounds_[5]) / 2.0};
  origin0[sliceAxis_] = pos0;
  slicePlane_->SetOrigin(origin0);

  sliceCutter_ = vtkSmartPointer<vtkCutter>::New();
  sliceCutter_->SetCutFunction(slicePlane_.Get());
  sliceCutter_->SetInputData(pointDataSource);
  sliceCutter_->Update();

  sliceMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
  // Probe cutter output with original pointDataSource to ensure scalars are
  // present on the slice polydata (some pipeline steps may lose point arrays).
  // Create a persistent probe filter that samples the original point data onto
  // the slice geometry. Use SetInputConnection so the probe follows the cutter
  // output when the slice plane origin changes.
  sliceProbe_ = vtkSmartPointer<vtkProbeFilter>::New();
  sliceProbe_->SetInputConnection(sliceCutter_->GetOutputPort());
  sliceProbe_->SetSourceData(pointDataSource);
  sliceMapper_->SetInputConnection(sliceProbe_->GetOutputPort());
  sliceMapper_->SetLookupTable(lookupTable_);
  sliceMapper_->SetScalarModeToUsePointFieldData();
  sliceMapper_->SelectColorArray(arrNameStd.c_str());
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
    vtkRenderWindow *rw = vtkWidget_->renderWindow();
    // ensure renderer is present (remove any previous instance then add)
    rw->RemoveRenderer(renderer_.Get());
    rw->AddRenderer(renderer_.Get());
    rw->Render();
  }

  statusBar()->showMessage(tr("Loaded mesh: %1").arg(path));
  currentMeshPath_ = path;

  return true;
}

ConductivityViewer::ConductivityViewer(QWidget *parent, const QString &initialMesh)
    : QMainWindow(parent)
{
  setupUi();
  if (!initialMesh.isEmpty())
    loadMeshFile(initialMesh);
}

void ConductivityViewer::setupUi()
{
  // central widget and layout
  QWidget *central = new QWidget(this);
  setCentralWidget(central);
  QHBoxLayout *mainLayout = new QHBoxLayout(central);

  // left controls
  QWidget *left = new QWidget(central);
  QVBoxLayout *leftLayout = new QVBoxLayout(left);
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
  QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
  QAction *openAct = new QAction(tr("Open..."), this);
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
  double minv = dataBounds_[2 * sliceAxis_];
  double maxv = dataBounds_[2 * sliceAxis_ + 1];
  double pos = minv + t * (maxv - minv);

  double origin[3] = {(dataBounds_[0] + dataBounds_[1]) / 2.0, (dataBounds_[2] + dataBounds_[3]) / 2.0, (dataBounds_[4] + dataBounds_[5]) / 2.0};
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