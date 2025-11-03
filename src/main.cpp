#include "ConductivityViewer.h"

#include <QApplication>
#include <QCommandLineParser>

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  QCoreApplication::setApplicationName("eitviewer");
  QCoreApplication::setApplicationVersion("0.1");

  QCommandLineParser parser;
  parser.setApplicationDescription("Conductivity mesh viewer");
  parser.addHelpOption();
  parser.addVersionOption();

  QCommandLineOption meshOpt(QStringList() << "m" << "mesh",
                             "Initial mesh file to open.",
                             "meshfile");
  parser.addOption(meshOpt);

  QCommandLineOption targetOpt(QStringList() << "t" << "target",
                               "Initial target mesh file to open (left).",
                               "meshfile");
  parser.addOption(targetOpt);

  QCommandLineOption predictedOpt(QStringList() << "p" << "predicted",
                                  "Initial predicted mesh file to open (right).",
                                  "meshfile");
  parser.addOption(predictedOpt);

  parser.process(app);

  QString meshPath;
  if (parser.isSet(meshOpt)) {
    meshPath = parser.value(meshOpt);
  }
  QString targetPath;
  if (parser.isSet(targetOpt)) {
    targetPath = parser.value(targetOpt);
  }
  QString predictedPath;
  if (parser.isSet(predictedOpt)) {
    predictedPath = parser.value(predictedOpt);
  }

  // Prefer explicit target/predicted options; fall back to single -m for left view
  QString leftInit = targetPath.isEmpty() ? meshPath : targetPath;
  ConductivityViewer w(nullptr, leftInit, predictedPath);
  w.show();

  return app.exec();
}
