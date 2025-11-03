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

  parser.process(app);

  QString meshPath;
  if (parser.isSet(meshOpt)) {
    meshPath = parser.value(meshOpt);
  }

  ConductivityViewer w(nullptr, meshPath);
  w.show();

  return app.exec();
}
