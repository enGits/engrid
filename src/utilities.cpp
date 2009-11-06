#include "utilities.h"

#include <cmath>

#include <QApplication>
#include <QPushButton>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QtDebug>
#include <QMainWindow>
#include <QHBoxLayout>
#include <QString>
#include <QStringList>
#include <QDir>
#include <QtDebug>
#include <QFileInfo>
// #include <QFileDialogArgs>

#include "error.h"
#include "engrid.h"
#include "math/mathvector.h"
#include "math/smallsquarematrix.h"
#include "math/linsolve.h"

double toDouble(QString str)
{
  str.replace(QString(QObject::tr(",")), QString(QObject::tr(".")));
  return(str.toDouble());
}

QString toString(double x, QString separator)
{
  QString ret;
  ret.setNum(x);
  ret.replace(QString(QObject::tr(".")), separator);
  return(ret);
}

QString addSuffix(QString filename, QString suffix, bool remove_old_suffix)
{
  QFileInfo fileinfo(filename);
  if (remove_old_suffix) {
    return fileinfo.completeBaseName() + QObject::tr(".") + suffix;
  }
  else {
    if (fileinfo.suffix().toLower() == suffix) {
      return fileinfo.absoluteFilePath();
    }
    else {
      return fileinfo.absoluteFilePath() + QObject::tr(".") + suffix;
    }
  }
}

Qt::CheckState int2CheckState(int a)
{
  if(a==0) return(Qt::Unchecked);
  if(a==1) return(Qt::PartiallyChecked);
  else return(Qt::Checked);
}

int CheckState2int(Qt::CheckState a)
{
  if(a==Qt::Unchecked) return(0);
  if(a==Qt::PartiallyChecked) return(1);
  else return(2);
}

QString vector2csv(QVector <double> vector)
{
  QString csv;
  for (int i = 0; i < vector.size(); i++) {
    csv += QString::number(vector[i]);
    if (i < vector.size() - 1) csv += QObject::tr(";");
  }
  return csv;
}

QVector <double> csv2vector(QString csv)
{
  QVector <double> vector;
  if (csv.isEmpty()) return vector;
  QStringList list = csv.split(QObject::tr(";"));
  foreach(QString str, list) {
    vector.push_back(str.toDouble());
  }
  return vector;
}

mat3_t rotationMatrix_X(double a_rad)
{
  mat3_t Rx;
  Rx[0][0] = 1; Rx[0][1] = 0;           Rx[0][2] = 0;
  Rx[1][0] = 0; Rx[1][1] = cos(a_rad);  Rx[1][2] = sin(a_rad);
  Rx[2][0] = 0; Rx[2][1] = -sin(a_rad); Rx[2][2] = cos(a_rad);
  return Rx;
}

mat3_t rotationMatrix_Y(double a_rad)
{
  mat3_t Ry;
  Ry[0][0] = cos(a_rad); Ry[0][1] = 0; Ry[0][2] = -sin(a_rad);
  Ry[1][0] = 0;          Ry[1][1] = 1; Ry[1][2] = 0;
  Ry[2][0] = sin(a_rad); Ry[2][1] = 0; Ry[2][2] = cos(a_rad);
  return Ry;
}

mat3_t rotationMatrix_Z(double a_rad)
{
  mat3_t Rz;
  Rz[0][0] = cos(a_rad);  Rz[0][1] = sin(a_rad); Rz[0][2] = 0;
  Rz[1][0] = -sin(a_rad); Rz[1][1] = cos(a_rad); Rz[1][2] = 0;
  Rz[2][0] = 0;           Rz[2][1] = 0;          Rz[2][2] = 1;
  return Rz;
}

mat3_t rotateRelativeZXZ(double angle_1_rad, double angle_2_rad, double angle_3_rad)
{
  return rotationMatrix_Z(angle_3_rad)*rotationMatrix_X(angle_2_rad)*rotationMatrix_Z(angle_1_rad);
}

mat3_t rotateAbsoluteZXY(double angle_1_rad, double angle_2_rad, double angle_3_rad)
{
  /*  cout<<"rotationMatrix_Z(angle_1_rad)="<<rotationMatrix_Z(angle_1_rad)<<endl;
    cout<<"rotationMatrix_X(angle_2_rad)="<<rotationMatrix_X(angle_2_rad)<<endl;
    cout<<"rotationMatrix_Y(angle_3_rad)="<<rotationMatrix_Y(angle_3_rad)<<endl;*/
  return rotationMatrix_Z(angle_1_rad)*rotationMatrix_X(angle_2_rad)*rotationMatrix_Y(angle_3_rad);
}

double getGamma(vec3_t V)
{
  return V[0] == 0.0 && V[1] == 0.0 && V[2] == 0.0 ? 0.0 : atan2(sqrt(V[0]*V[0] + V[1]*V[1]), V[2]);
}

double getPhi(vec3_t V)
{
  return V[0] == 0.0 && V[1] == 0.0 ? 0.0 : atan2(V[1], V[0]);
}

QString getDirectory(QWidget * parent, const QString & caption, const QString & selected)
{
  qDebug() << "selected=" << selected;
//   QFileDialog filedialog;
//   filedialog.selectFile(selected);
  QFileInfo fileinfo(selected);
  QString dir = fileinfo.absolutePath();
  qDebug() << "dir=" << dir;

  /*  QFileDialogArgs args;
    args.parent = parent;
    args.caption = caption;
    args.directory = QFileDialogPrivate::workingDirectory(dir);
    args.mode = (options & ShowDirsOnly ? DirectoryOnly : Directory);
    args.options = options;*/

  // create a qt dialog
  QFileDialog dialog(parent, caption, dir);//(args);
  dialog.setFileMode(QFileDialog::Directory);
//   dialog.setFileMode(QFileDialog::DirectoryOnly);
  dialog.setOption(QFileDialog::ShowDirsOnly, true);
  /*  args.parent = parent;
    args.caption = caption;
    args.directory = QFileDialogPrivate::workingDirectory(dir);*/
//   args.mode = (options & ShowDirsOnly ? DirectoryOnly : Directory);
//   args.options = options;

  dialog.selectFile(selected);
  if (dialog.exec() == QDialog::Accepted) {
    return dialog.selectedFiles().value(0);
  }

  return QString();
//   return filedialog.getExistingDirectory (parent, caption, dir);
}
