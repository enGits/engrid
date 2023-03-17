// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

#include "vtkIntArray.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "egvtkobject.h"
#include "vtkCell.h"

double toDouble(QString str) {
  str.replace(QString(QObject::tr(",")), QString(QObject::tr(".")));
  return(str.toDouble());
}

QString toString(double x, QString separator) {
  QString ret;
  ret.setNum(x);
  ret.replace(QString(QObject::tr(".")), separator);
  return(ret);
}

QString addSuffix(QString filename, QString suffix, bool remove_old_suffix) {
  QFileInfo fileinfo(filename);
  if (remove_old_suffix) {
    return fileinfo.completeBaseName() + QObject::tr(".") + suffix;
  } else {
    if (fileinfo.suffix().toLower() == suffix) {
      return fileinfo.absoluteFilePath();
    } else {
      return fileinfo.absoluteFilePath() + QObject::tr(".") + suffix;
    }
  }
}

Qt::CheckState int2CheckState(int a) {
  if (a == 0) return(Qt::Unchecked);
  if (a == 1) return(Qt::PartiallyChecked);
  else return(Qt::Checked);
}

int CheckState2int(Qt::CheckState a) {
  if (a == Qt::Unchecked) return(0);
  if (a == Qt::PartiallyChecked) return(1);
  else return(2);
}

QString vector2csv(QVector <double> vector) {
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
  return rotationMatrix_Z(angle_1_rad)*rotationMatrix_X(angle_2_rad)*rotationMatrix_Y(angle_3_rad);
}

double getGamma(vec3_t V) {
  return V[0] == 0.0 && V[1] == 0.0 && V[2] == 0.0 ? 0.0 : atan2(sqrt(V[0]*V[0] + V[1]*V[1]), V[2]);
}

double getPhi(vec3_t V) {
  return V[0] == 0.0 && V[1] == 0.0 ? 0.0 : atan2(V[1], V[0]);
}

QString getDirectory(QWidget * parent, const QString & caption, const QString & selected)
{
  qDebug() << "selected=" << selected;
  QFileInfo fileinfo(selected);
  QString dir = fileinfo.absolutePath();
  qDebug() << "dir=" << dir;

  // create a qt dialog
  QFileDialog dialog(parent, caption, dir);//(args);
  dialog.setFileMode(QFileDialog::Directory);
#if QT_VERSION >= 0x040500
  dialog.setOption(QFileDialog::ShowDirsOnly, true);
#endif
  dialog.selectFile(selected);

  if (dialog.exec() == QDialog::Accepted) {
    return dialog.selectedFiles().value(0);
  }
  return QString();
}

int cout_grid(ostream &stream, vtkUnstructuredGrid *grid, bool npoints, bool ncells, bool points, bool cells) {
  stream << "=============" << endl;
  if (npoints) stream << "grid->GetNumberOfPoints()=" << grid->GetNumberOfPoints() << endl;
  if (ncells) stream << "grid->GetNumberOfCells()=" << grid->GetNumberOfCells() << endl;
  if (points) {
    for (vtkIdType i = 0; i < grid->GetNumberOfPoints(); ++i) {
      vec3_t x;
      grid->GetPoint(i, x.data());
      stream << "Vertex " << i << " = " << x << endl;
    }
  }
  if (cells) {
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
      vtkCell *C = (vtkCell *) vtkCell::New();
      C = grid->GetCell(i);
      EG_GET_CELL(i, grid);
      stream << "Cell " << i << " = ";
      for (int j = 0; j < num_pts; j++)
        stream << pts[j] << " ";
      stream << "boundary_code=" << cell_code->GetValue(i);
      stream << endl;
    }
  }
  stream << "=============" << endl;
  return 0;
}

///////////////////////////////////////////
//Warning: Untested
int addCell(vtkUnstructuredGrid* a_grid, vtkIdType A, vtkIdType B, vtkIdType C, int bc) {
  vtkIdType npts = 3;
  vtkIdType pts[3];
  pts[0] = A;
  pts[1] = B;
  pts[2] = C;
  vtkIdType newCellId = a_grid->InsertNextCell(VTK_TRIANGLE, npts, pts);
  EG_VTKDCC(vtkIntArray, cell_code, a_grid, "cell_code");
  cell_code->SetValue(newCellId, bc);
  return(0);
}

///////////////////////////////////////////

int getShortestSide(vtkIdType a_id_cell, vtkUnstructuredGrid* a_grid) {
  EG_GET_CELL(a_id_cell, a_grid);
  vec3_t* x = new vec3_t[num_pts];
  for (int i = 0; i < num_pts; i++) a_grid->GetPoints()->GetPoint(pts[i], x[i].data());
  int id_minlen = 0;
  double minlen = (x[1] - x[0]).abs();
  for (int i = 1; i < num_pts; i++) {
    double len = (x[(i+1)%num_pts] - x[i]).abs();
    if (len < minlen) {
      minlen = len;
      id_minlen = i;
    }
  }
  delete x;
  return(id_minlen);
}

int getLongestSide(vtkIdType a_id_cell, vtkUnstructuredGrid* a_grid) {
  EG_GET_CELL(a_id_cell, a_grid);
  vec3_t* x = new vec3_t[num_pts];
  for (int i = 0; i < num_pts; i++) a_grid->GetPoints()->GetPoint(pts[i], x[i].data());
  int id_maxlen = 0;
  double maxlen = (x[1] - x[0]).abs();
  cout << "maxlen=" << maxlen << endl;
  for (int i = 1; i < num_pts; i++) {
    double len = (x[(i+1)%num_pts] - x[i]).abs();
    cout << "len[" << i << "]=" << len << endl;
    if (len > maxlen) {
      maxlen = len;
      id_maxlen = i;
    }
  }
  delete x;
  return(id_maxlen);
}

int getSide(vtkIdType a_id_cell, vtkUnstructuredGrid* a_grid, vtkIdType a_id_node1, vtkIdType a_id_node2) {
  EG_GET_CELL(a_id_cell, a_grid);
  QVector <vtkIdType> edge(2);

  int n = 0;
  for (int i = 0; i < num_pts; i++) {
    if (pts[i] == a_id_node1) { edge[0] = i; n++;}
    if (pts[i] == a_id_node2) { edge[1] = i; n++;}
  }
  if (n != 2) {
    EG_BUG;
    return(-1);
  }
  qSort(edge.begin(), edge.end());
  if (edge[0] == 0 && edge[1] == num_pts - 1) return(num_pts - 1);
  else return(edge[0]);
}
///////////////////////////////////////////

QString cell2str(vtkIdType id_cell, vtkUnstructuredGrid* grid) {
  QString tmp;
  tmp.setNum(id_cell);
  QString txt = "id_cell=" + tmp;

  EG_GET_CELL(id_cell, grid);

  txt += " [";
  for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
    tmp.setNum(pts[i_pts]);
    txt += tmp;
    if (i_pts < num_pts - 1) {
      txt += ",";
    }
  }
  txt += "]";
  return(txt);
}



///////////////////////////////////////////

pair<vtkIdType, vtkIdType> OrderedPair(vtkIdType a, vtkIdType b) {
  vtkIdType x = min(a, b);
  vtkIdType y = max(a, b);
  return(pair<vtkIdType, vtkIdType>(x, y));
}

const char* VertexType2Str(char T) {
  if (T == EG_SIMPLE_VERTEX)         return("EG_SIMPLE_VERTEX");
  if (T == EG_FIXED_VERTEX)          return("EG_FIXED_VERTEX");
  if (T == EG_FEATURE_EDGE_VERTEX)   return("EG_FEATURE_EDGE_VERTEX");
  if (T == EG_FEATURE_CORNER_VERTEX) return("EG_FEATURE_CORNER_VERTEX");
  if (T == EG_BOUNDARY_EDGE_VERTEX)  return("EG_BOUNDARY_EDGE_VERTEX");
  else return("Unknown vertex type");
}

char Str2VertexType(QString S) {
  if (S == "EG_SIMPLE_VERTEX")         return(EG_SIMPLE_VERTEX);
  if (S == "EG_FIXED_VERTEX")          return(EG_FIXED_VERTEX);
  if (S == "EG_FEATURE_EDGE_VERTEX")   return(EG_FEATURE_EDGE_VERTEX);
  if (S == "EG_FEATURE_CORNER_VERTEX") return(EG_FEATURE_CORNER_VERTEX);
  if (S == "EG_BOUNDARY_EDGE_VERTEX")  return(EG_BOUNDARY_EDGE_VERTEX);
  else return((char) - 1);
}

bool checkVector(vec3_t V)
{
  for (int i = 0; i < 3; i++) {
    if (isnan(V[i])) {
      return false;
    }
    if (isinf(V[i])) {
      return false;
    }
  }
  return true;
}

bool checkVector(vec2_t V)
{
  for (int i = 0; i < 2; i++) {
    if (isnan(V[i])) {
      return false;
    }
    if (isinf(V[i])) {
      return false;
    }
  }
  return true;
}

bool checkReal(double v)
{
  if (isnan(v)) {
    return false;
  }
  if (isinf(v)) {
    return false;
  }
  return true;
}

QDebug operator<<(QDebug dbg, const vec3_t &v)
{
  dbg.nospace() << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return dbg.space();
}

QDebug operator<<(QDebug dbg, const vec2_t &v)
{
  dbg.nospace() << "(" << v[0] << ", " << v[1] << ")";
  return dbg.space();
}

QDebug operator<<(QDebug dbg, const dcmplx &c)
{
  dbg.nospace() << real(c)<<" + "<<imag(c)<<" *i";
  return dbg.space();
}

dcmplx complex_pow(dcmplx base, double power)
{
  dcmplx i(0,1);
  double alpha = atan2(imag(base),real(base));
  return pow(abs(base),power)*exp(power*alpha*i);
}

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

int poly_solve_cubic(double a, double b, double c, double * x0, double * x1, double * x2)
{
  double q = (a * a - 3 * b);
  double r = (2 * a * a * a - 9 * a * b + 27 * c);
  
  double Q = q / 9;
  double R = r / 54;
  
  double Q3 = Q * Q * Q;
  double R2 = R * R;
  
  double CR2 = 729 * r * r;
  double CQ3 = 2916 * q * q * q;
  
  if (R == 0 && Q == 0)
  {
    *x0 = - a / 3 ;
    *x1 = - a / 3 ;
    *x2 = - a / 3 ;
    return 3 ;
  }
  else if (CR2 == CQ3) 
  {
    // this test is actually R2 == Q3, written in a form suitable
    // for exact computation with integers
    
    // Due to finite precision some double roots may be missed, and
    // considered to be a pair of complex roots z = x +/- epsilon i
    // close to the real axis.
    
    double sqrtQ = sqrt (Q);
    
    if (R > 0)
    {
      *x0 = -2 * sqrtQ  - a / 3;
      *x1 = sqrtQ - a / 3;
      *x2 = sqrtQ - a / 3;
    } else {
      *x0 = - sqrtQ  - a / 3;
      *x1 = - sqrtQ - a / 3;
      *x2 = 2 * sqrtQ - a / 3;
    }
    return 3;
  } else if (CR2 < CQ3) { // equivalent to R2 < Q3
    double sqrtQ = sqrt (Q);
    double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
    double theta = acos (R / sqrtQ3);
    double norm = -2 * sqrtQ;
    *x0 = norm * cos (theta / 3) - a / 3;
    *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
    *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
    
    // Sort *x0, *x1, *x2 into increasing order
    
    if (*x0 > *x1) {
      SWAP(*x0, *x1);
    }
    if (*x1 > *x2)
    {
      SWAP(*x1, *x2);
      if (*x0 > *x1) {
        SWAP(*x0, *x1);
      }
    }
    
    return 3;
  } else {
    double sgnR = (R >= 0 ? 1 : -1);
    double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
    double B = Q / A ;
    *x0 = A + B - a / 3;
    return 1;
  }
}

// a x^2 + b x + c = 0
int poly_solve_quadratic(double a, double b, double c, double * x0, double * x1)
{
  if (a == 0) {
    if (b == 0) {
      return(0);
    } else {
      *x0 = -c / b;
      return(1);
    }
  } else {
    double delta = pow(b, 2) - 4 * a * c;
    if (delta < 0) {
      return(0);
    } else {
      *x0 = (-b + sqrt(delta)) / (2 * a);
      *x1 = (-b - sqrt(delta)) / (2 * a);
      if (*x1 < *x0) {
        double tmp = *x0;
        *x0 = *x1;
        *x1 = tmp;
      }
    }
  }
  return(2);
}
