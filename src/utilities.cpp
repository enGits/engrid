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

int cout_grid(ostream &stream, vtkUnstructuredGrid *grid, bool npoints, bool ncells, bool points, bool cells)
{
  stream<<"============="<<endl;
  if(npoints) stream << "grid->GetNumberOfPoints()=" << grid->GetNumberOfPoints() << endl;
  if(ncells) stream << "grid->GetNumberOfCells()=" << grid->GetNumberOfCells() << endl;
  if(points) {
    for (vtkIdType i = 0; i < grid->GetNumberOfPoints(); ++i) {
      vec3_t x;
      grid->GetPoint(i, x.data());
      stream << "Vertex " << i << " = " << x << endl;
    }
  }
  if(cells) {
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
      vtkCell *C = (vtkCell *) vtkCell::New();
      C=grid->GetCell(i);
      vtkIdType npts=C->GetNumberOfPoints();
      vtkIdType* pts;
      grid->GetCellPoints(i, npts, pts);
      stream << "Cell " << i << " = ";
      for(int j=0;j<npts;j++) stream << pts[j] << " ";
      stream << "boundary_code=" << cell_code->GetValue(i);
      stream << endl;
    }
  }
  stream<<"============="<<endl;
  return 0;
}

///////////////////////////////////////////
//Warning: Untested
int addCell(vtkUnstructuredGrid* a_grid, vtkIdType A, vtkIdType B, vtkIdType C, int bc)
{
  vtkIdType npts=3;
  vtkIdType pts[3];
  pts[0]=A;
  pts[1]=B;
  pts[2]=C;
  vtkIdType newCellId = a_grid->InsertNextCell(VTK_TRIANGLE,npts,pts);
  EG_VTKDCC(vtkIntArray, cell_code, a_grid, "cell_code");
  cell_code->SetValue(newCellId, bc);
  return(0);
}

///////////////////////////////////////////

int getShortestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid)
{
  vtkIdType N_pts, *pts;
  a_grid->GetCellPoints(a_id_cell, N_pts, pts);
  vec3_t* x=new vec3_t[N_pts];
  for(int i=0;i<N_pts;i++) a_grid->GetPoints()->GetPoint(pts[i], x[i].data());
  int id_minlen=0;
  double minlen=(x[1]-x[0]).abs();
  for(int i=1;i<N_pts;i++)
  {
    double len=(x[(i+1)%N_pts]-x[i]).abs();
    if(len<minlen){
      minlen=len;
      id_minlen=i;
    }
  }
  delete x;
  return(id_minlen);
}

int getLongestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid)
{
  vtkIdType N_pts, *pts;
  a_grid->GetCellPoints(a_id_cell, N_pts, pts);
  vec3_t* x=new vec3_t[N_pts];
  for(int i=0;i<N_pts;i++) a_grid->GetPoints()->GetPoint(pts[i], x[i].data());
  int id_maxlen=0;
  double maxlen=(x[1]-x[0]).abs();
  cout<<"maxlen="<<maxlen<<endl;
  for(int i=1;i<N_pts;i++)
  {
    double len=(x[(i+1)%N_pts]-x[i]).abs();
    cout<<"len["<<i<<"]="<<len<<endl;
    if(len>maxlen){
      maxlen=len;
      id_maxlen=i;
    }
  }
  delete x;
  return(id_maxlen);
}

int getSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid,vtkIdType a_id_node1,vtkIdType a_id_node2)
{
  vtkIdType N_pts, *pts;
  a_grid->GetCellPoints(a_id_cell, N_pts, pts);
  QVector <vtkIdType> edge(2);
  
  int n=0;
  for(int i=0;i<N_pts;i++)
  {
    if(pts[i]==a_id_node1) { edge[0]=i;n++;}
    if(pts[i]==a_id_node2) { edge[1]=i;n++;}
  }
  if(n!=2){
    EG_BUG;
    return(-1);
  }
  qSort(edge.begin(),edge.end());
  if(edge[0]==0 && edge[1]==N_pts-1) return(N_pts-1);
  else return(edge[0]);
}
///////////////////////////////////////////

QString cell2str(vtkIdType id_cell,vtkUnstructuredGrid* grid)
{
  QString tmp;
  tmp.setNum(id_cell);
  QString txt = "id_cell=" + tmp;
  
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(id_cell, N_pts, pts);
  
  txt += " [";
  for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
    tmp.setNum(pts[i_pts]);
    txt += tmp;
    if (i_pts < N_pts-1) {
      txt += ",";
    }
  }
  txt += "]";
  return(txt);
}



///////////////////////////////////////////

pair<vtkIdType,vtkIdType> OrderedPair(vtkIdType a, vtkIdType b)
{
  vtkIdType x=min(a,b);
  vtkIdType y=max(a,b);
  return(pair<vtkIdType,vtkIdType>(x,y));
}

const char* VertexType2Str(char T)
{
  if(T==VTK_SIMPLE_VERTEX) return("VTK_SIMPLE_VERTEX");
  if(T==VTK_FIXED_VERTEX) return("VTK_FIXED_VERTEX");
  if(T==VTK_FEATURE_EDGE_VERTEX) return("VTK_FEATURE_EDGE_VERTEX");
  if(T==VTK_BOUNDARY_EDGE_VERTEX) return("VTK_BOUNDARY_EDGE_VERTEX");
  else return("Unknown vertex type");
}

char Str2VertexType(QString S)
{
  if(S=="VTK_SIMPLE_VERTEX") return(VTK_SIMPLE_VERTEX);
  if(S=="VTK_FIXED_VERTEX") return(VTK_FIXED_VERTEX);
  if(S=="VTK_FEATURE_EDGE_VERTEX") return(VTK_FEATURE_EDGE_VERTEX);
  if(S=="VTK_BOUNDARY_EDGE_VERTEX") return(VTK_BOUNDARY_EDGE_VERTEX);
  else return((char)-1);
}

bool checkVector(vec3_t V)
{
  for(int i=0;i<3;i++) {
    if(isnan(V[i])) {
      EG_ERR_RETURN("NAN");
      return false;
    }
    if(isinf(V[i])) {
      EG_ERR_RETURN("INFINITY");
      return false;
    }
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
