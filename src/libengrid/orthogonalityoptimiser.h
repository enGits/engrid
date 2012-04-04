#ifndef ORTHOGONALITYOPTIMISER_H
#define ORTHOGONALITYOPTIMISER_H

class OrthogonalityOptimiser;
class SurfaceProjection;

#include "optimisation.h"
#include "gridsmoother.h"

class OrthogonalityOptimiser : public GridSmoother, public Optimisation
{

protected: // types

  struct face_t {
    vec3_t c1, c2;
    QVector<vtkIdType> nodes;
  };

protected: // attributes

  QVector<face_t>             m_Faces;
  QVector<QList<int> >        m_Node2Face;
  QVector<bool>               m_Fixed;
  QVector<SurfaceProjection*> m_Projection;

  vtkIdType m_IdNode;

protected: // methods

  void getFace(vtkIdType id_cell, int subidx, QSet<vtkIdType> &nodes);
  int  getNumberOfFaces(vtkIdType id_cell);
  void buildNode2Face();
  void fixBoundaryNodes();
  double faceError(const face_t& face, vec3_t x);
  vec3_t newPosition(vtkIdType id_node);

  virtual double func(vec3_t x);
  virtual void operate();

public: // methods

  OrthogonalityOptimiser();

};

#endif // ORTHOGONALITYOPTIMISER_H
