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
//
#ifndef BOUNDARYLAYEROPERATION_H
#define BOUNDARYLAYEROPERATION_H

#include "surfaceoperation.h"
#include "edgelengthsourcemanager.h"
#include "cadinterface.h"

class BoundaryLayerOperation;


class BoundaryLayerOperation : public SurfaceOperation
{

protected: // data types

  enum nodetype_t { NormalNode, EdgeNode, CornerNode };


protected: // attributes

  QVector<vec3_t>           m_BoundaryLayerVectors;
  QVector<int>              m_BoundaryLayerCodes;
  QVector<bool>             m_BoundaryLayerNode;
  QVector<nodetype_t>       m_NodeTypes;
  QVector<QSet<vtkIdType> > m_SnapPoints;
  QVector<double>           m_Height;
  QSet<int>                 m_LayerAdjacentBoundaryCodes;
  double                    m_FeatureAngle;
  double                    m_StretchingRatio;
  double                    m_FarfieldRatio;
  double                    m_RadarAngle;
  double                    m_MaxHeightInGaps;
  double                    m_FaceSizeLowerLimit;
  double                    m_FaceSizeUpperLimit;
  double                    m_FaceAngleLimit;
  bool                      m_UseGrouping;
  double                    m_GroupingAngle;
  int                       m_NumBoundaryLayerVectorRelaxations;
  int                       m_NumBoundaryLayerHeightRelaxations;
  int                       m_NumShellRelaxations;
  int                       m_NumLayers;
  EdgeLengthSourceManager   m_ELSManagerBLayer;
  EdgeLengthSourceManager   m_ELSManagerSurface;


protected: // methods

  void readSettings();
  void correctBoundaryLayerVectors();
  void computeBoundaryLayerVectors();
  void addToSnapPoints(vtkIdType id_node, vtkIdType id_snap);
  void computeNodeTypes();
  void smoothBoundaryLayerVectors(int n_iter, double w_iso = 1.0, double w_dir = 0.0, QVector<bool> *node_fixed = NULL);
  void writeBoundaryLayerVectors(QString file_name, int counter = -1);
  void computeDesiredHeights();
  bool faceFine(vtkIdType id_face, double scale);
  void computeHeights();

  void   createSmoothShell(vtkUnstructuredGrid *shell_grid, int num_iter);
  void   fixBoundaryLayerVectors(const QList<vtkIdType> &bad_cells, int num_smooth_iter);
  double largestAngle(vtkIdType id_node1, vtkIdType id_node2);
  void   smoothUsingBLVectors();
  void   writeWallGrid(QString file_name, int counter = -1);

  void laplacianIntersectSmoother(const QVector<bool>& on_boundary);
  void weightedSmoother(const QVector<bool>& on_boundary);
  void angleSmoother(const QVector<bool>& on_boundary, const QVector<bool>& is_convex, QVector<vec3_t>& grid_pnts);
  void intersectSmoother(const QVector<bool>& on_boundary, const QVector<bool>& is_convex, QVector<vec3_t>& grid_pnts);
  void laplacianSmoother();
  void pushOut(const QVector<bool>& on_boundary, const QVector<bool>& is_convex);
  int  limitHeights(double safety_factor);
  void snapCornerVectorsToShell(vtkUnstructuredGrid* shell_grid);
  void newHeightFromShellIntersect(vtkUnstructuredGrid* shell_grid);
  void limitSizeAndAngleErrors();
  bool swapRequired(stencil_t stencil, CadInterface *cad, double threshold_angle);
  void swapEdgesToMatchShell(vtkUnstructuredGrid *shell_grid, double threshold_angle);


public: // methods

  QSet<int> getLayerAdjacentBoundaryCodes();



};

#endif // BOUNDARYLAYEROPERATION_H
