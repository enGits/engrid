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
#ifndef CADINTERFACE_H
#define CADINTERFACE_H

#include "engrid.h"
#include "meshpartition.h"
#include "surfacealgorithm.h"

#include <QVector>
#include <QPair>

class CadInterface : public SurfaceAlgorithm
{

private: // attributes

  QString m_Name; ///< The name of the CAD interface


protected: // attributes

  vtkUnstructuredGrid* m_FGrid;               ///< the foreground grid to project
  MeshPartition        m_FPart;               ///< MeshPartition for the foreground grid
  QVector<vec3_t>      m_RayPoints;           ///< template for ray cloud to snap points
  vec3_t               m_LastNormal;          ///< last surface normal which was encountered
  double               m_LastRadius;          ///< last surface radius which was encountered
  bool                 m_Failed;              ///< did the last operation fail
  bool                 m_ShootRayImplemented; ///< flag to determine if shootRay has been implemented
  double               m_CriticalSnapLength;  ///< relative distance to decide between project or real snap


protected: // methods

  /**
   * @brief issue an error message about not implemented functionality
   * Not all functionality will be implemented in every derived class of CadInterface.
   * Certain meshing processes might still work. For this reason an error message will be displayed,
   * rather than preventing instantiation by using abstract methods (=0).
   */
  void notImplemented();

  /**
   * @brief set the name
   * The CAD interface should be given a name. This will be used to display messages.
   * A name could, for example, be "triangulated geometry CAD interface".
   * @param name the name to set
   */
  void setName(QString name) { m_Name = name; }


public: // data types

  enum HitType { Miss, HitIn, HitOut };


  //enum PositionType { Inside, Outside, Surface };


public:

  CadInterface();

  void setForegroundGrid(vtkUnstructuredGrid* grid);

  QString name() { return m_Name; }

  /**
   * @brief vital interface method to the geometry
   * This method is a major part of interfacing geometries. The idea is to provide something like
   * a virtual laser scanner. The concept of a ray-tracing interface has been taken from the
   * open-source CAD software BRL-CAD (see http://brlcad.org).
   * @param x the origin of the ray
   * @param v the direction of the ray
   * @param x_hit if the ray hits this contains the intersection point between ray and geometry
   * @param n_hit if the ray hits this contains the geometry normal vector at the intersection point
   * @param r if the ray hits this contains the surface radius at the intersection point
   * @return the result (Miss, HitIn, HitOut)
   */
  virtual HitType shootRay(vec3_t x, vec3_t v, vec3_t &x_hit, vec3_t &n_hit, double &r);

  /**
   * @brief compute all intersections of a ray and the CAD geometry.
   * @param x the origin of the ray
   * @param v the direction of the ray
   * @param intersections will hold all intersections with the goemtry
   */
  virtual void computeIntersections(vec3_t x, vec3_t v, QVector<QPair<vec3_t,vtkIdType> > &intersections) { notImplemented(); }

  /**
   * @brief check if shootRay is available (implemented)
   * @return true id shootRay is available
   */
  bool shootRayAvailable() { return m_ShootRayImplemented; }

  /**
   * @brief snap a point to the geometry
   * If you want to snap a node of the grid, consider using snapNode because it might be faster.
   * Curvature correction only applies to certain geometric models (e.g. STL).
   * @param x the position of the point to snap
   * @param correct_curvature flag to determine if corveture correction shall be used
   * @return the snapped position
   */
  virtual vec3_t snap(vec3_t x, bool correct_curvature = false);

  virtual vec3_t snapToEdge(vec3_t x) { notImplemented(); }
  virtual vec3_t snapToCorner(vec3_t x) { notImplemented(); }

  /**
   * @brief snap a node of the foreground grid
   * Curvature correction only applies to certain geometric models (e.g. STL).
   * @param id_node the node to snap
   * @param correct_curvature flag to determine if corveture correction shall be used
   * @return the snapped position
   */
  virtual vec3_t  snapNode(vtkIdType id_node, bool correct_curvature = false);

  /**
   * @brief snap a node of the foreground grid using an alternate position
   * Curvature correction only applies to certain geometric models (e.g. STL).
   * @param id_node the node to snap
   * @param x the position to snap
   * @param correct_curvature flag to determine if curvature correction shall be used
   * @return the snapped position
   */
  virtual vec3_t  snapNode(vtkIdType id_node, vec3_t x, bool correct_curvature = false);

  /**
   * @brief project a point onto the geometry
   * If you want to project a node of the mesh, consider using projectNode because it might be faster.
   * This methods project in forward and backward direction and takes the closer point.
   * Curvature correction only applies to certain geometric models (e.g. STL).
   * @param x the position of the point to project
   * @param v the direction of projection
   * @param strict_direction if set to true only the forward direction will be considered
   * @param correct_curvature set this to true if you want to use curvature correction
   * @return the projected point
   */
  virtual vec3_t project(vec3_t x, vec3_t v, bool strict_direction = false, bool correct_curvature = false);

  /**
   * @brief project a node of the foreground mesh onto the geometry
   * This methods project in forward and backward direction and takes the closer point.
   * Curvature correction only applies to certain geometric models (e.g. STL).
   * @param id_node the node to project
   * @param strict_direction if set to true only the forward direction will be considered
   * @param correct_curvature set this to true if you want to use curvature correction
   * @return the projected point
   */
  virtual vec3_t projectNode(vtkIdType id_node, bool strict_direction = false, bool correct_curvature = false);

  /**
   * @brief project a node of the foreground mesh onto the geometry using an alternate position
   * This methods project in forward and backward direction and takes the closer point.
   * Curvature correction only applies to certain geometric models (e.g. STL).
   * @param id_node the node to project
   * @param x the position to project
   * @param strict_direction if set to true only the forward direction will be considered
   * @param correct_curvature set this to true if you want to use curvature correction
   * @return the projected point
   */
  virtual vec3_t projectNode(vtkIdType id_node, vec3_t x, bool strict_direction = false, bool correct_curvature = false);

  /**
   * @brief project a node of the foreground mesh onto the geometry using an alternate position and direction
   * This methods project in forward and backward direction and takes the closer point.
   * Curvature correction only applies to certain geometric models (e.g. STL).
   * @param id_node the node to project
   * @param x the position to project
   * @param v the direction of projection
   * @param strict_direction if set to true only the forward direction will be considered
   * @param correct_curvature set this to true if you want to use curvature correction
   * @return the projected point
   */
  virtual vec3_t projectNode(vtkIdType id_node, vec3_t x, vec3_t v, bool strict_direction = false, bool correct_curvature = false);

  virtual vec3_t correctCurvature(vec3_t x) { return x; }

  virtual double getRadius(vtkIdType id_node);

  //vtkIdType lastProjTriangle() { return -1; } /// delete this ???

  bool      failed()           { return m_Failed; }
  vec3_t    getLastNormal()    { return m_LastNormal; }
  double    getLastRadius()    { return m_LastRadius; }

};

#endif // CADINTERFACE_H
