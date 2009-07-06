#include "tricoord.h"

TriCoord::TriCoord(vtkUnstructuredGrid* surf_grid, vtkIdType id_tri, vec3_t x)
{
  m_SurfGrid = surf_grid;
  setPosition(x);
}

void TriCoord::setPosition(vec3_t x)
{
  /*
  bool intersects_face = GeometryTools::intersectEdgeAndTriangle(T.a, T.b, T.c, x1, x2, xi, ri);
  if (!intersects_face) {
    double kab = GeometryTools::intersection(T.a, T.b - T.a, xp, T.b - T.a);
    double kac = GeometryTools::intersection(T.a, T.c - T.a, xp, T.c - T.a);
    double kbc = GeometryTools::intersection(T.b, T.c - T.b, xp, T.c - T.b);
    double dab = (T.a + kab*(T.b-T.a) - xp).abs();
    double dac = (T.a + kac*(T.c-T.a) - xp).abs();
    double dbc = (T.b + kbc*(T.c-T.b) - xp).abs();
    bool set = false;
    if ((kab >= 0) && (kab <= 1)) {
      if (dab < d) {
        xi = T.a + kab*(T.b-T.a);
        d = dab;
        set = true;
      }
    }
    if ((kac >= 0) && (kac <= 1)) {
      if (dac < d) {
        xi = T.a + kac*(T.c-T.a);
        d = dac;
        set = true;
      }
    }
    if ((kbc >= 0) && (kbc <= 1)) {
      if (dbc < d) {
        xi = T.b + kbc*(T.c-T.b);
        d = dbc;
        set = true;
      }
    }
    double da = (T.a - xp).abs();
    double db = (T.b - xp).abs();
    double dc = (T.c - xp).abs();
    if (da < d) {
      xi = T.a;
      d = da;
      set = true;
    }
    if (db < d) {
      xi = T.b;
      d = db;
    }
    if (dc < d) {
      xi = T.c;
      d = dc;
      set = true;
    }
    if (!set) {
      EG_BUG;
    }
  }
  */
}
