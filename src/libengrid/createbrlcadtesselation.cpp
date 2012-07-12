//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#include "createbrlcadtesselation.h"

vec3_t CreateBrlCadTesselation::m_XIn;
vec3_t CreateBrlCadTesselation::m_XOut;
vec3_t CreateBrlCadTesselation::m_InNormal;
vec3_t CreateBrlCadTesselation::m_OutNormal;
bool   CreateBrlCadTesselation::m_Hit;

CreateBrlCadTesselation::CreateBrlCadTesselation(QString file_name, QString object_name)
{
  m_Rtip = rt_dirbuild(qPrintable(file_name), m_IdBuf, sizeof(m_IdBuf));
  if (m_Rtip == RTI_NULL) {
    EG_ERR_RETURN("Unable to open BRL-CAD database!");
  }
  if (rt_gettree(m_Rtip, qPrintable(object_name)) < 0) {
    EG_ERR_RETURN("unable to access selected object");
  }

  application ap = {0};
  m_Ap = ap;
  m_Ap.a_rt_i   = m_Rtip;
  m_Ap.a_onehit = 1; // X-ray functionality

  rt_prep_parallel(m_Rtip, 1);
}

int CreateBrlCadTesselation::hit(application *ap, struct partition *PartHeadp, seg *segs)
{
  register struct partition *pp;
  register struct hit *hitp;
  register struct soltab *stp;
  struct curvature cur;
  point_t		pt;
  vect_t		inormal;
  vect_t		onormal;
  for( pp=PartHeadp->pt_forw; pp != PartHeadp; pp = pp->pt_forw )  {

    // inhit info
    hitp = pp->pt_inhit;
    stp  = pp->pt_inseg->seg_stp;

    VJOIN1(pt, ap->a_ray.r_pt, hitp->hit_dist, ap->a_ray.r_dir);

    // This macro takes care of the flip flag and all that
    RT_HIT_NORMAL(inormal, hitp, stp, &(ap->a_ray), pp->pt_inflip);

    m_XIn[0] = pt[0];
    m_XIn[1] = pt[1];
    m_XIn[2] = pt[2];
    m_InNormal[0] = inormal[0];
    m_InNormal[1] = inormal[1];
    m_InNormal[2] = inormal[2];

    /* This next macro fills in the curvature information
     * which consists on a principle direction vector, and
     * the inverse radii of curvature along that direction
     * and perpendicular to it.  Positive curvature bends
     * toward the outward pointing normal. */
    /*
    RT_CURVATURE( &cur, hitp, pp->pt_inflip, stp );
    VPRINT("PDir", cur.crv_pdir );
    bu_log(" c1=%g\n", cur.crv_c1);
    bu_log(" c2=%g\n", cur.crv_c2);
    */

    /* outhit info */
    hitp = pp->pt_outhit;
    stp  = pp->pt_outseg->seg_stp;
    VJOIN1(pt, ap->a_ray.r_pt, hitp->hit_dist, ap->a_ray.r_dir);
    RT_HIT_NORMAL( onormal, hitp, stp, &(ap->a_ray), pp->pt_outflip );

    //rt_pr_hit( "  Out", hitp );
    //VPRINT(    "  Opoint", pt );
    //VPRINT(    "  Onormal", onormal );
    m_XOut[0] = pt[0];
    m_XOut[1] = pt[1];
    m_XOut[2] = pt[2];
    m_OutNormal[0] = onormal[0];
    m_OutNormal[1] = onormal[1];
    m_OutNormal[2] = onormal[2];
  }
  m_Hit = true;
}

int CreateBrlCadTesselation::miss(application *ap)
{
  m_Hit = false;
}

bool CreateBrlCadTesselation::shootOneRay(vec3_t x, vec3_t v, vec3_t &x_in, vec3_t &x_out, vec3_t &n_in, vec3_t &n_out)
{
  if (!checkVector(x)) {
    EG_BUG;
  }
  if (!checkVector(v)) {
    EG_BUG;
  }
  bool first = true;
  VSET(m_Ap.a_ray.r_pt,  x[0], x[1], x[2]);
  VSET(m_Ap.a_ray.r_dir, v[0], v[1], v[2]);
  m_Hit = false;
  m_Ap.a_hit  = CreateBrlCadTesselation::hit;
  m_Ap.a_miss = CreateBrlCadTesselation::miss;
  rt_shootray(&m_Ap);
  x_in = m_XIn;
  n_in = m_InNormal;
  x_out = m_XOut;
  n_out = m_OutNormal;
  x = x_out;
  if (!checkVector(x_in)) {
    EG_BUG;
  }
  if (!checkVector(x_out)) {
    EG_BUG;
  }
  if (!checkVector(n_in)) {
    EG_BUG;
  }
  if (!checkVector(n_out)) {
    EG_BUG;
  }
  return m_Hit;
}

bool CreateBrlCadTesselation::shootRay(vec3_t x, vec3_t v, vec3_t &x_in, vec3_t &x_out, vec3_t &n_in, vec3_t &n_out)
{
  v.normalise();
  double epsilon = min(1e-2*m_SmallestFeatureSize, 1e-5*(m_X1 - m_X2).abs());
  vec3_t x1, x2, n1, n2;
  if (shootOneRay(x, v, x1, x2, n1, n2)) {
    x_in = x1;
    n_in = n1;
    x_out = x2;
    n_out = n2;
    bool finished = false;
    while (!finished) {
      if (!shootOneRay(x2 + epsilon*v, v, x1, x2, n1, n2)) {
        finished = true;
      } else {
        if ((x1-x_out).abs() > epsilon) {
          finished = true;
        } else {
          x_out = x2;
          n_out = n2;
        }
      }
    }
    return true;
  } else {
    return false;
  }
}
