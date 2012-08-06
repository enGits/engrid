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

#include "brlcadinterface.h"

vec3_t BrlCadInterface::m_XIn;
vec3_t BrlCadInterface::m_XOut;
vec3_t BrlCadInterface::m_InNormal;
vec3_t BrlCadInterface::m_OutNormal;
bool   BrlCadInterface::m_Hit;
double BrlCadInterface::m_InRadius;
double BrlCadInterface::m_OutRadius;

BrlCadInterface::BrlCadInterface()
{
  m_Epsilon = 1e-10;
}

void BrlCadInterface::setupBrlCad(QString file_name, QString object_name)
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

int BrlCadInterface::hit(application *ap, struct partition *PartHeadp, seg *segs)
{
  register struct partition *pp;
  register struct hit *hitp;
  register struct soltab *stp;
  struct curvature cur;
  point_t		pt;
  vect_t		inormal;
  vect_t		onormal;
  double curv;
  int N = 0;
  for (pp=PartHeadp->pt_forw; pp != PartHeadp; pp = pp->pt_forw) {
    ++N;

    hitp = pp->pt_inhit;
    stp  = pp->pt_inseg->seg_stp;

    VJOIN1(pt, ap->a_ray.r_pt, hitp->hit_dist, ap->a_ray.r_dir);

    RT_HIT_NORMAL(inormal, hitp, stp, &(ap->a_ray), pp->pt_inflip);
    RT_CURVATURE(&cur, hitp, pp->pt_inflip, stp);
    curv = max(fabs(cur.crv_c1), fabs(cur.crv_c2));
    m_InRadius = 1.0/max(1e-10, curv);

    m_XIn[0] = pt[0];
    m_XIn[1] = pt[1];
    m_XIn[2] = pt[2];
    m_InNormal[0] = inormal[0];
    m_InNormal[1] = inormal[1];
    m_InNormal[2] = inormal[2];

    hitp = pp->pt_outhit;
    stp  = pp->pt_outseg->seg_stp;
    VJOIN1(pt, ap->a_ray.r_pt, hitp->hit_dist, ap->a_ray.r_dir);
    RT_HIT_NORMAL( onormal, hitp, stp, &(ap->a_ray), pp->pt_outflip );
    RT_CURVATURE(&cur, hitp, pp->pt_inflip, stp);
    curv = max(fabs(cur.crv_c1), fabs(cur.crv_c2));
    m_OutRadius = 1.0/max(1e-10, curv);

    m_XOut[0] = pt[0];
    m_XOut[1] = pt[1];
    m_XOut[2] = pt[2];
    m_OutNormal[0] = onormal[0];
    m_OutNormal[1] = onormal[1];
    m_OutNormal[2] = onormal[2];
  }
  m_Hit = true;
}

int BrlCadInterface::miss(application *ap)
{
  m_Hit = false;
}

bool BrlCadInterface::shootOneRay(vec3_t x, vec3_t v, vec3_t &x_in, vec3_t &x_out, vec3_t &n_in, vec3_t &n_out, double &r_in, double &r_out)
{
  if (!checkVector(x)) {
    EG_BUG;
  }
  if (!checkVector(v)) {
    EG_BUG;
  }
  VSET(m_Ap.a_ray.r_pt,  x[0], x[1], x[2]);
  VSET(m_Ap.a_ray.r_dir, v[0], v[1], v[2]);
  m_Hit = false;
  m_Ap.a_hit  = BrlCadInterface::hit;
  m_Ap.a_miss = BrlCadInterface::miss;
  rt_shootray(&m_Ap);
  x_in = m_XIn;
  n_in = m_InNormal;
  r_in = m_InRadius;
  x_out = m_XOut;
  n_out = m_OutNormal;
  r_out = m_OutRadius;
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

BrlCadInterface::HitType BrlCadInterface::shootRay(vec3_t x, vec3_t v, vec3_t &x_hit, vec3_t &n_hit, double &r)
{
  HitType hit_type = Miss;
  v.normalise();
  vec3_t x_in, x_out, n_in, n_out;
  double r_in, r_out;
  if (shootOneRay(x, v, x_in, x_out, n_in, n_out, r_in, r_out)) {
    double d_in = (x_in - x)*v;
    if (d_in > 0) {
      x_hit = x_in;
      n_hit = n_in;
      r = r_in;
      hit_type = HitIn;
    }
    double d_out = (x_out - x)*v;
    if (d_out > 0) {
      if (hit_type == Miss || d_out < d_in) {
        x_hit = x_out;
        n_hit = n_out;
        hit_type = HitOut;
        r = r_out;
      }
    }
  }
  return hit_type;
}

BrlCadInterface::PositionType BrlCadInterface::position(vec3_t x, vec3_t n)
{
  vec3_t x_hit, n_hit;
  double r_hit;
  HitType hit_type = shootRay(x, vec3_t(1,0,0), x_hit, n_hit, r_hit);
  if (hit_type == HitOut) {
    return Inside;
  }
  if (hit_type == HitIn) {
    return Outside;
  }
  return Surface;
}
