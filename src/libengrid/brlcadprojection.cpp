#include "brlcadprojection.h"


// Load database.
// rt_dirbuild() returns an "instance" pointer which describes
// the database to be ray traced.  It also gives you back the
// title string in the header (ID) record.
void BrlCadProjection::load(QString file_name)
{
  if ((m_Rtip = rt_dirbuild(qPrintable(file_name), m_IdBuf, sizeof(m_IdBuf))) == RTI_NULL ) {
    EG_ERR_RETURN("cannot open BRL-CAD database");
  }
  m_App.a_rt_i = m_Rtip;	// your application uses this instance
  cout << "db title: " << m_IdBuf << endl;
}

// Walk trees.
// Here you identify any object trees in the database that you
// want included in the ray trace.
void BrlCadProjection::prepare(QString object_name)
{
  if (rt_gettree(m_Rtip, qPrintable(object_name)) < 0) {
    EG_ERR_RETURN("cannot prepare BRL-CAD database");
  }

  // This next call gets the database ready for ray tracing.
  // (it precomputes some values, sets up space partitioning, etc.)
  rt_prep_parallel(m_Rtip, 1);
}

/*
// Set the ray start point and direction
// rt_shootray() uses these two to determine what ray to fire.
// In this case we simply shoot down the z axis toward the
// origin from 10 meters away [librt assumes units of millimeters.
// not that is really maters here, but an MGED database made with
// units=mm will have the same values in the file (and thus in
// librt) that you see displayed by MGED.
VSET( ap.a_ray.r_pt, 0, 0, 10000 );
VSET( ap.a_ray.r_dir, 0, 0, -1 );

VPRINT( "Pnt", ap.a_ray.r_pt );
VPRINT( "Dir", ap.a_ray.r_dir );

//
ap.a_hit = hit;			// where to go on a hit
ap.a_miss = miss;		// where to go on a miss
(void)rt_shootray( &ap );	// do it


// A real application would probably set up another
// ray and fire again.

return(0);
}


// rt_shootray() was told to call this on a hit.  He gives up the
// application structure which describes the state of the world
// (see raytrace.h), and a circular linked list of partitions,
// each one describing one in and out segment of one region.
hit(register struct application *ap, struct partition *PartHeadp, struct seg *segs)
{
// see raytrace.h for all of these guys
register struct partition *pp;
register struct hit *hitp;
register struct soltab *stp;
struct curvature cur;
point_t		pt;
vect_t		inormal;
vect_t		onormal;

// examine each partition until we get back to the head
for( pp=PartHeadp->pt_forw; pp != PartHeadp; pp = pp->pt_forw )  {
    bu_log("\n--- Hit region %s (in %s, out %s)\n",
           pp->pt_regionp->reg_name,
           pp->pt_inseg->seg_stp->st_name,
           pp->pt_outseg->seg_stp->st_name );

    // inhit info
    hitp = pp->pt_inhit;
    stp = pp->pt_inseg->seg_stp;

    VJOIN1( pt, ap->a_ray.r_pt, hitp->hit_dist, ap->a_ray.r_dir );

    // This macro takes care of the flip flag and all that
    RT_HIT_NORMAL( inormal, hitp, stp, &(ap->a_ray), pp->pt_inflip );

    rt_pr_hit( "  In", hitp );
    VPRINT(    "  Ipoint", pt );
    VPRINT(    "  Inormal", inormal );

    // This next macro fills in the curvature information
    // which consists on a principle direction vector, and
    // the inverse radii of curvature along that direction
    // and perpendicular to it.  Positive curvature bends
    // toward the outward pointing normal.
    RT_CURVATURE( &cur, hitp, pp->pt_inflip, stp );
    VPRINT("PDir", cur.crv_pdir );
    bu_log(" c1=%g\n", cur.crv_c1);
    bu_log(" c2=%g\n", cur.crv_c2);

    // outhit info
    hitp = pp->pt_outhit;
    stp = pp->pt_outseg->seg_stp;
    VJOIN1( pt, ap->a_ray.r_pt, hitp->hit_dist, ap->a_ray.r_dir );
    RT_HIT_NORMAL( onormal, hitp, stp, &(ap->a_ray), pp->pt_outflip );

    rt_pr_hit( "  Out", hitp );
    VPRINT(    "  Opoint", pt );
    VPRINT(    "  Onormal", onormal );
}

// A more complicated application would probably fill in a
// new local application structure and describe say a reflected
// or refracted ray, and then call rt_shootray with it.

// This value is returned by rt_shootray
// a hit usually returns 1, miss 0.

return(1);
}

// rt_shootray() was told to call this on a miss.
miss(register struct application *ap)
{
bu_log("missed\n");
return(0);
}
*/

BrlCadProjection::BrlCadProjection(QString file_name)
{

}

BrlCadProjection::~BrlCadProjection()
{

}

vec3_t BrlCadProjection::projectFree(vec3_t x, vtkIdType id_node)
{

}

vec3_t BrlCadProjection::projectRestricted(vec3_t x, vtkIdType id_node)
{

}
