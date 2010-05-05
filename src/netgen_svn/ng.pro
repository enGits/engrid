QT       -= gui

win32 {
	TARGET = nglib
}

TEMPLATE = lib
CONFIG += staticlib

LANGUAGE     = C++
CONFIG      += release warn_off
INCLUDEPATH += netgen-mesher/netgen/libsrc/include
INCLUDEPATH += ./netgen-mesher/netgen/nglib

#INCLUDEPATH += .

win32 {
	DEFINES     += NO_PARALLEL_THREADS
}

win32 {
	DEFINES     += NGLIB_EXPORTS
}

win32 {
	DEFINES     += DLL_EXPORT
}

DEFINES     += vtkRendering_EXPORTS
#DEFINES     += PIC
#DEFINES     += HAVE_CONFIG_H

#PreprocessorDefinitions="WIN32;_DEBUG;_WINDOWS;_USRDLL;NGLIB_EXPORTS;MSVC_EXPRESS;WINVER=0x0600;NTDDI_VERSION=NTDDI_VISTA"
#PreprocessorDefinitions="WIN32;_DEBUG;_WINDOWS;_USRDLL;NGLIB_EXPORTS;MSVC_EXPRESS;WINVER=0x0600;NTDDI_VERSION=NTDDI_VISTA"
#PreprocessorDefinitions="WIN32;NDEBUG;_WINDOWS;_USRDLL;NGLIB_EXPORTS;MSVC_EXPRESS;WINVER=0x0600;NTDDI_VERSION=NTDDI_VISTA"
#PreprocessorDefinitions="WIN32;NDEBUG;_WINDOWS;_USRDLL;NGLIB_EXPORTS;MSVC_EXPRESS;WINVER=0x0600;NTDDI_VERSION=NTDDI_VISTA"
#PreprocessorDefinitions="WNT;WIN32;NDEBUG;_WINDOWS;_USRDLL;NGLIB_EXPORTS;MSVC_EXPRESS;OCCGEOMETRY;WINVER=0x0600;NTDDI_VERSION=NTDDI_VISTA"
#PreprocessorDefinitions="WNT;WIN32;NDEBUG;_WINDOWS;_USRDLL;NGLIB_EXPORTS;MSVC_EXPRESS;OCCGEOMETRY;_OCC64;WINVER=0x0600;NTDDI_VERSION=NTDDI_VISTA"

SOURCES  = \
./netgen-mesher/netgen/libsrc/occ/occmeshsurf.cpp \
./netgen-mesher/netgen/libsrc/occ/occgenmesh.cpp \
./netgen-mesher/netgen/libsrc/occ/occconstruction.cpp \
./netgen-mesher/netgen/libsrc/occ/occgeom.cpp \
./netgen-mesher/netgen/libsrc/interface/writefeap.cpp \
./netgen-mesher/netgen/libsrc/interface/writegmsh.cpp \
./netgen-mesher/netgen/libsrc/interface/writetet.cpp \
./netgen-mesher/netgen/libsrc/interface/writejcm.cpp \
./netgen-mesher/netgen/libsrc/interface/writedolfin.cpp \
./netgen-mesher/netgen/libsrc/interface/writetochnog.cpp \
./netgen-mesher/netgen/libsrc/interface/read_fnf_mesh.cpp \
./netgen-mesher/netgen/libsrc/interface/writefluent.cpp \
./netgen-mesher/netgen/libsrc/interface/readtetmesh.cpp \
./netgen-mesher/netgen/libsrc/interface/writetecplot.cpp \
./netgen-mesher/netgen/libsrc/interface/writeabaqus.cpp \
./netgen-mesher/netgen/libsrc/interface/writeuser.cpp \
./netgen-mesher/netgen/libsrc/interface/readuser.cpp \
./netgen-mesher/netgen/libsrc/interface/writediffpack.cpp \
./netgen-mesher/netgen/libsrc/interface/writepermas.cpp \
./netgen-mesher/netgen/libsrc/interface/writeelmer.cpp \
./netgen-mesher/netgen/libsrc/interface/writegmsh2.cpp \
./netgen-mesher/netgen/libsrc/interface/wuchemnitz.cpp \
./netgen-mesher/netgen/libsrc/general/ngexception.cpp \
./netgen-mesher/netgen/libsrc/general/parthreads.cpp \
./netgen-mesher/netgen/libsrc/general/profiler.cpp \
./netgen-mesher/netgen/libsrc/general/array.cpp \
./netgen-mesher/netgen/libsrc/general/sort.cpp \
./netgen-mesher/netgen/libsrc/general/optmem.cpp \
./netgen-mesher/netgen/libsrc/general/hashtabl.cpp \
./netgen-mesher/netgen/libsrc/general/spbita2d.cpp \
./netgen-mesher/netgen/libsrc/general/symbolta.cpp \
./netgen-mesher/netgen/libsrc/general/flags.cpp \
./netgen-mesher/netgen/libsrc/general/bitarray.cpp \
./netgen-mesher/netgen/libsrc/general/seti.cpp \
./netgen-mesher/netgen/libsrc/general/dynamicmem.cpp \
./netgen-mesher/netgen/libsrc/general/mystring.cpp \
./netgen-mesher/netgen/libsrc/general/table.cpp \
./netgen-mesher/netgen/libsrc/stlgeom/stlgeomchart.cpp \
./netgen-mesher/netgen/libsrc/stlgeom/stltool.cpp \
./netgen-mesher/netgen/libsrc/stlgeom/stlgeommesh.cpp \
./netgen-mesher/netgen/libsrc/stlgeom/stltopology.cpp \
./netgen-mesher/netgen/libsrc/stlgeom/meshstlsurface.cpp \
./netgen-mesher/netgen/libsrc/stlgeom/stlgeom.cpp \
./netgen-mesher/netgen/libsrc/stlgeom/stlline.cpp \
./netgen-mesher/netgen/libsrc/gprim/geom2d.cpp \
./netgen-mesher/netgen/libsrc/gprim/transform3d.cpp \
./netgen-mesher/netgen/libsrc/gprim/geomfuncs.cpp \
./netgen-mesher/netgen/libsrc/gprim/adtree.cpp \
./netgen-mesher/netgen/libsrc/gprim/geom3d.cpp \
./netgen-mesher/netgen/libsrc/gprim/geomtest3d.cpp \
./netgen-mesher/netgen/libsrc/csg/curve2d.cpp \
./netgen-mesher/netgen/libsrc/csg/singularref.cpp \
./netgen-mesher/netgen/libsrc/csg/meshsurf.cpp \
./netgen-mesher/netgen/libsrc/csg/explicitcurve2d.cpp \
./netgen-mesher/netgen/libsrc/csg/edgeflw.cpp \
./netgen-mesher/netgen/libsrc/csg/manifold.cpp \
./netgen-mesher/netgen/libsrc/csg/spline3d.cpp \
./netgen-mesher/netgen/libsrc/csg/extrusion.cpp \
./netgen-mesher/netgen/libsrc/csg/surface.cpp \
./netgen-mesher/netgen/libsrc/csg/revolution.cpp \
./netgen-mesher/netgen/libsrc/csg/algprim.cpp \
./netgen-mesher/netgen/libsrc/csg/identify.cpp \
./netgen-mesher/netgen/libsrc/csg/bspline2d.cpp \
./netgen-mesher/netgen/libsrc/csg/triapprox.cpp \
./netgen-mesher/netgen/libsrc/csg/genmesh.cpp \
./netgen-mesher/netgen/libsrc/csg/gencyl.cpp \
./netgen-mesher/netgen/libsrc/csg/csgparser.cpp \
./netgen-mesher/netgen/libsrc/csg/specpoin.cpp \
./netgen-mesher/netgen/libsrc/csg/csgeom.cpp \
./netgen-mesher/netgen/libsrc/csg/solid.cpp \
./netgen-mesher/netgen/libsrc/csg/brick.cpp \
./netgen-mesher/netgen/libsrc/csg/polyhedra.cpp \
./netgen-mesher/netgen/libsrc/linalg/densemat.cpp \
./netgen-mesher/netgen/libsrc/linalg/bfgs.cpp \
./netgen-mesher/netgen/libsrc/linalg/linopt.cpp \
./netgen-mesher/netgen/libsrc/linalg/linsearch.cpp \
./netgen-mesher/netgen/libsrc/linalg/polynomial.cpp \
./netgen-mesher/netgen/libsrc/meshing/validate.cpp \
./netgen-mesher/netgen/libsrc/meshing/adfront3.cpp \
./netgen-mesher/netgen/libsrc/meshing/localh.cpp \
./netgen-mesher/netgen/libsrc/meshing/bisect.cpp \
./netgen-mesher/netgen/libsrc/meshing/basegeom.cpp \
./netgen-mesher/netgen/libsrc/meshing/hprefinement.cpp \
./netgen-mesher/netgen/libsrc/meshing/parallelmesh.cpp \
./netgen-mesher/netgen/libsrc/meshing/improve3.cpp \
./netgen-mesher/netgen/libsrc/meshing/smoothing2.5.cpp \
./netgen-mesher/netgen/libsrc/meshing/msghandler.cpp \
./netgen-mesher/netgen/libsrc/meshing/triarls.cpp \
./netgen-mesher/netgen/libsrc/meshing/smoothing2.cpp \
./netgen-mesher/netgen/libsrc/meshing/meshfunc2d.cpp \
./netgen-mesher/netgen/libsrc/meshing/pyramid2rls.cpp \
./netgen-mesher/netgen/libsrc/meshing/prism2rls.cpp \
./netgen-mesher/netgen/libsrc/meshing/netrule3.cpp \
./netgen-mesher/netgen/libsrc/meshing/adfront2.cpp \
./netgen-mesher/netgen/libsrc/meshing/specials.cpp \
./netgen-mesher/netgen/libsrc/meshing/curvedelems.cpp \
./netgen-mesher/netgen/libsrc/meshing/meshfunc.cpp \
./netgen-mesher/netgen/libsrc/meshing/ruler3.cpp \
./netgen-mesher/netgen/libsrc/meshing/refine.cpp \
./netgen-mesher/netgen/libsrc/meshing/smoothing3.cpp \
./netgen-mesher/netgen/libsrc/meshing/meshing2.cpp \
./netgen-mesher/netgen/libsrc/meshing/delaunay.cpp \
./netgen-mesher/netgen/libsrc/meshing/pyramidrls.cpp \
./netgen-mesher/netgen/libsrc/meshing/improve2gen.cpp \
./netgen-mesher/netgen/libsrc/meshing/zrefine.cpp \
./netgen-mesher/netgen/libsrc/meshing/clusters.cpp \
./netgen-mesher/netgen/libsrc/meshing/quadrls.cpp \
./netgen-mesher/netgen/libsrc/meshing/boundarylayer.cpp \
./netgen-mesher/netgen/libsrc/meshing/parser3.cpp \
./netgen-mesher/netgen/libsrc/meshing/parser2.cpp \
./netgen-mesher/netgen/libsrc/meshing/meshing3.cpp \
./netgen-mesher/netgen/libsrc/meshing/global.cpp \
./netgen-mesher/netgen/libsrc/meshing/topology.cpp \
./netgen-mesher/netgen/libsrc/meshing/meshclass.cpp \
./netgen-mesher/netgen/libsrc/meshing/netrule2.cpp \
./netgen-mesher/netgen/libsrc/meshing/secondorder.cpp \
./netgen-mesher/netgen/libsrc/meshing/ruler2.cpp \
./netgen-mesher/netgen/libsrc/meshing/paralleltop.cpp \
./netgen-mesher/netgen/libsrc/meshing/tetrarls.cpp \
./netgen-mesher/netgen/libsrc/meshing/bcfunctions.cpp \
./netgen-mesher/netgen/libsrc/meshing/meshtool.cpp \
./netgen-mesher/netgen/libsrc/meshing/improve2.cpp \
./netgen-mesher/netgen/libsrc/meshing/meshtype.cpp \
./netgen-mesher/netgen/libsrc/meshing/geomsearch.cpp \
./netgen-mesher/netgen/libsrc/geom2d/splinegeometry.cpp \
./netgen-mesher/netgen/libsrc/geom2d/geom2dmesh.cpp \
./netgen-mesher/netgen/libsrc/geom2d/genmesh2d.cpp \
./netgen-mesher/netgen/libsrc/geom2d/spline.cpp \
./netgen-mesher/netgen/nglib/ng_stl.cpp \
./netgen-mesher/netgen/nglib/nglib.cpp \
./netgen-mesher/netgen/nglib/ng_vol.cpp
