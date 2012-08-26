;VTK libs
  Delete "$INSTDIR\QVTK.DLL"
  Delete "$INSTDIR\VTKALGLIB.DLL"
  Delete "$INSTDIR\VTKCOMMON.DLL"
  Delete "$INSTDIR\VTKDICOMPARSER.DLL"
  Delete "$INSTDIR\VTKEXOIIC.DLL"
  Delete "$INSTDIR\VTKEXPAT.DLL"
  Delete "$INSTDIR\VTKFILTERING.DLL"
  Delete "$INSTDIR\VTKFREETYPE.DLL"
  Delete "$INSTDIR\VTKFTGL.DLL"
  Delete "$INSTDIR\VTKGRAPHICS.DLL"
  Delete "$INSTDIR\VTKHYBRID.DLL"
  Delete "$INSTDIR\VTKIMAGING.DLL"
  Delete "$INSTDIR\VTKINFOVIS.DLL"
  Delete "$INSTDIR\VTKIO.DLL"
  Delete "$INSTDIR\VTKJPEG.DLL"
  Delete "$INSTDIR\VTKLIBXML2.DLL"
  Delete "$INSTDIR\VTKMETAIO.DLL"
  Delete "$INSTDIR\VTKNETCDF.DLL"
  Delete "$INSTDIR\VTKRENDERING.DLL"
  Delete "$INSTDIR\VTKSYS.DLL"
  Delete "$INSTDIR\VTKTIFF.DLL"
  Delete "$INSTDIR\VTKVERDICT.DLL"
  Delete "$INSTDIR\VTKVIEWS.DLL"
  Delete "$INSTDIR\VTKWIDGETS.DLL"
  Delete "$INSTDIR\ZLIB1.DLL"
  Delete "$INSTDIR\LIBPNG15.DLL"
  Delete "$INSTDIR\LIBMPI.DLL"
  Delete "$INSTDIR\LIBMPI_CXX.DLL"
  Delete "$INSTDIR\LIBOPEN-PAL.DLL"
  Delete "$INSTDIR\LIBOPEN-RTE.DLL"
  Delete "$INSTDIR\VTKNETCDF_CXX.DLL"
  Delete "$INSTDIR\VTKPARALLEL.DLL"
  Delete "$INSTDIR\VPIC.DLL"
  Delete "$INSTDIR\COSMO.DLL"

;Qt libs
  Delete "$INSTDIR\QTCORE4.DLL"
  Delete "$INSTDIR\QTGUI4.DLL"
  Delete "$INSTDIR\QTXML4.DLL"
  Delete "$INSTDIR\QTWEBKIT4.DLL"
  Delete "$INSTDIR\phonon4.dll"
  Delete "$INSTDIR\QtXmlPatterns4.dll"
  Delete "$INSTDIR\QtNetwork4.dll"
  Delete "$INSTDIR\QtSql4.dll"

  Delete "$INSTDIR\plugins\sqldrivers\qsqlite4.dll"
  RMDir "$INSTDIR\plugins\sqldrivers"
  RMDir "$INSTDIR\plugins"

;BRLCAD libs
  Delete "$INSTDIR\brlcad.dll"
  Delete "$INSTDIR\openNURBS.dll"
;BRLCAD legal docs
  Delete "$INSTDIR\COPYING.BRLCAD"
  Delete "$INSTDIR\readme.BRLCAD.txt"

;MinGW runtime
!ifdef USE_mingw32
  Delete "$INSTDIR\MINGWM10.DLL"
!endif

;MSVC 2008 Runtimes
!ifdef USE_VisualCppExpress2008
  Delete "$INSTDIR\msvcm90.dll"
  Delete "$INSTDIR\msvcp90.dll"
  Delete "$INSTDIR\msvcr90.dll"
  Delete "$INSTDIR\Microsoft.VC90.CRT.manifest"
!endif
