;Tutorials
  Delete "$INSTDIR\tutorials\01_Damper\Throttle.msh"
  Delete "$INSTDIR\tutorials\01_Damper\Damper.stl"
  Delete "$INSTDIR\tutorials\01_Damper\Damper.blend"
  RMDir "$INSTDIR\tutorials\01_Damper"

  Delete "$INSTDIR\tutorials\02_Sphere\Sphere_S1.egc.vtu"
  Delete "$INSTDIR\tutorials\02_Sphere\Sphere_S1.egc"
  Delete "$INSTDIR\tutorials\02_Sphere\Sphere_S1.blend"
  Delete "$INSTDIR\tutorials\02_Sphere\Sphere_S0.blend"
  RMDir "$INSTDIR\tutorials\02_Sphere"

  Delete "$INSTDIR\tutorials\03_TwoCubes\TwoCubes_S0.begc"
  Delete "$INSTDIR\tutorials\03_TwoCubes\TwoCubes_S0.egc"
  Delete "$INSTDIR\tutorials\03_TwoCubes\TwoCubes_S0.egc.vtu"
  Delete "$INSTDIR\tutorials\03_TwoCubes\TwoCubes_S0.blend"
  Delete "$INSTDIR\tutorials\03_TwoCubes\TwoCubes_S0.egc.geo.vtu"
  RMDir "$INSTDIR\tutorials\03_TwoCubes"

  RMDir "$INSTDIR\tutorials"

;Blender scripts
  Delete "$INSTDIR\blender_scripts\2.4\add_mesh_cylinderintersection.py"
  Delete "$INSTDIR\blender_scripts\2.4\engrid_export.py"
  Delete "$INSTDIR\blender_scripts\2.4\engrid_import.py"
  Delete "$INSTDIR\blender_scripts\2.4\gmsh2_export.py"
  Delete "$INSTDIR\blender_scripts\2.4\gmsh2_import.py"
  RMDir "$INSTDIR\blender_scripts\2.4"

  Delete "$INSTDIR\blender_scripts\2.59\io_export_engrid.py"
  Delete "$INSTDIR\blender_scripts\2.59\io_import_engrid.py"
  RMDir "$INSTDIR\blender_scripts\2.59"

  Delete "$INSTDIR\blender_scripts\2.63a\io_export_engrid.py"
  Delete "$INSTDIR\blender_scripts\2.63a\io_import_engrid.py"
  RMDir "$INSTDIR\blender_scripts\2.63a"

  RMDir "$INSTDIR\blender_scripts"

;Python module scripts
  Delete  "$INSTDIR\pymodules\EngitsPyOcc.py"
  Delete  "$INSTDIR\pymodules\pymged.py"

  RMDir "$INSTDIR\pymodules"
