;Tutorials
  CreateDirectory "$INSTDIR\tutorials"
  SetOutPath "$INSTDIR\tutorials\01_Damper"
  SetOverwrite try
  File "${SRC_ROOT}\..\tutorials\01_Damper\Throttle.msh"
  File "${SRC_ROOT}\..\tutorials\01_Damper\Damper.stl"
  File "${SRC_ROOT}\..\tutorials\01_Damper\Damper.blend"

  SetOutPath "$INSTDIR\tutorials\02_Sphere"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S0.blend"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S1.blend"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S1.egc"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S1.egc.vtu"

  SetOutPath "$INSTDIR\tutorials\03_TwoCubes"
  File "${SRC_ROOT}\..\tutorials\03_TwoCubes\TwoCubes_S0.begc"
  File "${SRC_ROOT}\..\tutorials\03_TwoCubes\TwoCubes_S0.egc"
  File "${SRC_ROOT}\..\tutorials\03_TwoCubes\TwoCubes_S0.egc.vtu"
  File "${SRC_ROOT}\..\tutorials\03_TwoCubes\TwoCubes_S0.blend"
  File "${SRC_ROOT}\..\tutorials\03_TwoCubes\TwoCubes_S0.egc.geo.vtu"

;Blender scripts
  CreateDirectory "$INSTDIR\blender_scripts"
  SetOutPath "$INSTDIR\blender_scripts"
  SetOverwrite try

  CreateDirectory "$INSTDIR\blender_scripts\2.4"
  SetOutPath "$INSTDIR\blender_scripts\2.4"
  File "${SRC_ROOT}\blender_scripts\2.4\add_mesh_cylinderintersection.py"
  File "${SRC_ROOT}\blender_scripts\2.4\engrid_export.py"
  File "${SRC_ROOT}\blender_scripts\2.4\engrid_import.py"
  File "${SRC_ROOT}\blender_scripts\2.4\gmsh2_export.py"
  File "${SRC_ROOT}\blender_scripts\2.4\gmsh2_import.py"

  CreateDirectory "$INSTDIR\blender_scripts\2.59"
  SetOutPath "$INSTDIR\blender_scripts\2.59"
  File "${SRC_ROOT}\blender_scripts\2.59\io_export_engrid.py"
  File "${SRC_ROOT}\blender_scripts\2.59\io_import_engrid.py"

  CreateDirectory "$INSTDIR\blender_scripts\2.63a"
  SetOutPath "$INSTDIR\blender_scripts\2.63a"
  File "${SRC_ROOT}\blender_scripts\2.63a\io_export_engrid.py"
  File "${SRC_ROOT}\blender_scripts\2.63a\io_import_engrid.py"

;Python module scripts
  CreateDirectory "$INSTDIR\pymodules"
  SetOutPath "$INSTDIR\pymodules"
  SetOverwrite try

  CreateDirectory "$INSTDIR\pymodules"
  SetOutPath "$INSTDIR\pymodules"
  File "${SRC_ROOT}\pymodules\EngitsPyOcc.py"
  File "${SRC_ROOT}\pymodules\pymged.py"
