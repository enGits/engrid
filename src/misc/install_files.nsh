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
