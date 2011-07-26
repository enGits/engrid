;Tutorials
  CreateDirectory "$INSTDIR\tutorials"
  SetOutPath "$INSTDIR\tutorials\01_Throttle"
  SetOverwrite try
  File "${SRC_ROOT}\..\tutorials\01_Throttle\Throttle.msh"
  SetOutPath "$INSTDIR\tutorials\02_Sphere"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S0.blend"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S1.blend"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S1.egc"
  File "${SRC_ROOT}\..\tutorials\02_Sphere\Sphere_S1.egc.vtu"
