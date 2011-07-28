
rem source the environment variables for the chosen NSIS installation
call %~dp0\nsisvars.bat %7

echo Processing NSIS file: %1
makensis.exe /DOUTPUTFOLDER="%2" /D%3 /DBUILDVERSION=%4 /DPRODUCT_VERSION=%5 /DPRODUCT_VERSIONx4=%6 %1
