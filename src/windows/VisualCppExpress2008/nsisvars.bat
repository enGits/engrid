
rem source the environment variables for the chosen NSIS installation
call %~dp0..\..\third_party\NSIS\nsisvars.bat

rem run any variables left to be used
%*
