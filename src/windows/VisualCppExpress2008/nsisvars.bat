
rem source the environment variables for the chosen NSIS installation
IF %1 == win32 call %~dp0..\..\third_party\NSIS\nsisvars.bat
IF %1 == x64 call %~dp0..\..\third_party64\NSIS\nsisvars.bat

rem throw the first parameter away
shift
set params=%1
:loop
shift
if [%1]==[] goto afterloop
set params=%params% %1
goto loop
:afterloop

rem run any variables left to be used
%params%
