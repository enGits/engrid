
rem source the environment variables for the chosen qt installation
IF %1 == Win32 call %~dp0\..\..\third_party\Qt\bin\qtvars.bat > NUL:
IF %1 == x64 call %~dp0\..\..\third_party64\Qt\bin\qtvars.bat > NUL:

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
