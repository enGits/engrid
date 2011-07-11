
rem source the environment variables for the chosen qt installation
call %~dp0\..\..\third_party\Qt\bin\qtvars.bat > NUL:

rem run any variables left to be used
%*
