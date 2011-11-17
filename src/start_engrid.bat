@echo off
%~d0
cd %~dp0

call windows\mingw32\qtvars.bat
call windows\mingw32\vtkvars.bat
call windows\mingw32\mingwvars.bat

rem set version=debug
set version=release
set PATH=%CD%\%version%;%CD%\libengrid\%version%;%CD%\netgen_svn\%version%;%PATH%

engrid
pause
