@echo off
%~d0
cd %~dp0

set QTDIR=Q:\2009.02\qt
set VTKLIBDIR=P:\OpenAPIs\VTK 5.4\VTK_msys\bin

set PATH=%QTDIR%\bin;Q:\2009.02\mingw\bin;Q:\2009.02\bin;%VTKLIBDIR%;%PATH%

set version=release
set PATH=%CD%\%version%;%CD%\libengrid\%version%;%PATH%

engrid
pause
