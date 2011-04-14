@echo off

%~d0
cd %~dp0

call qtvars.bat
call vtkvars.bat

echo Deploying response file for header path "includes"
echo /I "%QTDIR%\include" /I "%QTDIR%\include\Qt" /I "%QTDIR%\include\QtCore" /I "%QTDIR%\include\QtGui" /I "%QTDIR%\include\QtHelp" /I "%QTDIR%\include\QtNetwork" /I "%QTDIR%\include\QtOpenGL" /I "%QTDIR%\include\QtScript" /I "%QTDIR%\include\QtSql" /I "%QTDIR%\include\QtSvg" /I "%QTDIR%\include\QtTest" /I "%QTDIR%\include\QtUiTools" /I "%QTDIR%\include\QtWebKit" /I "%QTDIR%\include\QtXml" /I "%QTDIR%\include\QtXmlPatterns" /I "%QTDIR%\src" /I "%VTKINCDIR%" > qt_vtk_includes.rsp
