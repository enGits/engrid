@echo off

call %~dp0\qtvars.bat
call %~dp0\vtkvars.bat

@echo on

if not exist "%2" mkdir "%2"

@echo off

if "%1"=="Debug" goto debug
if "%1"=="Release" goto release

goto end


:debug
@echo on
copy "%QTDIR%\lib\Qt3Supportd4.lib" "%2"
copy "%QTDIR%\lib\QtAssistantClientd4.lib" "%2"
copy "%QTDIR%\lib\QtCLucened4.lib" "%2"
copy "%QTDIR%\lib\QtCored4.lib" "%2"
copy "%QTDIR%\lib\QtDesignerComponentsd4.lib" "%2"
copy "%QTDIR%\lib\QtDesignerd4.lib" "%2"
copy "%QTDIR%\lib\QtGuid4.lib" "%2"
copy "%QTDIR%\lib\QtHelpd4.lib" "%2"
copy "%QTDIR%\lib\QtNetworkd4.lib" "%2"
copy "%QTDIR%\lib\QtOpenGLd4.lib" "%2"
copy "%QTDIR%\lib\QtScriptToolsd4.lib" "%2"
copy "%QTDIR%\lib\QtScriptd4.lib" "%2"
copy "%QTDIR%\lib\QtSqld4.lib" "%2"
copy "%QTDIR%\lib\QtSvgd4.lib" "%2"
copy "%QTDIR%\lib\QtTestd4.lib" "%2"
copy "%QTDIR%\lib\QtUiToolsd.lib" "%2"
copy "%QTDIR%\lib\QtXmlPatternsd4.lib" "%2"
copy "%QTDIR%\lib\QtXmld4.lib" "%2"
copy "%QTDIR%\lib\qtmaind.lib" "%2"
copy "%VTKLIBDIR_DEBUG%\*.lib" "%2"
@echo off
goto end


:release
@echo on
copy "%QTDIR%\lib\Qt3Support4.lib" "%2"
copy "%QTDIR%\lib\QtAssistantClient4.lib" "%2"
copy "%QTDIR%\lib\QtCLucene4.lib" "%2"
copy "%QTDIR%\lib\QtCore4.lib" "%2"
copy "%QTDIR%\lib\QtDesigner4.lib" "%2"
copy "%QTDIR%\lib\QtDesignerComponents4.lib" "%2"
copy "%QTDIR%\lib\QtGui4.lib" "%2"
copy "%QTDIR%\lib\QtHelp4.lib" "%2"
copy "%QTDIR%\lib\QtNetwork4.lib" "%2"
copy "%QTDIR%\lib\QtOpenGL4.lib" "%2"
copy "%QTDIR%\lib\QtScript4.lib" "%2"
copy "%QTDIR%\lib\QtScriptTools4.lib" "%2"
copy "%QTDIR%\lib\QtSql4.lib" "%2"
copy "%QTDIR%\lib\QtSvg4.lib" "%2"
copy "%QTDIR%\lib\QtTest4.lib" "%2"
copy "%QTDIR%\lib\QtUiTools.lib" "%2"
copy "%QTDIR%\lib\QtXml4.lib" "%2"
copy "%QTDIR%\lib\QtXmlPatterns4.lib" "%2"
copy "%QTDIR%\lib\qtmain.lib" "%2"
copy "%VTKLIBDIR%\*.lib" "%2"
@echo off
goto end

:end
