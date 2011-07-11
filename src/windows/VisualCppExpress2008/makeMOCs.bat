@echo off

rem source the environment variables for the chosen qt installation
call %~dp0\qtvars.bat

set TARGETFILE=%3\moc_%2.cpp

moc.exe -D%4 -D_WINDOWS -DQT_LARGEFILE_SUPPORT -DQT_DLL -DQT_GUI_LIB -DQT_CORE_LIB -DQT_THREAD_SUPPORT -D_CRT_SECURE_NO_WARNINGS -D_UNICODE -D_MSC_VER=1500 -DWIN32 %5%2.h -o %TARGETFILE% 2> NUL:

IF NOT EXIST %TARGETFILE% goto QTLESS

FOR %%? IN (%TARGETFILE%) DO (
    if %%~z? == 0 del %TARGETFILE%
)

IF EXIST %TARGETFILE% goto QTEND

:QTLESS
echo %5%2.h: Skipping file.
exit 0

:QTEND
echo Generated MOC file: %TARGETFILE%
exit 0
