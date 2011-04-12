
rem source the environment variables for the chosen qt installation
call %~dp0\qtvars.bat

uic.exe %5%1 -o %5ui_%2.h
moc.exe -I"%QTDIR%\qt\include" -I"%QTDIR%\qt\include\Qt" -I"%QTDIR%\qt\include\QtCore" -I"%QTDIR%\qt\include\QtGui" -I"%QTDIR%\qt\include\QtHelp" -I"%QTDIR%\qt\include\QtNetwork" -I"%QTDIR%\qt\include\QtOpenGL" -I"%QTDIR%\qt\include\QtScript" -I"%QTDIR%\qt\include\QtSql" -I"%QTDIR%\qt\include\QtSvg" -I"%QTDIR%\qt\include\QtTest" -I"%QTDIR%\qt\include\QtUiTools" -I"%QTDIR%\qt\include\QtWebKit" -I"%QTDIR%\qt\include\QtXml" -I"%QTDIR%\qt\include\QtXmlPatterns" -I"%QTDIR%\qt\src" -D%4 -D_WINDOWS -DQT_LARGEFILE_SUPPORT -DQT_DLL -DQT_GUI_LIB -DQT_CORE_LIB -DQT_THREAD_SUPPORT -D_CRT_SECURE_NO_WARNINGS -D_UNICODE -D_MSC_VER=1500 -DWIN32 %5%2.h -o %3\moc_%2.cpp
