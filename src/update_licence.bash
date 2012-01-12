#!/bin/bash
export LD_LIBRARY_PATH=./libengrid:./netgen_svn:$LD_LIBRARY_PATH
./engrid -appendlic *.cpp libengrid/*.h libengrid/*.cpp scripts/*.sh scripts/*.bash scripts/*.py ../*.bash

