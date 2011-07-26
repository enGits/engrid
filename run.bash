#!/usr/bin/env bash
dir=`dirname $0`
export LD_LIBRARY_PATH=$dir/src/libengrid:$LD_LIBRARY_PATH
$dir/src/engrid                                                                                                                                                 
