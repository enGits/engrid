#!/usr/bin/env bash
dir=`dirname $0`
export LD_LIBRARY_PATH=$dir/src/libengrid:$dir/src/netgen_svn:$LD_LIBRARY_PATH
$dir/src/engrid                                                                                                                                                 
