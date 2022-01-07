#!/bin/sh
set -e
cmake -B build
cmake --build build
echo "------baseline------"
build/main
echo "------   my   ------"
build/work
