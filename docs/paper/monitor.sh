#!/bin/bash

while true
do
   ATIME=`stat -c %Z thesis_template.tex`

   if [[ "$ATIME" != "$LTIME" ]]
   then
       echo "Rebuilding..."
       ./build.sh
       LTIME=$ATIME
   fi
   sleep 1
done
