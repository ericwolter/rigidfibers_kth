#!/bin/sh
arara thesis.tex
rubber-info thesis.tex

#now=$(date +"%Y%m%d_%H%M")
#cp -f thesis.pdf ~/Dropbox/theresa/Eric/${now}_thesis.pdf
