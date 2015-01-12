#!/bin/bash

mkdir crop/
mkdir png/
rm -f crop/*
rm -f png/*
rm -f *.mp4

i=0
for f in state_*.pdf;
do
  #(pdfcrop --margins "10 10 10 10" "$f" "crop/$f") &
  (cp "$f" "crop/$f") &

  i=$((i + 1));
  if (( $i % 8 == 0 )); then wait; fi

done
wait

i=0
for f in state_*.pdf;
do
  (convert -bordercolor White -border 4%x4% -density 288 -transparent none "crop/$f" "png/$f.png") &

  i=$((i + 1));
  if (( $i % 8 == 0 )); then wait; fi

done
wait

ffmpeg -framerate 10 -i png/state_%05d.pdf.png -pix_fmt yuv420p -y anim_state.mp4
