#!/bin/bash
if [[ -z "$3" ]]; then
	echo "Usage: $0 outfile scale fps" 
else
	ffmpeg -r $3 -i img/%04d.png -q:v 2 -pix_fmt yuv444p -vf scale=iw/$2:-1 $1
fi
