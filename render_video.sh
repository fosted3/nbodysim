#!/bin/sh
if [[ -z "$1" ]]; then
	echo "Usage: $0 outfile"
else
	ffmpeg -r 30 -i img/%04d.png -q:v 2 -pix_fmt yuv444p $1
fi
