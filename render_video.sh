#!/bin/sh
if [[ -z "$1" ]]; then
	echo "Usage: $0 outfile"
else
	ffmpeg -r 30 -i img/%04d.png -qscale 2 $1
fi
