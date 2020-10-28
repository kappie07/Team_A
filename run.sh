#!/bin/sh

ffmpeg -framerate 10 -pattern_type glob -i "merge_plots/snap*.png" -s:v 640x480 -c:v libx264 -profile:v high -level 4.0 -crf 10 -tune animation -preset slow -pix_fmt yuv420p -r 25 -threads 0 -f mp4 galaxy_merge.mp4

ffmpeg -framerate 10 -pattern_type glob -i "star_plots/snapshots/snap*.png" -s:v 640x480 -c:v libx264 -profile:v high -level 4.0 -crf 10 -tune animation -preset slow -pix_fmt yuv420p -r 25 -threads 0 -f mp4 galaxy_merge_tracers.mp4
