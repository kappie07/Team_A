#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00

cd /data/s1504312/SMA_project/Team_A

source /usr/alice/python/3.6.0/setup.sh
python3 Running_bridge.py
