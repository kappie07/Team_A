#!/bin/bash

current_time=5050
while [ $current_time -le 5250 ]
do 
	python3 Running_mega_set.py
	#file="current_time.txt"
	#time=$(cat "$file")
	#addition=10
	current_time=$(( $current_time + 100 ))
done
