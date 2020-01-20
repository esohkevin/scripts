#!/bin/bash

for i in {21..97} {97..21}; do
	echo -en "\e[38;5;${i}m#\e[0m"
done
echo " "



