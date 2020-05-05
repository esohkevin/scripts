#!/usr/bin/env bash

#if [[ $# == 0 ]]; then
#	echo "Usage: ./case_ex.sh <name>"
#else
echo "Enter your name"
read name
#name=$1
#   while ; do
      case $(echo $name| tr [:upper:] [:lower:]) in
	esoh) echo "gentle birth";;
	nahsh) echo "apples and oranges";;
	arnold) echo "satelite song";;
	maddow) echo "moving stream";;
	cathy) echo "galaxy world";;
	eli) echo "figure bar";;
	frank) echo "music box";;
	kesoh) echo "editor"
	*) echo "Sorry '$name' not yet in database. Available names: Arnold, Cathy, Eli, Nahsh, Esoh, Frank, Val"
      esac
#   done
#fi
# 
# #read -p "Please enter your username/userid: " ID; 
# case $(echo $name | tr [:upper:] [:lower:] ) in 
# 	eso[aA0-zZ9]|0) echo "Welcome Esoh";; 
# 	isa*c|1) echo "Welcome Isaac";; 
# 	j?|2) echo "Welcome Jo";; 
# 	"") echo "Please specify an ID";; 
# 	*) echo "User $name doesn't exist in the database";; 
# esac

# counter=0;
# while [ $counter -lt 10 ]; do 
# 	(( counter++ )); 
# 	if [ $counter -eq 7 ]; then 
# 		echo "Terminator signal received: $counter"; 
# 		break; 
# 	elif [ $counter -eq 3 ]; then 
# 		echo "Count $counter is skpped!"; 
# 		continue; 
# 	fi; 
# 	echo "Counter position: $counter"; 
# done

#while getopts ":hi:o:l:" opt; do
#        usage(){
#		printf "\nUsage: %s: -i <input> -l <logfile> -o <output>\n" $(basename $0)
#	echo """
#	Usage: case_ex.sh -i [input] -l [logfile] -o [output]
#
#		-i [input]: File to process
#		-o [output]: File to save results to
#		-l [logfile]: File to contain processing information
#
#		-h [help]: Print this help message
#	"""
#	}
#	case "$opt" in 
#		i) input=$OPTARG; echo "Input: $input";; 
#		o) output=$OPTARG; echo "Output: $output";; 
#		l) logfile=$OPTARG; echo "Logfile: $logfile";;
#		h) usage; exit 0;;
#		'?') echo "Invalid Option!"; usage;;
#		:) echo "ERROR: Option -$OPTARG requires an argument"; exit 0;; 
#	esac;
#done
#shift $(( $OPTIND-1 ))	

