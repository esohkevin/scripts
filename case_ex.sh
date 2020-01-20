#!/usr/bin/env bash

if [[ $# == 0 ]]; then
	echo "Usage: ./case_ex.sh <name>"
else
#echo "Enter your name"
#read name
name=$1
#   while ; do
      case $(echo $name| tr [:upper:] [:lower:]) in
	esoh) echo "Gentle-Kind-Generous-Handsome birth";;
	nahsh) echo "Musician";;
	arnold) echo "Strongest man";;
	val) echo "Beauty";;
	cathy) echo "Strong woman";;
	eli) echo "Cool girl";;
	frank) echo "Genious";;
	*) echo "Sorry '$name' not yet in database. Available names: Arnold, Cathy, Eli, Nahsh, Esoh, Frank, Val"
      esac
#   done
fi

#read -p "Please enter your username/userid: " ID; 
case $(echo $name | tr [:upper:] [:lower:] ) in 
	esoh|0) echo "Welcome Esoh";; 
	isaac|1) echo "Welcome Isaac";; 
	jo|2) echo "Welcome Jo";; 
	"") echo "Please specify an ID";; 
	*) echo "User $name doesn't exist in the database";; 
esac

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
