#!/usr/bin/env bash

for i in *; do
   if [[ -e $i ]] && [[ ! -d $i ]]; then
      if [ -s $i ]; then 
         echo "$i is not empty"; 
      else echo "$i is empty" && read -p 'would you like to remove it? [y/n]: ' res;
         if [[ $res == [Yy] ]]; then 
            echo -e "You have chosen to remove $i\n"; 
            sleep 1; 
            rm $i; 
            continue; 
         else echo -e "You have chosen to do nothing with the empty file\n"; 
            sleep 1; 
            continue; 
         fi; 
      fi;
   elif [[ -d $i ]]; then 
	echo "$i is a directory! skipping..."
	sleep 1;
	continue
   fi
done
