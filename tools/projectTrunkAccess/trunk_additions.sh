#!/bin/bash
#Creates a trunk_additions file which contains the directories
#that appear in FLASH trunk but not in FLASH release.


#Parameters.
#-----------------------------------------------------------------------
#flash_trunk and flash_release must be full file names (i.e. including path).
flash_trunk="/home/cdaley/custom_tarball/trunk"
flash_release="/home/cdaley/custom_tarball/FLASH3.3_release"
work_dir="$(pwd)"


#Do not edit below this line.
#-----------------------------------------------------------------------
trees=( "${flash_trunk}" "${flash_release}" )
logs=( $(basename "${trees[0]}").out $(basename "${trees[1]}").out )
trunk_additions="trunk_additions.txt"
start_dir="$(pwd)"

#Record the directories within each source tree.
for i in 0 1
do
  cd "${trees[$i]}"
  find source -type d | grep -Fv .svn | sort > "${work_dir}"/"${logs[$i]}"
done

#Suppress lines that are unique to the release or that appear in both files.
#This is my list of trunk additions.
cd "${work_dir}"
comm -23 "${logs[0]}" "${logs[1]}" > "${trunk_additions}"
rm -f "${logs[0]}" "${logs[1]}"
cd "${start_dir}"
