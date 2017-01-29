#!/bin/bash
#This is a small script that compares the output files of FLASH runs
#produced using MeshDump I/O implementation.
#I use the script to ensure that the Chombo UG Grid implementation
#gives identical results to the FLASH UG Grid implementation.

if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` dir1 dir2"
  exit ${E_BADARGS}
else
    dir1="$1"
    dir2="$2"
fi
err=0


echo "Comparing Mesh dump files between directories ${dir1} and ${dir2}"
for f in $(find ${dir1} -name "CENTER*.bin")
  do
  outfile=$(basename "${f}")
  diff -q "${dir1}/${outfile}" "${dir2}/${outfile}"
  if [ $? == 0 ]
      then
      rtn="SUCCESS"
  else
      rtn="FAILURE"
      err=$(expr ${err} + 1)
  fi
  echo "${outfile} - ${rtn}"
done


if [ ${err} == 0 ]
then
    echo "SUCCESS - all files match"
else
    echo "FAILURE - ${err} files do not match"
fi
