#!/bin/bash

# Default Options
dim="-2d"
nxb=32
nyb=32
nzb=32
opt="-opt"
grd="+rg"
dir="RBCAV"
flg=""
hpc="manwe.ddns.net"

# Default Advanced Options
advn=false
comp=false
ncmp=4
over=false
coll=false
h5pp=false

# Output Usage Information
usage ()
{
    echo "usage: ./Rayleigh_Benard.sh [[[-nxb 24] [-dim 2] [-grd rg] ...] | [-h]]"
    echo ""
}

overwrite() { echo -e "\r\033[1A\033[0K$@"; }

# Process Commandline Options
while [ "$1" != '' ]; do
    case $1 in
        -nnb | --block_xyz )    shift
                                nxb=$1
                                nyb=$1
                                nzb=$1
                                ;;
        -nxb | --block_x )      shift
                                nxb=$1
                                ;;
        -nyb | --block_y )      shift
                                nyb=$1
                                ;;
        -nzb | --block_z )      shift
                                nzb=$1
                                ;;
        -dim | --dimension )    shift
                                dim="-$1d"
                                ;;
        -opt | --optimize )     shift
                                opt="-$1"
                                ;;
        -grd | --grid )         shift
                                grd="+$1"
                                ;;
        -dir | --name )         shift
                                dir=$i1
                                ;;
        -flg | --note )         shift
                                flg=$1
                                ;;
        -hpc | --site )         shift
                                hpc=$1
                                ;;
        -advn | --advanced )    advn=true
                                ;;
        -comp | --compile )     shift
                                comp=true
                                ncmp=$1
                                ;;
        -over | --overwrite )   over=true
                                ;;
        -h5pp | --h5_parallel ) h5pp=true
                                ;;
        -coll | --collection )  coll=true
                                ;;
        -all  | --advn_coll  )  coll=true
                                advn=true
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
        esac
        shift
done

# Identify Simulation Directory
if [ "${grd}" == "+rg" ]; then
    simdir="INavierStokes/zoso/Rayleigh_Benard_Stretch/"
elif [ "${grd}" == "+ug" ]; then
    simdir="INavierStokes/zoso/Rayleigh_Benard/"
elif [ "${grd}" == "+pg" ]; then
    simdir="INavierStokes/zoso/Poisson/"
    grd="+rg"
    dir="POISS"
fi

# Identify Parallel or Serial IO
if ${h5pp}; then
	h5f=" +parallelIO"
else
	h5f=""
fi

# Create Setup Command From User Options
if [ "${dim}" == "-2d" ]; then
    objdir="${grd:1}${dir}_2D${nxb}-${nyb}"
    setcom="./setup ${simdir} ${dim} -auto ${opt}${h5f} -nxb=${nxb} -nyb=${nyb} ${grd} -objdir=${objdir}${flg} -site=${hpc}"
else
    objdir="${grd:1}${dir}_3D${nxb}-${nyb}-${nzb}"  
    setcom="./setup ${simdir} ${dim} -auto ${opt}${h5f} -nxb=${nxb} -nyb=${nyb} -nzb=${nzb} ${grd} -objdir=${objdir}${flg} -site=${hpc}"
fi

# Run the Setup Script
echo ""
echo "-------------------------------------------------------------"
echo "--  Create FLASH Compilation Directory -- Rayleigh Benard  --"
echo "-------------------------------------------------------------"
echo ""
echo "Using the following FLASH setup command:"
echo "----------------------------------------"
echo "${setcom}"
echo ""

# Advanced Options Setup Logic
if $advn && ! $coll; then

    # Do not overwrite setup directory
    if [ -d "${objdir}" ] &&  ! $over; then
        echo ""
        echo "Setup Directory ${objdir} already exists -- use -over to overwrite"
        echo ""

    # Overwrite setup directory
    else
      ${setcom}
      echo "----------------------------------------"
      echo ""    

    fi

    # Compile the setup directory
    if $comp; then
        echo "Compile the FLASH setup directory w/ ${ncmp} cores:"
        echo "-----------------------------------------------------"
        cd ${objdir}
         
        screen -d -m -S makeFLASH bash -c "make -j ${ncmp} > out"
        count=1
        visual="."
        while screen -list | grep -q "makeFLASH"
        do
            overwrite "Compiling FLASH code ${visual}"
            count=$(($count + 1))
            visual="${visual} ."
            if (( $count >= 15 )); then
                count=1
                visual="."
            fi
            sleep 1
        done

        grep -i success out
        cd ..
        echo "-----------------------------------------------------"
        echo ""
        echo ""

    # Do not compile the setup directory
    else
        echo ""

    fi

# Process Batch of Setup Commands
elif $advn && $coll; then
    echo ""
    echo ""
    echo "//////////////////////////////////////////////////////////"
    echo "//                                                      //"
    echo "// --  -- Batch Process FLASH Setup Directories  -- --  //"
    echo "//                                                      //"
    echo "//////////////////////////////////////////////////////////"
    echo ""
    echo ""
    echo ""
    echo ""

    while IFS= read -r line
    do
        if [ "${line:0:1}" = "-" ]; then
            if $over; then
                ./Rayleigh_Benard.sh $line -over
            else
                ./Rayleigh_Benard.sh $line
           fi 
        fi
    done < Rayleigh_Benard_Collection

# Defualt Options Setup Logic
else
    ${setcom}
    echo "----------------------------------------"
    echo ""    

fi

# Clear Changes to bin/ and lib/
git checkout bin/*
git checkout lib/*
