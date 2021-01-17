site="manwe.ddns.net"
build="pm" 

while [ "$1" != '' ]; do
    case $1 in
        -pm | --paramesh )   shift
                             build="pm"
			     ;;
        -ug | --uniform )    shift
                             build="ug"
			     ;;
        -rg | --regular )    shift
                             build="rg"
			     ;;
        * )                  echo "Unkown Argument!"
	                     exit 1
        esac
	shift
done

if [ "${build}" == "pm" ]; then
    setup="./setup INavierStokes/zoso/LidDrivenCavity/ -2d -auto -opt +parallelIO -nxb=16 -nyb=16 -maxblocks=3000 -gridinterpolation=native +pm4dev PfftSolver=HomBcTrigSolver -objdir=pmLIDCAV -site=${site}"
elif [ "${build}" == "ug" ]; then
    setup="./setup INavierStokes/zoso/LidDrivenCavity/ -2d -auto -opt -nxb=64 -nyb=64 +ug -objdir=ugLIDCAV -site=${site}"
fi

${setup}

# Clear Changes to bin/ and lib/
git checkout bin/*
git checkout lib/*
