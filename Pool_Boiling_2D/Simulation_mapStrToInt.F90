!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_mapStrToInt(str,key,map)
    use Grid_interface, only: Grid_parseNonRep
    use Driver_interface, only: Driver_getMype, Driver_getNumProcs
    use RuntimeParameters_interface, only: RuntimeParameters_get
    implicit none
#include "constants.h"
#include "Flash.h"
    character(len=*), intent(in) :: str
    integer, intent(out) :: key 
    integer, intent(in) :: map

    integer, parameter :: locunk1(0:NONREP_COUNT) = NONREP_LOCUNK1
    character(*), parameter :: rpcount_flat = NONREP_RPCOUNT_FLAT
    integer, parameter :: rpcount_start(NONREP_COUNT+1) = NONREP_RPCOUNT_START
    
    integer :: mesh, meshes
    character(len=MAX_STRING_LENGTH) :: strlwr
    integer :: nonrep, glob, nglob
    
    call Driver_getMype(MESH_ACROSS_COMM, mesh)
    call Driver_getNumProcs(MESH_ACROSS_COMM, meshes)
    key = NONEXISTENT
    strlwr = str
    call makeLowercase(strlwr)
    
    call Grid_parseNonRep(strlwr(1:len(str)), nonrep, glob)
    if(nonrep .gt. 0) then
        call RuntimeParameters_get(rpcount_flat(rpcount_start(nonrep):rpcount_start(nonrep+1)-1), nglob)
        if(glob .gt. nglob .or. mesh .ne. NONREP_MESHOFGLOB(glob, meshes)) return ! NONEXISTENT
        key = locunk1(nonrep)-1 + NONREP_GLOB2LOC(glob, mesh, meshes)
        return
    end if
    
    select case(strlwr)
    case("curv")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10
        end select
    case("velo")
        select case(map)
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  9)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  9
        end select
    case("xi")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6
        end select
    case("veli")
        select case(map)
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  8)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  8
        end select
    case("velc")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  7)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  7
        end select
    case("mflg")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15
        end select
    case("rhst")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34
        end select
    case("sigp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36
        end select
    case("pres")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32
        end select
    case("velz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47
        end select
    case("vely")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46
        end select
    case("velx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45
        end select
    case("aajunk")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8
        end select
    case("tnlq")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41
        end select
    case("dupo")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3
        end select
    case("rhso")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33
        end select
    case("preo")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31
        end select
    case("dfun")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12
        end select
    case("ptes")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4
        end select
    case("nrmy")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24
        end select
    case("pfun")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30
        end select
    case("xipo")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7
        end select
    case("nrmx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23
        end select
    case("rold")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35
        end select
    case("sigm")
        select case(map)
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  6)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  6
        end select
    case("rtes")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5
        end select
    case("told")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43
        end select
    case("flxmc")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1
        end select
    case("tmic")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40
        end select
    case("con_angle")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1
        end select
    case("delp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11
        end select
    case("tnvp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42
        end select
    case("dust")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13
        end select
    case("divv")
        select case(map)
        case(((MAPBLOCK_SCRATCH * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_SCRATCH * MAPBLOCKSIZE)+  1
        end select
    case("divu")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2
        end select
    case("mgw9")
        select case(map)
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  2
        end select
    case("tvis")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44
        end select
    case("mdot")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14
        end select
    case("temp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39
        end select
    case("smrh")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38
        end select
    case("visc")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48
        end select
    case("omgy")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28
        end select
    case("omgx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27
        end select
    case("omgz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29
        end select
    case("rhds")
        select case(map)
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  5)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  5
        end select
    case("smhv")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37
        end select
    case("mgw1")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16
        end select
    case("nrmz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25
        end select
    case("mgw3")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18
        end select
    case("mgw2")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17
        end select
    case("mgw5")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20
        end select
    case("mgw4")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19
        end select
    case("mgw7")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22
        end select
    case("mgw6")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21
        end select
    case("rh1f")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  3
        end select
    case("mgw8")
        select case(map)
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  1
        end select
    case("rh2f")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3
        case(((MAPBLOCK_FACES * MAPBLOCKSIZE)+  4)/MAPBLOCKSIZE); key = (MAPBLOCK_FACES * MAPBLOCKSIZE)+  4
        end select
    case("alph")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9
        end select
    case("omgm")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26
        end select
    end select

    if(key .ne. NONEXISTENT) key = mod(key,MAPBLOCKSIZE)
end subroutine Simulation_mapStrToInt
