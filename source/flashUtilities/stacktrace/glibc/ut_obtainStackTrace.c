/* 
   Obtain a stack trace at any time during an application run.  The stack 
   trace is printed to a process unique file.  Note that the stack trace is 
   only useful when there is debugging information in the binary.

   Call init_stacktrace() to obtain a file descriptor to a process unique file.
   Call finalise_stacktrace() to close the file associated with file descriptor.
   Call dump_stacktrace() to write to the file associated with file descriptor.

   Use NFRAMES to control the number of stack frames we output:

   NFRAMES 4
      Example output when binary compiled using -g option:
      [0x94e97d]
      [0x91d6cd]
      [0x5eedef]
      [0x58b5d8]

      (Obtain File:Line information using: addr2line 0x5eedef -e flash3)

   
      Example output when binary compiled using -ggdb option:
      ./flash3(dump_stacktrace_and_int_+0x56)[0x5f0e35]
      ./flash3(mpi_ssend_+0x32)[0x5e1169]
      ./flash3(fill_old_loc_+0x6d5)[0x4d94c5]
      ./flash3(amr_migrate_tree_data_+0x579)[0x4b325d]

   Possible usage shown at the bottom of this file.
*/

#include <execinfo.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mangle_names.h"

#define PERMS 0666
#define NFRAMES 4
static int s_fd;

void FTOC(init_stacktrace)(const int * const pMyPE)
{
  int myPE = *pMyPE;
  char baseName[] = "stacktrace.out.";
  char procExtension[10], fileName[100];

  /* Construct a unique file name */
  sprintf(procExtension, "%05d", myPE);
  strcpy(fileName, baseName);
  strcat(fileName, procExtension);

  s_fd = creat(fileName, PERMS);
  if (s_fd == -1) {
    fprintf(stderr, "Creating output file failed\n");
    exit(EXIT_FAILURE);
  }  
}


void FTOC(dump_stacktrace)()
{
  void *frameArray[NFRAMES];
  int nptrs, bytesWritten;

  nptrs = backtrace(frameArray, NFRAMES);

  /* Write the frame strings to file descriptor s_fd. */
  backtrace_symbols_fd(frameArray, nptrs, s_fd);
}


void FTOC(dump_stacktrace_and_int)(const int * const pSomeInt)
{
  void *frameArray[NFRAMES];
  int nptrs, bytesWritten;
  int someInt = *pSomeInt;
  char buf[10];
  buf[0]='\n'; buf[1]='\n';
  sprintf(&buf[2], "%06d", someInt);
  buf[8]='\n'; buf[9]='\0';

  nptrs = backtrace(frameArray, NFRAMES);

  /* Write the buf and the frame strings to file descriptor s_fd. */
  bytesWritten = write(s_fd, buf, sizeof(char)*9);
  backtrace_symbols_fd(frameArray, nptrs, s_fd);
}


void FTOC(finalise_stacktrace)()
{
  close(s_fd);
}


/*   
     Possible usage:  Suppose you want to know the MPI tags of every 
     single executed point-to-point communication.  

     1.  Create a Fortran file that contains an implementation for 
     each MPI point-to-point communications.  In the implementation 
     we add the appropriate calls to the stack trace routines:
     

subroutine MPI_INIT(IERROR)
  implicit none
  include "mpif.h"
  INTEGER :: IERROR
  integer :: myPE, ierr
  
  call PMPI_INIT(IERROR)
  call pmpi_comm_rank(MPI_COMM_WORLD, myPE, ierr)
  call init_stacktrace(myPE)
end subroutine MPI_INIT


subroutine MPI_FINALIZE(IERROR)
  implicit none
  include "mpif.h"
  INTEGER :: IERROR

  call finalise_stacktrace()
  call PMPI_Finalize(IERROR)
end subroutine MPI_Finalize


subroutine MPI_ABORT(COMM, ERRORCODE, IERROR)
  implicit none
  include "mpif.h"
  INTEGER :: COMM, ERRORCODE, IERROR

  call finalise_stacktrace()
  call PMPI_Abort(COMM, ERRORCODE, IERROR)
end subroutine MPI_ABORT


subroutine MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  implicit none
  include "mpif.h"
  REAL :: BUF(*)
  INTEGER :: COUNT, DATATYPE, DEST, TAG, COMM, IERROR 
  integer :: myPE, ierr

  call dump_stacktrace_and_int(TAG)
  call PMPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  call pmpi_comm_rank(COMM, myPE, ierr)
  call Driver_checkMPIErrorCode(myPE, IERROR)
end subroutine MPI_SEND


subroutine MPI_ISEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, REQUEST, IERROR)
  implicit none
  include "mpif.h"
  REAL :: BUF(*)
  INTEGER :: COUNT, DATATYPE, DEST, TAG, COMM, REQUEST, IERROR 
  integer :: myPE, ierr

  call dump_stacktrace_and_int(TAG)
  call PMPI_ISEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, REQUEST, IERROR)
  call pmpi_comm_rank(COMM, myPE, ierr)
  call Driver_checkMPIErrorCode(myPE, IERROR)
end subroutine MPI_ISEND


subroutine MPI_SSEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  implicit none
  include "mpif.h"
  REAL :: BUF(*)
  INTEGER :: COUNT, DATATYPE, DEST, TAG, COMM, IERROR 
  integer :: myPE, ierr  

  call dump_stacktrace_and_int(TAG)
  call PMPI_SSEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
  call pmpi_comm_rank(COMM, myPE, ierr)
  call Driver_checkMPIErrorCode(myPE, IERROR)
end subroutine MPI_SSEND


subroutine MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)
  implicit none
  include "mpif.h"
  REAL :: BUF(*)
  INTEGER :: COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS(MPI_STATUS_SIZE), IERROR 
  integer :: myPE, ierr

  call dump_stacktrace_and_int(TAG)
  call PMPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)
  call pmpi_comm_rank(COMM, myPE, ierr)
  call Driver_checkMPIErrorCode(myPE, IERROR)
end subroutine MPI_RECV


subroutine MPI_IRECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR)
  implicit none
  include "mpif.h"
  REAL :: BUF(*)
  INTEGER :: COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR 
  integer :: myPE, ierr  

  call dump_stacktrace_and_int(TAG)
  call PMPI_IRECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR)
  call pmpi_comm_rank(COMM, myPE, ierr)
  call Driver_checkMPIErrorCode(myPE, IERROR)
end subroutine MPI_IRECV


subroutine MPI_SENDRECV(SENDBUF, SENDCOUNT, SENDTYPE, DEST, SENDTAG, &
     RECVBUF, RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM, STATUS, IERROR)
  implicit none
  include "mpif.h"
  REAL :: SENDBUF(*), RECVBUF(*)
  INTEGER :: SENDCOUNT, SENDTYPE, DEST, SENDTAG, RECVCOUNT, RECVTYPE, &
       SOURCE, RECVTAG, COMM, STATUS(MPI_STATUS_SIZE), IERROR 
  integer :: myPE, ierr

  call dump_stacktrace_and_int(SENDTAG)
  call PMPI_SENDRECV(SENDBUF, SENDCOUNT, SENDTYPE, DEST, SENDTAG, &
     RECVBUF, RECVCOUNT, RECVTYPE, SOURCE, RECVTAG, COMM, STATUS, IERROR)
  call pmpi_comm_rank(COMM, myPE, ierr)
  call Driver_checkMPIErrorCode(myPE, IERROR)
end subroutine MPI_SENDRECV


subroutine MPI_WAIT(REQUEST, STATUS, IERROR)
  implicit none
  include "mpif.h"  
  integer REQUEST, STATUS(MPI_STATUS_SIZE), IERROR
  integer :: myPE, ierr

  call PMPI_WAIT(REQUEST, STATUS, IERROR)
  call Driver_checkMPIErrorCode(-1, IERROR)
end subroutine MPI_WAIT


subroutine MPI_WAITANY(COUNT, ARRAY_OF_REQUESTS, INDEX, STATUS, IERROR)
  implicit none
  include "mpif.h"  
  INTEGER COUNT, ARRAY_OF_REQUESTS(*), INDEX, STATUS(MPI_STATUS_SIZE), IERROR 
  integer :: myPE, ierr

  call PMPI_WAITANY(COUNT, ARRAY_OF_REQUESTS, INDEX, STATUS, IERROR)
  call Driver_checkMPIErrorCode(-1, IERROR)
end subroutine MPI_WAITANY



     2.  Post-process the output files to obtain the MPI tags for each 
     point-to-point communication.  Below is a python script to print 
     every MPI tag along with the function (file and line number) making 
     the MPI call:


#!/usr/bin/python
import commands

def get_stack_trace(fnAddress, pwdStr):
    addr2lineStr = "addr2line " + fnAddress + " -e flash3"
    fullName = commands.getoutput(addr2lineStr)
    fnName = fullName.replace(pwdStr+'/','')
    return fnName   


count = 0
dict_tag_fn_map = {}

fileHandle = open('stacktrace.out.00000')
fileList = fileHandle.readlines()

for fileLine in fileList:

    count = count + 1
    line = fileLine.strip() #Remove '\n'.
    if count == 3:
        tag = line   #3rd line down stores the tag.
    if count == 6:   #6th line down stores the MPI call address.
        caller = line.strip('[]')  #Remove the square brackets.    
    if count == 7:   #7th line down stores the parent to the MPI call address.
        parent = line.strip('[]')

        #If key already exists we take a copy of the key's current value list.
        if tag in dict_tag_fn_map:
            list_callers = dict_tag_fn_map[tag]
        else:
            list_callers = []

        #Add an address only if it does not already appear in the list.
        if (caller,parent) not in list_callers:
            list_callers.append( (caller,parent) )

        dict_tag_fn_map[tag] = list_callers


        count = 0 #Reset our counter.


#At this point I have the fn addresses in a dictionary of list.
#Translate each fn address to a name and print out.
pwdStr = commands.getoutput('pwd')
list_tags_with_list_fns = dict_tag_fn_map.items()
list_tags_with_list_fns.sort()

for (tag, fns) in list_tags_with_list_fns:
    for (fn1,fn2) in fns:
        fnName1 = get_stack_trace(fn1, pwdStr)
        fnName2 = get_stack_trace(fn2, pwdStr)
        print "tag:", tag, "Calling fn:", fnName1, "  Parent fn:", fnName2

fileHandle.close()



OUTPUT:

tag: 000100 Calling fn: mpi_lib.F90:859   Parent fn: mpi_amr_comm_setup.F90:528
tag: 000100 Calling fn: mpi_lib.F90:891   Parent fn: mpi_amr_comm_setup.F90:528
tag: 000120 Calling fn: process_fetch_list.F90:166   Parent fn: mpi_morton_bnd.F90:381
tag: 000120 Calling fn: process_fetch_list.F90:188   Parent fn: mpi_morton_bnd.F90:381
tag: 000130 Calling fn: process_fetch_list.F90:166   Parent fn: mpi_morton_bnd_prolong.F90:384
tag: 000130 Calling fn: process_fetch_list.F90:188   Parent fn: mpi_morton_bnd_prolong.F90:384
tag: 000140 Calling fn: process_fetch_list.F90:166   Parent fn: mpi_morton_bnd_fluxcon.F90:335
tag: 000140 Calling fn: process_fetch_list.F90:188   Parent fn: mpi_morton_bnd_fluxcon.F90:335
tag: 000150 Calling fn: process_fetch_list.F90:188   Parent fn: mpi_morton_bnd_restrict.F90:168



For bulk processing of stack trace files without MPI tags (Here NFRAMES = 3):
addr2line $(cat $(for i in $(ls stacktrace.out.*); do echo $i; done) | 
sed -n '3~3p' | tr -d [ | tr -d ]) -e flash3 | sort | uniq

*/
