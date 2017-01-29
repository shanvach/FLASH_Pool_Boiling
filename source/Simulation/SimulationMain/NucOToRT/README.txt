This setup generates a utility for mapping the nucleosynthesis output from
trajectory data back onto a grid, and then generating a file that us meant
as input for radiation tranfer processiong by the likes of Phoenix and Sedona.

To use this, there is a pre-configuration step to generate an appropriate
Config file, which will then by used by the FLASH setup script as usual to
build the application. The pre-configuration step has to be done manually
and needs on of the nucleosynthesis output files as input, because it
reads the number and names of isotopes (and some other particle properties)
for there.  This way, the utility adapts semi-automatically to changes in
the number of isotopes (and, just maybe, to other minor changes).

How to preconfigure:

1.     cd .../SimulationMain/NucOToRT

2.     cp Config.top Config

3.     bash CustConfig NUCOUTPUT_FILE >> Config

     (where NUCOUTPUT_FILE stand for name and location of a nucleosynthesis
     output file)

4.   Review the Config file, hopefully it makes sense.

5.   The LAST non-empty line in Config should looks something like
     this:

       PARAMETER particle_attribute_250 STRING "he4"

     Take note of the number NNN in the particle_attribute_NNN runtime
     parameter name, you will need it in the next step.

6.   Run setup (after switching back to the appropriate directory), with some
     additional command line options.  I use something like this:

       ./setup -auto NucOToRT -objdir=NucOToRT -noc -debug -2d +curvilinear +pm4dev --unit=IO/IOParticles +parallelIO -defines=PT_MAX_ATTRIBUTES=NNN Convolve=True -nxb=64 -nyb=128 -parfile=Set-2D-2-noconvo.par -maxblocks=1 -gridinterpolation=native

     where NNN is the number from the previous step.  Modify as needed.


7.   change into the object directory and run make.  This should generate
     flash3, as usual.

8.   Modify one of the sample runtime parameter files to get an appropriate
     flash.par file for your purposes.
8a.  You may want to modify the runtime parameters tinitial and tmax (both)
     to the appropriate time.  This isn't important for correct processing
     (as long as they are the same), but it's the time that will show up
     on diagnostic plots (see 10a).



9.   Make sure the appropriate input files are in the run directory.

10.  Run flash3.

10a. You'll get a number of checkpoint and plot files showing various stages
     of processing. Note that fidlr (i.e., xflash3) chokes on the checkpoint
     files but should be able to display the plot files if the number of
     variables chosen for plotting is not too large.  Visit should also
     be able to display the files.
     Note that checkpoint files will not contain all the particles if they
     were read in from several files (sim_ptNumPartFiles > 1), omly the
     set from the last file read is retained. But UNK data will be shown
     as accumulated from all input files.

11.  If all goes well, enjoy your output file named "PhoenixInputData".

     Note that currently, the output file still contains 1 additional layer of
     output cells in each direction; cells of the PIKN FLASH grid that do not
     fall into any of the inner output grid cells contribute to these pseudo-
     cells.  To get rid of then, you can use something like

       awk '{if ($1 != 0 && $2 != 0 && $1 != 65 && $2 != 129) print;}'  PhoenixInputData > PhoenixInputData.inner

     (Adjust the numbers for upper indices appropriately if your output
     grid is not sized 64 x 128.)


CAVEATS
=======
Things aren't quite as automatic and as configurable as maybe they should.
In particular, the list of variables to write to the main output file is
currently hardwired in sim_writeOutputGrid.F90.  And if the number of
variables changes - evne if NUNK_VARS changes, without changes to the
list of columns to be output - the dimensioning in following lines has to
be changed MANUALLY to be consistent:

sim_initOutputGrid.F90:
  nvarsOgOut = NUNK_VARS - 12           !for now - no velz,gaus,gpot,RPV{1,2,3},gam{c,e},eint,grac,_numc_,_nup{0,1}_,pden,entr

sim_writeOutputGrid.F90:
  integer :: outVars(NUNK_VARS-12)



