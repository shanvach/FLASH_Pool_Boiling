./setup unitTest/XYZpoisson_SuperLU -auto -3d -debug -nxb=8 -nyb=8 -nzb=8 -maxblocks=400 -gridinterpolation=native  -parfile=/Users/mvanella/Documents/ADAPTIVE/FLASH/test/INS_IB/source/Simulation/SimulationMain/unitTest/XYZpoisson_SuperLU/flash_pm_3d.par  +pm4dev  -objdir=XYZ_POISSON_SUPERLU -site=SEAS10927.gwu.edu -noclobber 


#./setup unitTest/PFFT_XYdir_2D_MG  -auto -2d -debug -nxb=8 -nyb=8 -maxblocks=4000 -parfile=/Users/mvanella/Documents/ADAPTIVE/FLASH/FLASH3_TRUNK/source/Simulation/SimulationMain/unitTest/PFFT_XYdir_2D_MG/flash_pm_2d.par  -objdir=XYDIR_MULTI_2D_MC -site=SEAS10927.gwu.edu -noclobber +noio -gridinterpolation=native Grid=PM4DEV ParameshLibraryMode=True PfftSolver=Generic_Direct

#Grid=PM4DEV ParameshLibraryMode=True
#-gridinterpolation=native
#PfftSolver=Generic_Direct


#./setup unitTest/PFFT_XYdir_2D_MG +ug -auto -2d -debug -nxb=64 -nyb=64 -maxblocks=4000 -parfile=/Users/mvanella/Documents/ADAPTIVE/FLASH/FLASH3_TRUNK/source/Simulation/SimulationMain/unitTest/PFFT_XYdir_2D_MG/flash_pm_2d.par  -objdir=XYDIR_MULTI_2D_PFFT -site=SEAS10927.gwu.edu -noclobber +noio PfftSolver=Generic_Direct

