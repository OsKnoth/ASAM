###############
### Linux

F90 = gfortran -g
#F90 = /opt/irun/bin/gfortran -g 
#F90 = ifort -g -traceback
#F90  = ifort -assume byterecl      # wenn SizeOfReal=4 gesetzt
F90FLAGS = 
# export LC_ALL=C  zuvor setzen


#----------- end of user configuration parameters ------------


READ: Kind_Mod.f90 Geometry_Mod.f90 Shape_Mod.f90 DBase_Mod.f90   OSD_Mod.f90 Main.f90   
	$(F90) $(F90FLAGS) -o READ $(MODS) Kind_Mod.f90 Geometry_Mod.f90 Shape_Mod.f90 DBase_Mod.f90   OSD_Mod.f90 Main.f90

COMPOSE: Kind_Mod.f90 Geometry_Mod.f90 Shape_Mod.f90 DBase_Mod.f90   OSD_Mod.f90 MainCompose.f90   
	$(F90) $(F90FLAGS) -o COMPOSE $(MODS) Kind_Mod.f90 Geometry_Mod.f90 Shape_Mod.f90 DBase_Mod.f90   OSD_Mod.f90 MainCompose.f90
