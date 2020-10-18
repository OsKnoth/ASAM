################################################################################
###                                                                          ###
###                          M A K E F I L E                                 ###
###                                                                          ###
###              for the All Scale Atmospheric Model (ASAM)                  ###
###                                                                          ###
###      Institute for Tropospheric Research (TROPOS), Leipzig, Germany      ###
###                                                                          ###
################################################################################


IN = M_DEF_$(MACH)
include ${IN}

GRIDLIB = -L$(LIBGRID) -lgrid

SHAPELIB = -L$(LIBSHAPE) -lshape

CHEMSYSLIB = -L$(LIBCHEMSYS) -lchemsys

MODLIB = -L$(LIB) -lmod

CGLIB = -L$(LIB) -lcg

WINDFARMLIB =-L$(LIB) -lwind

SOILLIB = -L$(LIB) -lsoil

CANOPYLIB = -L$(LIB) -lcanopy

ADVECTLIB = -L$(LIB) -ladvect

ADJUSTLIB = -L$(LIB) -ladjust

INITLIB = -L$(LIB) -linit

MULTILIB = -L$(LIB) -lmulti

MATRIXLIB = -L$(LIB) -lmatrix

SOLULIB = -L$(LIB) -lsolu

PITZERLIB = -L$(LIB) -lpitzer

CHEMIELIB = -L$(LIB) -lchemie

EXAMPLELIB = -L$(LIB) -lexample

RADIATIONLIB = -L$(LIB) -lradiation

MICROPHYSLIB = -L$(LIB) -lmicrophys

SPECSLIB = -L$(LIB) -lspecs

OBJ =  $(MODLIB) $(WINDFARMLIB) $(SOILLIB) $(CANOPYLIB) $(EXAMPLELIB) $(PITZERLIB) $(ADVECTLIB) $(CGLIB) $(ADJUSTLIB) $(INITLIB) \
       $(MULTILIB) $(MATRIXLIB) $(SOLULIB) $(CHEMIELIB) \
       $(RADIATIONLIB) $(MICROPHYSLIB) $(SPECSLIB)

OBJGRID = $(GRIDLIB)       

OBJSHAPE = $(SHAPELIB)       

OBJCHEMSYS = $(CHEMSYSLIB)       
       

Prog: CgD ModD WindFarmD ChemieD SoilD CanopyD PitzerD MicroPhysD SpecsD AdjustD InitD MultiD MatrixD SoluD RadiationD ExampleD AdvectD   
	$(LINK) $(OPT) $(INCL) $(KINCL) -o $(EX) \
        ADVECT/MainProg.f90 EXAMPLE/$(EX).f90 $(OP) $(OBJ) $(PARMETIS) $(MUMPS) $(SCALAP) $(NETCDF) $(LMPI) $(CL)  ; \
        mv $(EX) $(COMPPATH)/$(EX)/$(EX).$(MACH); \
        rm -f ~/bin/$(EX).$(MACH);\
        ln -s $(COMPPATH)/$(EX)/$(EX).$(MACH)  ~/bin/$(EX).$(MACH);\
        rm *.mod  
Read: CgD ModD WindFarmD ChemieD SoilD CanopyD PitzerD MicroPhysD SpecsD AdjustD InitD MultiD MatrixD SoluD RadiationD ExampleD AdvectD   
	$(LINK) $(OPT) $(INCL) $(KINCL) -o Read \
        ADVECT/MainRead.f90 $(OP) $(OBJ) $(PARMETIS) $(NETCDF) $(LMPI) $(CL)  ; \
        mv Read ../bin/Read.$(MACH); \
        rm -f ~/bin/Read.$(MACH);\
        ln -s $(DEFPATH)/bin/Read.$(MACH) ~/bin/Read.$(MACH);
ReadProf: CgD ModD WindFarmD ChemieD SoilD CanopyD PitzerD MicroPhysD SpecsD AdjustD InitD MultiD MatrixD SoluD RadiationD ExampleD AdvectD   
	$(LINK) $(OPT) $(INCL) $(KINCL) -o $(EX) \
        ADVECT/MainReadProf.f90 EXAMPLE/$(EX).f90 $(OP) $(OBJ) $(PARMETIS) $(NETCDF) $(LMPI) $(CL)  ; \
        mv $(EX) $(COMPPATH)/$(EX)/ReadProf.$(MACH); \
        rm *.mod  
TestMPI: CgD ModD WindFarmD ChemieD SoilD CanopyD PitzerD MicroPhysD SpecsD AdjustD InitD MultiD MatrixD SoluD RadiationD ExampleD AdvectD   
	$(LINK) $(OPT) $(INCL) $(KINCL) -o TestMPI \
        ADVECT/MainProgMPI_Test.f90 $(OP) $(OBJ) $(PARMETIS) $(NETCDF) $(LMPI) $(CL)  ; \
        mv TestMPI ~/bin/TestMPI.$(MACH); 
Average: CgD ModD WindFarmD ChemieD SoilD CanopyD PitzerD MicroPhysD SpecsD AdjustD InitD MultiD MatrixD SoluD RadiationD ExampleD AdvectD   
	$(LINK) $(OPT) $(INCL) $(KINCL) -o Average \
        ADVECT/MainTimeAverage.f90 $(OP) $(OBJ) $(PARMETIS) $(NETCDF) $(LMPI) $(CL)  ; \
        mv Average ../bin/Average.$(MACH); \
        rm -f ~/bin/Average.$(MACH);\
        ln -s $(DEFPATH)/bin/Average.$(MACH) ~/bin/Average.$(MACH);\

GridGen: GridD 
	$(LINK) $(OPT) $(INCL) $(KGRIDINCL) -o GridGen\
        GRID/MainProg.f90 $(OP) $(OBJGRID) $(CURFIT) $(SURFIT) $(GK) $(CL) ;\
        mv GridGen $(DEFPATH)/LIBGRID_$(MACH)/Grid.$(MACH);\
        rm -f ~/bin/Grid.$(MACH);\
        ln -s $(DEFPATH)/LIBGRID_$(MACH)/Grid.$(MACH) ~/bin/Grid.$(MACH);

GridCompose:
	$(LINK) -o GridCompose GRID/GridCompose.f90;\
        mv GridCompose ../LIBGRID_$(MACH)/GridCompose.$(MACH);\
        rm -f ~/bin/GridCompose.$(MACH);\
        ln -s $(ASAMPATH)/LIBGRID_$(MACH)/GridCompose.$(MACH) ~/bin/GridCompose.$(MACH);

PreChem: ChemSysD 
	$(LINK) $(OPT) $(INCL) $(KGRIDINCL) -o PreChem\
        $(OP) $(OBJCHEMSYS) ;\
        mv PreChem ../LIBCHEMSYS_$(MACH)/PreChem.$(MACH);\
        rm -f ~/bin/PreChem.$(MACH);\
        ln -s $(ASAMPATH)/LIBCHEMSYS_$(MACH)/PreChem.$(MACH) ~/bin/PreChem.$(MACH);

ShapeRead: ShapeD 
	$(LINK) $(OPT) $(INCL) $(KSHAPEINCL) -o ShapeRead\
        SHAPE/Main.f90 $(OP) $(OBJSHAPE) $(CL) ;\
        mv ShapeRead ../LIBGRID_$(MACH)/ShapeRead.$(MACH);\
        rm -f ~/bin/ShapeRead.$(MACH);\
        ln -s $(DEFPATH)/LIBGRID_$(MACH)/ShapeRead.$(MACH) ~/bin/ShapeRead.$(MACH);

ShapeOSD: ShapeD 
	$(LINK) $(OPT) $(INCL) $(KSHAPEINCL) -o ShapeOSD\
        SHAPE/MainOSD.f90 $(OP) $(OBJSHAPE) ;\
        mv ShapeOSD ../LIBGRID_$(MACH)/ShapeOSD.$(MACH);\
        rm -f ~/bin/ShapeOSD.$(MACH);\
        ln -s $(DEFPATH)/LIBGRID_$(MACH)/ShapeOSD.$(MACH) ~/bin/ShapeOSD.$(MACH);

Adjust: dassl CgD ModD
	$(LINK) $(OPT) -fallow-argument-mismatch $(INCL) $(KINCL) -o Adjust AdjustProfile/Profile_Mod.f90 AdjustProfile/ProfileMain.f90 AdjustProfile/DASSL/*.o\
        $(CGLIB) $(MODLIB);\
        mv Adjust ../bin/AdjustProfile.$(MACH);\
        rm -f ~/bin/AdjustProfile.$(MACH);\
        ln -s $(DEFPATH)/bin/AdjustProfile.$(MACH) ~/bin/AdjustProfile.$(MACH);\

DataWRF: dassl CgD ModD
	$(LINK) $(OPT) $(INCL) $(KINCL) -o DataWRF AdjustProfile/Profile_Mod.f90 AdjustProfile/TestReadWRF.f90 AdjustProfile/ReadWRF.f90 AdjustProfile/DASSL/*.o $(NETCDF)\
        $(CGLIB) $(MODLIB);\
        mv DataWRF ../bin/DataWRF.$(MACH);\
        rm -f ~/bin/DataWRF.$(MACH);\
        ln -s $(DEFPATH)/bin/DataWRF.$(MACH) ~/bin/DataWRF.$(MACH);\

DataWRFD: dassl CgD ModD AdjustD
	$(LINK) $(OPT) $(INCL) $(KINCL) -o DataWRF AdjustProfile/Profile_Mod.f90 INIT/ReadWRFncD.f90 AdjustProfile/DASSL/*.o $(NETCDF)\
        $(CGLIB) $(MODLIB);\
        mv DataWRF ../bin/DataWRFD.$(MACH);\
        rm -f ~/bin/DataWRFD.$(MACH);\
        ln -s $(DEFPATH)/bin/DataWRFD.$(MACH) ~/bin/DataWRFD.$(MACH);\

dassl:
	@cd AdjustProfile/DASSL; make

CgD: 
	@cd CG; make -f Make_src "IN=$(IN)"
ModD: 
	@cd MODULE; make -f Make_src "IN=$(IN)"
WindFarmD:
	@cd WINDFARM; make -f Make_src "IN=$(IN)"
SoilD: 
	@cd SOIL; make -f Make_src "IN=$(IN)"
CanopyD: 
	@cd CANOPY; make -f Make_src "IN=$(IN)"
ChemieD: 
	@cd CHEMIE; make -f Make_src "IN=$(IN)"
PitzerD: 
	@cd PITZER; make -f Make_src "IN=$(IN)"
MicroPhysD:
	@cd MICROPHYS; make -f Make_src "IN=$(IN)"
AdjustD: 
	@cd AdjustProfile; make -f Make_src "IN=$(IN)"
InitD: 
	@cd INIT; make -f Make_src "IN=$(IN)"
AdvectD: 
	@cd ADVECT; make -f Make_src "IN=$(IN)"
MultiD: 
	@cd MULTI; make -f Make_src "IN=$(IN)"
MatrixD: 
	@cd MATRIX; make -f Make_src "IN=$(IN)"
SoluD: 
	@cd SOLU; make -f Make_src "IN=$(IN)"
ExampleD: 
	@cd EXAMPLE; make -f Make_src "IN=$(IN)"
RadiationD:
	@cd STRAHLUNG; make -f Make_src "IN=$(IN)"
SpecsD:
	@cd SPECS; make -f Make_src "IN=$(IN)"
		
ChemSysD: 
	@cd CHEMSYS; make -f Make_src "IN=$(IN)"
GridD: 
	@cd GRID; make -f Make_src "IN=$(IN)"
ShapeD: 
	@cd SHAPE; make -f Make_src "IN=$(IN)"
clean:
	rm ../LIB_Mac/*
