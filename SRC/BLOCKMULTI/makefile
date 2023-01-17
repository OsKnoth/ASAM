FC=                                     mpif90

FFLAGS=                                 -c -fallow-argument-mismatch -O3
LN=                                     -L/usr/local/lib -lnetcdf -lnetcdff
LB=                                     -llapack -lblas

OBJS=                                   test.o Parallel_Mod.o Kind_Mod_Intel.o MatSpRowCol_Mod.o Index_Mod.o \
                                        Subdom_Mod.o DomDecomp_Mod.o Block_Mod.o Part_Mod.o \
                                        CoarseGrids_Mod.o List_Mod.o Sort_Mod.o Transfer_Mod.o \
                                        Boundary_Mod.o Laplace_block_Mod.o Gradient_block_Mod.o \
                                        Gradient_block_interpolation_Mod.o Crop_Mod.o Smooth_Mod.o \
                                        Output_Mod.o Compressible_Mod.o Constants_Mod.o

boundary_mod_dpds =                     Boundary_Mod.f90 kind_mod.mod block_mod.mod \
                                        index_mod.mod matsprowcol_mod.mod parallel_mod.mod

block_mod_dpds =                        Block_Mod.f90 kind_mod.mod

coarsegrids_mod_dpds =                  CoarseGrids_Mod.f90 subdom_mod.mod domdecomp_mod.mod \
                                        matsprowcol_mod.mod crop_mod.mod block_mod.mod \
                                        list_mod.mod sort_mod.mod transfer_mod.mod index_mod.mod \
                                        kind_mod.mod parallel_mod.mod laplace_block_mod.mod \
                                        gradient_block_interpolation_mod.mod
constants_mod_dpds = 			Constants_Mod.f90 kind_mod.mod

compressible_mod_dpds =                 Compressible_Mod.f90 subdom_mod.mod matsprowcol_mod.mod \
                                        kind_mod.mod constants_mod.mod parallel_mod.mod

crop_mod_dpds =                         Crop_Mod.f90 kind_mod.mod block_mod.mod matsprowcol_mod.mod \
                                        subdom_mod.mod

domdecomp_mod_dpds =                    DomDecomp_Mod.f90 index_mod.mod block_mod.mod \
                                        boundary_mod.mod subdom_mod.mod matsprowcol_mod.mod part_mod.mod

gradient_block_mod_dpds =               Gradient_block_Mod.f90 index_mod.mod kind_mod.mod \
                                        matsprowcol_mod.mod block_mod.mod gradient_block_interpolation_mod.mod

gradient_block_interpolation_mod_dpds = Gradient_block_interpolation_Mod.f90 kind_mod.mod \
                                        block_mod.mod

index_mod_dpds =                        Index_Mod.f90 kind_mod.mod block_mod.mod

kind_mod_dpds =                         Kind_Mod_Intel.f90

laplace_block_mod_dpds =                Laplace_block_Mod.f90 kind_mod.mod index_mod.mod \
                                        block_mod.mod crop_mod.mod subdom_mod.mod \
                                        gradient_block_mod.mod gradient_block_interpolation_mod.mod \
                                        crop_mod.mod smooth_mod.mod

list_mod_dpds =                         List_Mod.f90 kind_mod.mod sort_mod.mod parallel_mod.mod

matsprowcol_mod_dpds =                  MatSpRowCol_Mod.f90 kind_mod.mod index_mod.mod

parallel_mod_dpds =                     Parallel_Mod.f90 kind_mod.mod

output_mod_dpds =                       Output_Mod.f90 kind_mod.mod domdecomp_mod.mod

part_mod_dpds =                         Part_Mod.f90 kind_mod.mod parallel_mod.mod

subdom_mod_dpds =                       Subdom_Mod.f90 kind_mod.mod index_mod.mod block_mod.mod \
                                        boundary_mod.mod matsprowcol_mod.mod kind_mod.mod smooth_mod.mod

smooth_mod_dpds =                       Smooth_Mod.f90 matsprowcol_mod.mod kind_mod.mod boundary_mod.mod \
                                        Block_Mod.f90

sort_mod_dpds =                         Sort_Mod.f90 kind_mod.mod 

test_dpds =                             test.f90 parallel_mod.mod kind_mod.mod matsprowcol_mod.mod index_mod.mod \
                                        subdom_mod.mod domdecomp_mod.mod block_mod.mod coarsegrids_mod.mod \
                                        list_mod.mod laplace_block_mod.mod smooth_mod.mod output_mod.mod \
                                        compressible_mod.mod

transfer_mod_dpds =                     Transfer_Mod.f90 index_mod.mod kind_mod.mod matsprowcol_mod.mod \
                                        block_mod.mod boundary_mod.mod parallel_mod.mod boundary_mod.mod \
                                        list_mod.mod 

test_run: $(OBJS)
	$(FC) -o $@ $^ $(LB) $(LN)
Block_Mod.o: $(block_mod_dpds)
	$(FC) $(FFLAGS) $<
block_mod.mod: $(block_mod_dpds)
	$(FC) $(FFLAGS) $<
Boundary_Mod.o: $(boundary_mod_dpds)
	$(FC) $(FFLAGS) $<
boundary_mod.mod: $(boundary_mod_dpds)
	$(FC) $(FFLAGS) $<
coarsegrids_mod.mod: $(coarsegrids_mod_dpds)
	$(FC) $(FFLAGS) $<
CoarseGrids_Mod.o: $(coarsegrids_mod_dpds)
	$(FC) $(FFLAGS) $<
Constants_Mod.o: $(constants_mod_dpds)
	$(FC) $(FFLAGS) $<
constants_mod.mod: $(constants_mod_dpds)
	$(FC) $(FFLAGS) $<
compressible_mod.mod: $(compressible_mod_dpds)
	$(FC) $(FFLAGS) $<
Compressible_Mod.o: $(compressible_mod_dpds)
	$(FC) $(FFLAGS) $<
crop_mod.mod: $(crop_mod_dpds)
	$(FC) $(FFLAGS) $<
Crop_Mod.o: $(crop_mod_dpds)
	$(FC) $(FFLAGS) $<
DomDecomp_Mod.o: $(domdecomp_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
domdecomp_mod.mod: $(domdecomp_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
gradient_block_mod.mod: $(gradient_block_mod_dpds)
	$(FC) $(FFLAGS) $<
Gradient_block_Mod.o: $(gradient_block_mod_dpds)
	$(FC) $(FFLAGS) $<
Gradient_block_interpolation_Mod.o: $(gradient_block_interpolation_mod_dpds)
	$(FC) $(FFLAGS) $< $(LB)
gradient_block_interpolation_mod.mod: $(gradient_block_interpolation_mod_dpds)
	$(FC) $(FFLAGS) $< $(LB)
Index_Mod.o: $(index_mod_dpds)
	$(FC) $(FFLAGS) $<
index_mod.mod: $(index_mod_dpds)
	$(FC) $(FFLAGS) $<
Kind_Mod_Intel.o: $(kind_mod_dpds)
	$(FC) $(FFLAGS) $<
kind_mod.mod: $(kind_mod_dpds)
	$(FC) $(FFLAGS) $<
Laplace_block_Mod.o: $(laplace_block_mod_dpds)
	$(FC) $(FFLAGS) $<
laplace_block_mod.mod: $(laplace_block_mod_dpds)
	$(FC) $(FFLAGS) $<
list_mod.mod: $(list_mod_dpds)
	$(FC) $(FFLAGS) $<
List_Mod.o: $(list_mod_dpds)
	$(FC) $(FFLAGS) $<
MatSpRowCol_Mod.o: $(matsprowcol_mod_dpds)
	$(FC) $(FFLAGS) $<
matsprowcol_mod.mod:  $(matsprowcol_mod_dpds)
	$(FC) $(FFLAGS) $<
parallel_mod.mod: $(parallel_mod_dpds)
	$(FC) $(FFLAGS) $<
Parallel_Mod.o: $(parallel_mod_dpds)
	$(FC) $(FFLAGS) $<
Output_Mod.o: $(output_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
output_mod.mod: $(output_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
Part_Mod.o: $(part_mod_dpds)
	$(FC) $(FFLAGS) $<
part_mod.mod: $(part_mod_dpds)
	$(FC) $(FFLAGS) $<
Subdom_Mod.o: $(subdom_mod_dpds)
	$(FC) $(FFLAGS) $<
subdom_mod.mod: $(subdom_mod_dpds)
	$(FC) $(FFLAGS) $<
Smooth_Mod.o: $(smooth_mod_dpds)
	$(FC) $(FFLAGS) $<
smooth_mod.mod: $(smooth_mod_dpds)
	$(FC) $(FFLAGS) $<
sort_mod.mod: $(sort_mod_dpds)
	$(FC) $(FFLAGS) $<
Sort_Mod.o: $(sort_mod_dpds)
	$(FC) $(FFLAGS) $<
test.o: $(test_dpds)
	$(FC) $(FFLAGS) $<
Transfer_Mod.o: $(transfer_mod_dpds)
	$(FC) $(FFLAGS) $<
transfer_mod.mod: $(transfer_mod_dpds)
	$(FC) $(FFLAGS) $<

.PHONY : clean

clean:
	rm *.mod
	rm *.o
	rm *.smod
