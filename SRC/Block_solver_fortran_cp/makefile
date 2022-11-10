FC=                                     mpif90

FFLAGS=                                 -c -fallow-argument-mismatch -O3
LN=                                     -L/usr/local/lib -lnetcdf -lnetcdff
LB=                                     -llapack -lblas

OBJS=                                   test.o mpi_mod.o kind_mod_intel.o matsprowcol_mod.o index_mod.o \
                                        subdom_mod.o domdecomp_mod.o block_mod.o part_man_mod.o \
                                        coarsegrids_mod.o list_mod.o sort_mod.o transfer_mod.o \
                                        boundary_mod.o laplace_block_mod.o gradient_block_mod.o \
                                        gradient_block_interpolation_mod.o crop_mod.o smooth_mod.o \
                                        output_mod.o compressible_mod.o

boundary_mod_dpds =                     Boundary_Mod.f90 kind_mod.mod block_mod.mod \
                                        index_mod.mod matsprowcol_mod.mod mpi_mod.mod

block_mod_dpds =                        Block_Mod.f90 kind_mod.mod

coarsegrids_mod_dpds =                  CoarseGrids_Mod.f90 subdom_mod.mod domdecomp_mod.mod \
                                        matsprowcol_mod.mod crop_mod.mod block_mod.mod \
                                        list_mod.mod sort_mod.mod transfer_mod.mod index_mod.mod \
                                        kind_mod.mod mpi_mod.mod laplace_block_mod.mod \
                                        gradient_block_interpolation_mod.mod

compressible_mod_dpds =                 Compressible_Mod.f90 subdom_mod.mod matsprowcol_mod.mod \
                                        kind_mod.mod

crop_mod_dpds =                         Crop_Mod.f90 kind_mod.mod block_mod.mod matsprowcol_mod.mod \
                                        subdom_mod.mod

domdecomp_mod_dpds =                    DomDecomp_Mod.f90 index_mod.mod block_mod.mod \
                                        boundary_mod.mod subdom_mod.mod matsprowcol_mod.mod

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

list_mod_dpds =                         List_Mod.f90 kind_mod.mod sort_mod.mod

matsprowcol_mod_dpds =                  MatSpRowCol_Mod.f90 kind_mod.mod index_mod.mod

mpi_mod_dpds =                          MPI_Mod.f90 kind_mod.mod

output_mod_dpds =                       Output_Mod.f90 kind_mod.mod domdecomp_mod.mod

part_man_mod_dpds =                     Part_Man_Mod.f90 kind_mod.mod

subdom_mod_dpds =                       Subdom_Mod.f90 kind_mod.mod index_mod.mod block_mod.mod \
                                        boundary_mod.mod matsprowcol_mod.mod kind_mod.mod smooth_mod.mod

smooth_mod_dpds =                       Smooth_Mod.f90 matsprowcol_mod.mod kind_mod.mod boundary_mod.mod \
                                        Block_Mod.f90

sort_mod_dpds =                         Sort_Mod.f90 kind_mod.mod 

test_dpds =                             test.f90 mpi_mod.mod kind_mod.mod matsprowcol_mod.mod index_mod.mod \
                                        subdom_mod.mod domdecomp_mod.mod block_mod.mod part_man_mod.mod \
                                        coarsegrids_mod.mod list_mod.mod laplace_block_mod.mod smooth_mod.mod \
                                        output_mod.mod compressible_mod.mod

transfer_mod_dpds =                     Transfer_Mod.f90 index_mod.mod kind_mod.mod matsprowcol_mod.mod \
                                        block_mod.mod boundary_mod.mod mpi_mod.mod boundary_mod.mod \
                                        list_mod.mod 

test_run: $(OBJS)
	$(FC) -o $@ $^ $(LB) $(LN)
block_mod.o: $(block_mod_dpds)
	$(FC) $(FFLAGS) $<
block_mod.mod: $(block_mod_dpds)
	$(FC) $(FFLAGS) $<
boundary_mod.o: $(boundary_mod_dpds)
	$(FC) $(FFLAGS) $<
boundary_mod.mod: $(boundary_mod_dpds)
	$(FC) $(FFLAGS) $<
coarsegrids_mod.mod: $(coarsegrids_mod_dpds)
	$(FC) $(FFLAGS) $<
coarsegrids_mod.o: $(coarsegrids_mod_dpds)
	$(FC) $(FFLAGS) $<
compressible_mod.mod: $(compressible_mod_dpds)
	$(FC) $(FFLAGS) $<
compressible_mod.o: $(compressible_mod_dpds)
	$(FC) $(FFLAGS) $<
crop_mod.mod: $(crop_mod_dpds)
	$(FC) $(FFLAGS) $<
crop_mod.o: $(crop_mod_dpds)
	$(FC) $(FFLAGS) $<
domdecomp_mod.o: $(domdecomp_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
domdecomp_mod.mod: $(domdecomp_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
gradient_block_mod.mod: $(gradient_block_mod_dpds)
	$(FC) $(FFLAGS) $<
gradient_block_mod.o: $(gradient_block_mod_dpds)
	$(FC) $(FFLAGS) $<
gradient_block_interpolation_mod.o: $(gradient_block_interpolation_mod_dpds)
	$(FC) $(FFLAGS) $< $(LB)
gradient_block_interpolation_mod.mod: $(gradient_block_interpolation_mod_dpds)
	$(FC) $(FFLAGS) $< $(LB)
index_mod.o: $(index_mod_dpds)
	$(FC) $(FFLAGS) $<
index_mod.mod: $(index_mod_dpds)
	$(FC) $(FFLAGS) $<
kind_mod_intel.o: $(kind_mod_dpds)
	$(FC) $(FFLAGS) $<
kind_mod.mod: $(kind_mod_dpds)
	$(FC) $(FFLAGS) $<
laplace_block_mod.o: $(laplace_block_mod_dpds)
	$(FC) $(FFLAGS) $<
laplace_block_mod.mod: $(laplace_block_mod_dpds)
	$(FC) $(FFLAGS) $<
list_mod.mod: $(list_mod_dpds)
	$(FC) $(FFLAGS) $<
list_mod.o: $(list_mod_dpds)
	$(FC) $(FFLAGS) $<
matsprowcol_mod.o: $(matsprowcol_mod_dpds)
	$(FC) $(FFLAGS) $<
matsprowcol_mod.mod:  $(matsprowcol_mod_dpds)
	$(FC) $(FFLAGS) $<
mpi_mod.mod: $(mpi_mod_dpds)
	$(FC) $(FFLAGS) $<
mpi_mod.o: $(mpi_mod_dpds)
	$(FC) $(FFLAGS) $<
output_mod.o: $(output_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
output_mod.mod: $(output_mod_dpds)
	$(FC) $(FFLAGS) $< `nf-config --fflags --flibs`
part_man_mod.o: $(part_man_mod_dpds)
	$(FC) $(FFLAGS) $<
part_man_mod.mod: $(part_man_mod_dpds)
	$(FC) $(FFLAGS) $<
subdom_mod.o: $(subdom_mod_dpds)
	$(FC) $(FFLAGS) $<
subdom_mod.mod: $(subdom_mod_dpds)
	echo compile mod
	$(FC) $(FFLAGS) $<
smooth_mod.o: $(smooth_mod_dpds)
	$(FC) $(FFLAGS) $<
smooth_mod.mod: $(smooth_mod_dpds)
	$(FC) $(FFLAGS) $<
sort_mod.mod: $(sort_mod_dpds)
	$(FC) $(FFLAGS) $<
sort_mod.o: $(sort_mod_dpds)
	$(FC) $(FFLAGS) $<
test.o: $(test_dpds)
	$(FC) $(FFLAGS) $<
transfer_mod.o: $(transfer_mod_dpds)
	$(FC) $(FFLAGS) $<
transfer_mod.mod: $(transfer_mod_dpds)
	$(FC) $(FFLAGS) $<

.PHONY : clean

clean:
	rm *.mod
	rm *.o
	rm *.smod
