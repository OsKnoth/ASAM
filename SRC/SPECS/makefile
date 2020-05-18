include ../$(IN)
#
# Erstellt von Harald Heinrich 31.03.2005
# geaendert von Verena Gruetzun 06.12.2005
#                   Angelegt wird die Bibliothek liblm.a
#
#
SCRIPTS = *.f90
OBJECTS = *.o
#
FC =     mpxlf90 -qsuffix=f=f90:cpp=f90
FFLAGS = -c -O3 -q64 -qarch=auto -qtune=auto -qstrict -qspillsize=32704 -bmaxdata:0x80000000 -qmaxmem=-1 -qinitauto
#
MIKROXD.a : $(SCRIPTS)
	$(FC) $(FFLAGS) $(KINCL) $?
	ar -X64 rsv libmp.a $(?:.f90=.o)
#find . -name '*.o' -exec rm -i '{}'\;
