
#include Makefile.${ARCH}
include Makefile.hydre

#Modules
	MODULES=module_structure.o module_entree.o module_decoupe.o \
		module_MPI.o module_calcul.o module_solve.o module_embed.o \
		 module_sortie.o 
	#module_sparse.o

all :  visudecoup visuclusters clusters gmsh2cluster

clusters : $(MODULES) clusters.o
	$(MPIF90) $(MODULES) clusters.o  -o clusters $(LIBLAPACK) $(LIBPATH) $(LIBARP)
	mv clusters $(BINDIR)

clean :
	rm -vf $(MODULES) *.o clusters visudecoup visuclusters gmsh2cluster
	rm -vf *.mod
	rm -vf $(BINDIR)/clusters $(BINDIR)/visuclusters

clusters.o : $(MODULES) clusters.f90
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c clusters.f90

module_structure.o : module_structure.f90
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c module_structure.f90

module_solve.o : module_solve.f90
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(INCARP) $(LIBPATH) -c module_solve.f90

module_embed.o : module_structure.o module_embed.f90
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c module_embed.f90

module_entree.o : module_structure.o module_entree.f90
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c module_entree.f90

module_decoupe.o : module_structure.o module_decoupe.f90  module_sortie.o
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c module_decoupe.f90

module_MPI.o : module_structure.o module_MPI.f90
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c module_MPI.f90

module_sortie.o : module_structure.o module_sortie.f90
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c module_sortie.f90

module_calcul.o : module_structure.o module_calcul.f90 module_solve.o \
	 module_embed.o
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(LIBPATH) -c module_calcul.f90

module_sparse.o : module_structure.o module_sparse.f90 module_solve.o \
	 module_embed.o
	$(MPIF90) $(FLAGS) $(OPTIONS) $(INCLUDE) $(INCARP) $(LIBPATH) -c module_sparse.f90
#****************************************************************************
# convertisseur gmsh > format d'entree de clusters
#****************************************************************************

gmsh2cluster : gmsh2cluster.f90
	$(F90) $(OPTIONS) gmsh2cluster.f90 -o gmsh2cluster
	@echo
	@echo " *** Programme GMSH2CLUSTER compile ! ***"
	@echo

#****************************************************************************
# outils de visualisation
#****************************************************************************

visudecoup : visudecoup.f90
	$(F90) $(OPTIONS) visudecoup.f90 -o visudecoup
	@echo
	@echo " *** Programme VISUDECOUP compile ! ***"
	@echo

visuclusters : visuclusters.f90 module_visuclusters.o \
	module_visuclusters_structure.o
	$(F90) $(OPTIONS) visuclusters.f90 module_visuclusters.o \
	module_visuclusters_gmsh.o module_visuclusters_structure.o \
	module_visuclusters_paraview.o -o visuclusters
	mv visuclusters $(BINDIR)
	@echo
	@echo " *** Programme VISUCLUSTERS compile ! ***"
	@echo

module_visuclusters.o : module_visuclusters.f90 module_visuclusters_gmsh.o \
	module_visuclusters_structure.o module_visuclusters_paraview.o
	$(F90) $(OPTIONS) -c module_visuclusters.f90

module_visuclusters_gmsh.o : module_visuclusters_gmsh.f90 \
	module_visuclusters_structure.o
	$(F90) $(OPTIONS) -c module_visuclusters_gmsh.f90

module_visuclusters_paraview.o : module_visuclusters_paraview.f90 \
	module_visuclusters_structure.o
	$(F90) $(OPTIONS) -c module_visuclusters_paraview.f90

module_visuclusters_structure.o : module_visuclusters_structure.f90
	$(F90) $(OPTIONS) -c module_visuclusters_structure.f90


#****************************************************************************
# outils de test automatique
#****************************************************************************

teste_clusters : teste_clusters.f90 module_teste_clusters.o
	$(F90) $(OPTIONS) teste_clusters.f90 module_teste_clusters.o \
	 -o teste_clusters
	@echo
	@echo " *** Programme TESTE_CLUSTERS compile ! ***"
	@echo

module_teste_clusters.o : module_teste_clusters.f90
	$(F90) $(OPTIONS) -c module_teste_clusters.f90
