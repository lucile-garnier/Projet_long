PROGRAM gmsh2cluster
  IMPLICIT NONE
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  CHARACTER(30) :: file_name
  INTEGER :: dim
  INTEGER :: i
  INTEGER :: j
  INTEGER :: nb_nodes
  DOUBLE PRECISION :: x
  DOUBLE PRECISION :: y
  DOUBLE PRECISION :: z

  !###########################################
  ! INSTRUCTIONS
  !###########################################
  PRINT *
  PRINT *, 'Conversion from gmsh to input cluster'
  PRINT *
  PRINT *, 'Name of gmsh file?'
  READ *, file_name
  OPEN(FILE=file_name,UNIT=1)
  PRINT *, 'Name of output file?'
  READ *, file_name
  OPEN(FILE=file_name,UNIT=2)
  PRINT *, 'DIMENSION ?'
  READ *, dim
  READ (1,*)
  READ (1,*)  
  READ (1,*)
  READ (1,*)  
  READ (1,*)  nb_nodes
  PRINT *, 'Number of nodes :', nb_nodes
  WRITE(2,*) nb_nodes, dim
  DO i=1,nb_nodes
     READ (1,*) j,x,y,z
     IF (dim==2) THEN
        WRITE(2,*) x,y
     ELSE
        WRITE (2,*) x,y,z
     ENDIF
  ENDDO
  CLOSE(1)
  CLOSE(2)
  STOP
END PROGRAM gmsh2cluster
