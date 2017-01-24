PROGRAM visudecoup
  IMPLICIT NONE
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  CHARACTER (LEN=30) :: files
  CHARACTER (LEN=30) :: num
  INTEGER :: i
  INTEGER :: j
  INTEGER :: n
  INTEGER :: nb
  DOUBLE PRECISION :: coord(2)
  DOUBLE PRECISION :: x_max
  DOUBLE PRECISION :: x_min
  DOUBLE PRECISION :: y_max
  DOUBLE PRECISION :: y_min

  !###########################################
  ! INSTRUCTIONS
  !###########################################
  PRINT *
  PRINT *,'visualisation of parallele partitionings in 2D'
  PRINT *

  PRINT *,'nb partitionings +1 ? (=nbproc)'
  READ *,n
  PRINT *
  
  ! Geometry of partitioning
  OPEN(FILE='fort.2',UNIT=2)
  OPEN(FILE='decoupe.geo',UNIT=10)
  DO i=1,n-1
     READ(2,*) x_min,y_min,num,x_max,y_max
     WRITE(10,*) 'Point(',4*(i-1)+1,')={',x_min,',',y_min,',0.};'
     WRITE(10,*) 'Point(',4*(i-1)+2,')={',x_max,',',y_min,',0.};'
     WRITE(10,*) 'Point(',4*(i-1)+3,')={',x_max,',',y_max,',0.};'
     WRITE(10,*) 'Point(',4*(i-1)+4,')={',x_min,',',y_max,',0.};'
     WRITE(10,*) 'Line(',4*(i-1)+1,')={',4*(i-1)+1,',',4*(i-1)+2,'};'
     WRITE(10,*) 'Line(',4*(i-1)+2,')={',4*(i-1)+2,',',4*(i-1)+3,'};'
     WRITE(10,*) 'Line(',4*(i-1)+3,')={',4*(i-1)+3,',',4*(i-1)+4,'};'
     WRITE(10,*) 'Line(',4*(i-1)+4,')={',4*(i-1)+4,',',4*(i-1)+1,'};'
  ENDDO
  CLOSE(10)

  ! Output file
  OPEN(FILE='decoupe.visu',UNIT=1)
  WRITE(1,*) 'View "MPI" {'
  ! Reads the files
  DO i=0,n-1
     ! File name
     IF (i<10) THEN
        WRITE(num,'(i1)') i
     ELSEIF (i<100) THEN
        WRITE(num,'(i2)') i
     ELSEIF (i<1000) THEN
        WRITE(num,'(i3)') i
     ENDIF
     files='decoupe.'//trim(adjustl(num))
     OPEN(FILE=files,UNIT=10)
     READ(10,*) nb
     PRINT *,'  > ',i,' :',nb
     DO j=1,nb
        READ(10,*) coord(:)
        WRITE(1,*) 'SP(',coord(1),',',coord(2),',',0.,'){',i,'};'
     ENDDO
     CLOSE(10)

  ENDDO
  WRITE(1,*) '};'
  CLOSE(1)
  PRINT *
  PRINT *,'gmsh decoupe.visu'
  PRINT *,'gmsh decoupe.geo'
  PRINT *,'gmsh decoupe.visu decoupe.geo'
  STOP
END PROGRAM visudecoup
