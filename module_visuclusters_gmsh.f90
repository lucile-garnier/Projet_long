!>Contains methods enabling writing results in file specifically formatted for GMSH
MODULE module_visuclusters_gmsh
  USE module_visuclusters_structure
CONTAINS



!>Writes the geometry partitioning for Gmsh visualization
!!@details This methods extracts the domain definitions
!!from the <em>fort.2</em> file and writes to <em>decoupe.geo</em> file,
!!a source code file. We refer you to a <a href="linkURL"> Gmsh tutorial</a>.
!!@see module_calcul::write_partitioning()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_partitioning_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params

    !#### Variables  ####
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    DOUBLE PRECISION :: x0
    DOUBLE PRECISION :: x1
    DOUBLE PRECISION :: x_max
    DOUBLE PRECISION :: x_min
    DOUBLE PRECISION :: y0
    DOUBLE PRECISION :: y1
    DOUBLE PRECISION :: y_max
    DOUBLE PRECISION :: y_min
    DOUBLE PRECISION :: z0
    DOUBLE PRECISION :: z1
    DOUBLE PRECISION :: zmax
    DOUBLE PRECISION :: zmin

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *,'-> decoupe.geo'
    OPEN(FILE='fort.2',UNIT=2)
    OPEN(FILE='decoupe.geo',UNIT=10)
    DO i=1,params%nb_proc-params%is_interfacing
       IF (((params%coords==1).AND.(params%dim==2)) &
          .OR.((params%is_image==1).AND.(params%image_dim==2)) &
          .OR.((params%is_threshold==1).AND.(params%image_dim==2)) &
          .OR.((params%is_geom==1).AND.(params%image_dim+params%image_times==2))) THEN
          !2D
          IF (params%coords==1) THEN
             READ(2,*) x0,y0,num,x1,y1
             x_min=x0
             y_min=y0
             x_max=x1
             y_max=y1
          ELSEIF (params%is_threshold==1) THEN
             READ(2,*) x0,num,x1
             x_min=1
             y_min=-1
             x_max=params%partitioning(2)
             y_max=-params%partitioning(1)
             zmin=x0
             zmax=x1
          ELSEIF ((params%is_image==1).OR.(params%is_geom==1)) THEN
             READ(2,*) x0,y0,num,x1,y1
             x_min=y0
             y_min=-x0
             x_max=y1
             y_max=-x1
             zmin=x0
          ENDIF
          IF (params%is_threshold==1) THEN
             WRITE(10,*) 'Point(',8*(i-1)+1,')={',x_min,',',y_min,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+2,')={',x_max,',',y_min,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+3,')={',x_max,',',y_max,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+4,')={',x_min,',',y_max,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+5,')={',x_min,',',y_min,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+6,')={',x_max,',',y_min,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+7,')={',x_max,',',y_max,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+8,')={',x_min,',',y_max,',',zmax,'};'
             WRITE(10,*) 'Line(',12*(i-1)+1,')={',8*(i-1)+1,',',8*(i-1)+2,'};'
             WRITE(10,*) 'Line(',12*(i-1)+2,')={',8*(i-1)+2,',',8*(i-1)+3,'};'
             WRITE(10,*) 'Line(',12*(i-1)+3,')={',8*(i-1)+3,',',8*(i-1)+4,'};'
             WRITE(10,*) 'Line(',12*(i-1)+4,')={',8*(i-1)+4,',',8*(i-1)+1,'};'
             WRITE(10,*) 'Line(',12*(i-1)+5,')={',8*(i-1)+5,',',8*(i-1)+6,'};'
             WRITE(10,*) 'Line(',12*(i-1)+6,')={',8*(i-1)+6,',',8*(i-1)+7,'};'
             WRITE(10,*) 'Line(',12*(i-1)+7,')={',8*(i-1)+7,',',8*(i-1)+8,'};'
             WRITE(10,*) 'Line(',12*(i-1)+8,')={',8*(i-1)+8,',',8*(i-1)+5,'};'
             WRITE(10,*) 'Line(',12*(i-1)+9,')={',8*(i-1)+1,',',8*(i-1)+5,'};'
             WRITE(10,*) 'Line(',12*(i-1)+10,')={',8*(i-1)+2,',',8*(i-1)+6,'};'
             WRITE(10,*) 'Line(',12*(i-1)+11,')={',8*(i-1)+3,',',8*(i-1)+7,'};'
             WRITE(10,*) 'Line(',12*(i-1)+12,')={',8*(i-1)+4,',',8*(i-1)+8,'};'
          ELSE
             WRITE(10,*) 'Point(',4*(i-1)+1,')={',x_min,',',y_min,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+2,')={',x_max,',',y_min,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+3,')={',x_max,',',y_max,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+4,')={',x_min,',',y_max,',0.};'
             WRITE(10,*) 'Line(',4*(i-1)+1,')={',4*(i-1)+1,',',4*(i-1)+2,'};'
             WRITE(10,*) 'Line(',4*(i-1)+2,')={',4*(i-1)+2,',',4*(i-1)+3,'};'
             WRITE(10,*) 'Line(',4*(i-1)+3,')={',4*(i-1)+3,',',4*(i-1)+4,'};'
             WRITE(10,*) 'Line(',4*(i-1)+4,')={',4*(i-1)+4,',',4*(i-1)+1,'};'
          ENDIF
       ELSEIF (((params%coords==1).AND.(params%dim==3)) &
          .OR.((params%is_image==1).AND.(params%image_dim==3)) &
          .OR.((params%is_geom==1).AND.(params%image_dim+params%image_times==3))) THEN
          !3D
          IF (params%coords==1) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             x_min=x0
             y_min=y0
             zmin=z0
             x_max=x1
             y_max=y1
             zmax=z1
          ELSEIF ((params%is_image==1).OR.(params%is_geom==1)) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             x_min=y0
             y_min=-x0
             zmin=z0
             x_max=y1
             y_max=-x1
             zmax=z1
          ENDIF
          WRITE(10,*) 'Point(',8*(i-1)+1,')={',x_min,',',y_min,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+2,')={',x_max,',',y_min,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+3,')={',x_max,',',y_max,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+4,')={',x_min,',',y_max,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+5,')={',x_min,',',y_min,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+6,')={',x_max,',',y_min,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+7,')={',x_max,',',y_max,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+8,')={',x_min,',',y_max,',',zmax,'};'
          WRITE(10,*) 'Line(',12*(i-1)+1,')={',8*(i-1)+1,',',8*(i-1)+2,'};'
          WRITE(10,*) 'Line(',12*(i-1)+2,')={',8*(i-1)+2,',',8*(i-1)+3,'};'
          WRITE(10,*) 'Line(',12*(i-1)+3,')={',8*(i-1)+3,',',8*(i-1)+4,'};'
          WRITE(10,*) 'Line(',12*(i-1)+4,')={',8*(i-1)+4,',',8*(i-1)+1,'};'
          WRITE(10,*) 'Line(',12*(i-1)+5,')={',8*(i-1)+5,',',8*(i-1)+6,'};'
          WRITE(10,*) 'Line(',12*(i-1)+6,')={',8*(i-1)+6,',',8*(i-1)+7,'};'
          WRITE(10,*) 'Line(',12*(i-1)+7,')={',8*(i-1)+7,',',8*(i-1)+8,'};'
          WRITE(10,*) 'Line(',12*(i-1)+8,')={',8*(i-1)+8,',',8*(i-1)+5,'};'
          WRITE(10,*) 'Line(',12*(i-1)+9,')={',8*(i-1)+1,',',8*(i-1)+5,'};'
          WRITE(10,*) 'Line(',12*(i-1)+10,')={',8*(i-1)+2,',',8*(i-1)+6,'};'
          WRITE(10,*) 'Line(',12*(i-1)+11,')={',8*(i-1)+3,',',8*(i-1)+7,'};'
          WRITE(10,*) 'Line(',12*(i-1)+12,')={',8*(i-1)+4,',',8*(i-1)+8,'};'
       ENDIF
    ENDDO
    CLOSE(10)
    CLOSE(2)
    RETURN
  END SUBROUTINE write_partitioning_gmsh



!>Initializes the file of the partitioning for Gmsh visualization
!!@details This method extracts details on partitioning from the
!!<em>decoupe.x</em> files and writes them in <em>decoupe.visu</em>
!!file.
!!@note <em>x</em> is the identifier (number) of a process.
!!@see write_partitioning()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_assignment_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: nb
    INTEGER :: nb_slaves
    INTEGER :: offset
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE='decoupe.visu',UNIT=1)
    WRITE(1,*) 'View "MPI" {'
    !Reading files
    IF (params%nb_proc==1) THEN
       offset=1
       nb_slaves=1
    ELSE
       offset=0
       nb_slaves=params%nb_proc-1
    ENDIF
    ALLOCATE(coords(1,params%dim))
    DO i=offset,nb_slaves
       ! File name
       WRITE(num,*) i
       files='decoupe.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       READ(10,*) nb
       PRINT *, '> ', i, ' : ', nb
       DO j=1,nb
          IF (params%coords==1) THEN
             ! Partitioning by coordinates
             coords(:,:)=0.
             READ(10,*) coords(1,:)
             CALL write_point_coord_format(params%dim, i, 1, 1, coords)
          ELSE
             ! Partitioning 1D pictures
             READ (10,*) k
             CALL write_point_picture_format(params, i, k, 1)
          ENDIF
       ENDDO
       CLOSE(10)
    ENDDO
    DEALLOCATE(coords)
    WRITE(1,*) '};'
    CLOSE(1)
    PRINT *, '-> decoupe.visu'
    RETURN
  END SUBROUTINE write_assignment_gmsh



!>Writes the clusters before grouping for Gmsh visualization
!!@details This methods extracts details on computed clusters
!!on each domain from <em>cluster.partiel.x</em> files and writes
!!them in corresponding <em>cluster.partiel.x.visu</em> files.
!!@note <em>x</em> is the identifier (number) of a process.
!!@see module_calcul::write_partial_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_partial_clusters_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: id
    INTEGER :: j
    INTEGER :: k
    INTEGER :: length
    INTEGER :: nb
    INTEGER, DIMENSION(:), POINTER :: matchings
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    DO i=0,params%nb_proc-1
       ! File name
       WRITE(num,*) i
       num=adjustl(num)
       length=len(trim(num))
       files='cluster.partiel.'//trim(num)
       OPEN(FILE=files,UNIT=10)
       ! Output
       files='cluster.partiel.'//num(1:length)//'.visu'
       OPEN(FILE=files,UNIT=1)
       WRITE(1,*) 'View "Clusters for part '//num(1:length)//'" {'
       READ(10,*) nb,k
       PRINT *, '> ', i, ' : ', nb, ' -> ', files
       ALLOCATE(coords(1,k))
       IF ((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) THEN
          ! Reading the matchings
          files='decoupe.'//trim(num)
          OPEN(FILE=files,UNIT=11)
          READ(11,*)
          ALLOCATE(matchings(nb))
          DO j=1,nb
             READ(11,*) matchings(j)
          ENDDO
          CLOSE(11)
       ENDIF
       DO j=1,nb
          IF (params%coords==1) THEN
             READ(10,*) coords(1,:),k
             CALL write_point_coord_format(params%dim, k, 1, 1, coords)
          ELSE
             READ(10,*) k,id
             CALL write_point_picture_format(params, id, matchings(k), 1)
          ENDIF
       ENDDO
       CLOSE(10) 
       WRITE(1,*) '};'
       CLOSE(1)
       DEALLOCATE(coords)
       IF ((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) THEN
          DEALLOCATE(matchings)
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE write_partial_clusters_gmsh



!>Writes the clusters after grouping for Gmsh visualization
!!@details This methods extracts details on computed clusters
!!from <em>cluster.final.x</em> files and writes them in 
!!<em>cluster.final.visu</em> file.
!!@note <em>x</em> is the identifier (number) of a cluster.
!!@see module_calcul::write_final_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_final_clusters_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    !=== IN/OUT ===
    !====  OUT ====

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: nb
    INTEGER :: nb0
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=params%input_file,UNIT=1)
    IF (params%coords==1) THEN
       ! Reading with coordinates format
       READ(1,*) j,k
       PRINT *, 'Reading mesh file...', j, k
       ALLOCATE(coords(j,k))
       nb0=0
       DO i=1,j
          READ(1,*,END=100) coords(i,:)
          nb0=nb0+1
       ENDDO
100    PRINT *, 'Number of points : ', nb0
       j=nb0
       CLOSE(1)
    ENDIF

    nb0=k ! Points dimension
    OPEN(FILE='cluster.final.visu',UNIT=1)
    WRITE(1,*) 'View "Clusters" {'
    ! Reading files
    DO i=1,params%nb_clusters
       ! File name
       WRITE(num,*) i
       files='cluster.final.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       READ(10,*) nb
       PRINT *,'  > ',i,' :',nb
       DO j=1,nb
          READ(10,*) k
          IF (params%coords==1) THEN
             ! Classic coordinates
             CALL write_point_coord_format(nb0, i, k, 1, coords)
          ELSE
             ! Pictures reassembly
             CALL write_point_picture_format(params, i, k, 1)
          ENDIF
       ENDDO
       CLOSE(10)
    ENDDO
    WRITE(1,*) '};'
    CLOSE(1)
    PRINT *,'-> cluster.final.visu'
    IF (params%coords==1) DEALLOCATE(coords)
    RETURN
  END SUBROUTINE write_final_clusters_gmsh



!>Writes point coordinates from coordinates format for Gmsh visualization
  SUBROUTINE write_point_coord_format(dim, id, k, unit, coords)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: dim
    INTEGER :: id
    INTEGER :: k
    INTEGER :: unit
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF (dim==2) THEN
       !2D
       WRITE(unit,*) 'SP(',coords(k,1),',',coords(k,2),',',0.,'){',id,'};' 
    ELSEIF (dim==3) THEN
       !3D
       WRITE(unit,*) 'SP(',coords(k,1),',',coords(k,2),',',coords(k,3),'){',id,'};'
    ENDIF
    RETURN
  END SUBROUTINE write_point_coord_format



!>Writes point coordinates from picture format for Gmsh visualization
!!@details The following is written in the file :</br>
!!<ol>
!!<li><b> For 2D case</b>: <em> SP(x,y){id} </em></li>
!!<li><b> For 3D case</b>: <em> SP(x,y,z){id} </em></li>
!!</ol>
!!@note x,y and z are the coordinates and id is the identifier of the process

!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] id the process identifier
!! @param[in] k 
!! @param[in] unit the written file identifier
  SUBROUTINE write_point_picture_format(params, id, k, unit)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    INTEGER :: id
    INTEGER :: k
    INTEGER :: unit

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: ix
    INTEGER :: iy
    DOUBLE PRECISION :: kx
    DOUBLE PRECISION :: ky
    DOUBLE PRECISION :: kz
    DOUBLE PRECISION, DIMENSION(:), POINTER :: data

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF (((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) &
         .AND.(params%image_dim==2)) THEN
       ALLOCATE(data(params%nb_points))
       OPEN(FILE=params%input_file,UNIT=50)
       READ(50,*)
       READ(50,*)
       DO i=1,params%nb_points
          READ(50,*) data(i)
       ENDDO
       CLOSE(50)
    ENDIF
    ! Coordinates
    ix=params%image_ref(k,1)
    iy=params%image_ref(k,2)
    IF (params%image_dim==2) THEN
       ! 2D points
       IF (params%is_geom==1) THEN
          kx=iy*params%step(2)
          ky=-ix*params%step(1)
       ELSE
          kx=float(iy)
          ky=-float(ix)
       ENDIF
       IF (((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) &
            .AND.(params%image_dim==2)) THEN
          kz=data(k)
       ELSE
          kz=0.
       ENDIF
    ELSEIF (params%image_dim==3) THEN
       ! 3D points
       IF (params%is_geom==1) THEN
          kx=iy*params%step(2)
          ky=-ix*params%step(1)
          kz=float(params%image_ref(k,3))*params%step(3)
       ELSE
          kx=float(iy)
          ky=-float(ix)
          kz=float(params%image_ref(k,3))
       ENDIF
    ENDIF
    ! Writing
    WRITE(unit,*) 'SP(',kx,',',ky,',',kz,'){',id,'};'
    RETURN
  END SUBROUTINE write_point_picture_format



!>Lists the commands related to Gmsh visualization
  SUBROUTINE list_commands_gmsh
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *,'gmsh decoupe.visu'
    PRINT *,'gmsh decoupe.geo'
    PRINT *,'gmsh decoupe.visu decoupe.geo'
    PRINT *,'gmsh cluster.partiel.*.visu'
    PRINT *,'gmsh cluster.final.visu'
    PRINT *,'gmsh decoupe.geo cluster.partiel.*.visu'
    PRINT *,'gmsh decoupe.geo cluster.final.visu'
    RETURN
  END SUBROUTINE list_commands_gmsh


END MODULE module_visuclusters_gmsh
