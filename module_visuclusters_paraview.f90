!>Contains methods enabling writing results in file specifically formatted for Paraview
MODULE module_visuclusters_paraview
  USE module_visuclusters_structure
CONTAINS



!>Writes the geometry partitioning for Paraview visualization
!!@details This methods extracts the domain definitions
!!from the <em>fort.2</em> file and writes to <em>visu/decoupe.*</em>
!!files.
!!@see module_calcul::write_partitioning()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_partitioning_paraview(params)
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
    DOUBLE PRECISION :: y0
    DOUBLE PRECISION :: y1
    DOUBLE PRECISION :: z0
    DOUBLE PRECISION :: z1
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x_min
    DOUBLE PRECISION, DIMENSION(:), POINTER :: y_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: y_min
    DOUBLE PRECISION, DIMENSION(:), POINTER :: z_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: z_min
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    ! Reading
    OPEN(FILE='fort.2',UNIT=2)
    ALLOCATE(x_min(params%nb_proc))
    ALLOCATE(x_max(params%nb_proc))
    ALLOCATE(y_min(params%nb_proc))
    ALLOCATE(y_max(params%nb_proc))
    ALLOCATE(z_min(params%nb_proc))
    ALLOCATE(z_max(params%nb_proc))
    DO i=1,params%nb_proc-params%is_interfacing
       IF (((params%coords==1).AND.(params%dim==2)) &
          .OR.((params%is_image==1).AND.(params%image_dim==2)) &
          .OR.((params%is_threshold==1).AND.(params%image_dim==2)) &
          .OR.((params%is_geom==1).AND.(params%image_dim+params%image_times==2))) THEN
          !2D
          IF (params%coords==1) THEN
             READ(2,*) x0,y0,num,x1,y1
             x_min(i)=x0
             y_min(i)=y0
             x_max(i)=x1
             y_max(i)=y1
          ELSEIF (params%is_threshold==1) THEN
             READ(2,*) x0,num,x1
             x_min(i)=1
             y_min(i)=-1
             x_max(i)=params%partitioning(2)
             y_max(i)=-params%partitioning(1)
             z_min(i)=x0
             z_max(i)=x1
          ELSEIF ((params%is_image==1).OR.(params%is_geom==1)) THEN
             READ(2,*) x0,y0,num,x1,y1
             x_min(i)=y0
             y_min(i)=-x0
             x_max(i)=y1
             y_max(i)=-x1
             z_min(i)=x0
          ENDIF
       ELSEIF (((params%coords==1).AND.(params%dim==3)) &
          .OR.((params%is_image==1).AND.(params%image_dim==3)) &
          .OR.((params%is_geom==1).AND.(params%image_dim+params%image_times==3))) THEN
          !3D
          IF (params%coords==1) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             x_min(i)=x0
             y_min(i)=y0
             z_min(i)=z0
             x_max(i)=x1
             y_max(i)=y1
             z_max(i)=z1
          ELSEIF ((params%is_image==1).OR.(params%is_geom==1)) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             x_min(i)=y0
             y_min(i)=-x0
             z_min(i)=z0
             x_max(i)=y1
             y_max(i)=-x1
             z_max(i)=z1
          ENDIF
       ENDIF
    ENDDO
    CLOSE(2)
    ! Writing
    PRINT *,'-> visu/decoupe.geo'
    PRINT *,'-> visu/decoupe.indices'
    OPEN(FILE='visu/decoupe.geo',UNIT=10)
    OPEN(FILE='visu/decoupe.indices',UNIT=11)
    WRITE(10,*) '** Ouput of  visuclusters **'
    WRITE(10,*) '** Partitioning of subclusters **'
    WRITE(10,'(a)') 'node id assign'
    WRITE(10,'(a)') 'element id assign'
    WRITE(11,*) '** Indexes of processes **'
    WRITE(10,'(a)') 'part'
    WRITE(10,*) 1
    WRITE(10,*) '** Partitionings **'
    WRITE(10,'(a)') 'Coordinates'
    IF (((params%coords==1).AND.(params%dim==2)) &
         .OR.((params%is_image==1).AND.(params%image_dim==2)) &
         .OR.((params%is_threshold==1).AND.(params%image_dim==2)) &
         .OR.((params%is_geom==1).AND.(params%image_dim+params%image_times==2))) THEN
       !2D
       IF ((params%coords==1).OR.(params%is_geom==1).OR.(params%is_image==1)) THEN
          WRITE(10,*) 4*(params%nb_proc-params%is_interfacing)
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) x_min(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_min(i)
          ENDDO
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_max(i)
             WRITE(10,*) y_max(i)
          ENDDO
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) 0.
             WRITE(10,*) 0.
             WRITE(10,*) 0.
             WRITE(10,*) 0.
          ENDDO
          WRITE(10,'(a)') 'quad4'
          WRITE(10,*) params%nb_proc-params%is_interfacing
          WRITE(11,'(a)') 'part'
          WRITE(11,*) 1
          WRITE(11,'(a)') 'quad4'
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) 4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*(i-1)+4
             WRITE(11,*) i
          ENDDO
       ELSEIF (params%is_threshold==1) THEN
          WRITE(10,*) 8*(params%nb_proc-params%is_interfacing)
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) x_min(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_min(i)
             WRITE(10,*) x_min(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_min(i)
          ENDDO
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_max(i)
             WRITE(10,*) y_max(i)
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_max(i)
             WRITE(10,*) y_max(i)
          ENDDO
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_max(i)
             WRITE(10,*) z_max(i)
             WRITE(10,*) z_max(i)
             WRITE(10,*) z_max(i)
          ENDDO
          WRITE(10,'(a)') 'hexa8'
          WRITE(10,*) params%nb_proc-params%is_interfacing
          WRITE(11,'(a)') 'part'
          WRITE(11,*) 1
          WRITE(11,'(a)') 'hexa8'
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) 8*(i-1)+1,8*(i-1)+2,8*(i-1)+3,8*(i-1)+4,&
                  8*(i-1)+5,8*(i-1)+6,8*(i-1)+7,8*(i-1)+8
             WRITE(11,*) i
          ENDDO
       ENDIF
    ELSEIF (((params%coords==1).AND.(params%dim==3)) &
         .OR.((params%is_image==1).AND.(params%image_dim==3)) &
         .OR.((params%is_geom==1).AND.(params%image_dim+params%image_times==3))) THEN
       !3D
       IF ((params%coords==1).OR.(params%is_geom==1).OR.(params%is_image==1)) THEN
          WRITE(10,*) 8*(params%nb_proc-params%is_interfacing)
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) x_min(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_min(i)
             WRITE(10,*) x_min(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_max(i)
             WRITE(10,*) x_min(i)
          ENDDO
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_max(i)
             WRITE(10,*) y_max(i)
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_min(i)
             WRITE(10,*) y_max(i)
             WRITE(10,*) y_max(i)
          ENDDO
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_min(i)
             WRITE(10,*) z_max(i)
             WRITE(10,*) z_max(i)
             WRITE(10,*) z_max(i)
             WRITE(10,*) z_max(i)
          ENDDO
          WRITE(10,'(a)') 'hexa8'
          WRITE(10,*) params%nb_proc-params%is_interfacing
          WRITE(11,'(a)') 'part'
          WRITE(11,*) 1
          WRITE(11,'(a)') 'hexa8'
          DO i=1,params%nb_proc-params%is_interfacing
             WRITE(10,*) 8*(i-1)+1,8*(i-1)+2,8*(i-1)+3,8*(i-1)+4,&
                  8*(i-1)+5,8*(i-1)+6,8*(i-1)+7,8*(i-1)+8
             WRITE(11,*) i
          ENDDO
       ENDIF
    ENDIF
    CLOSE(10)
    CLOSE(11)
    DEALLOCATE(x_min)
    DEALLOCATE(x_max)
    DEALLOCATE(y_min)
    DEALLOCATE(y_max)
    DEALLOCATE(z_min)
    DEALLOCATE(z_max)
    ! Master file
    PRINT *,'-> decoupe.CASE'
    OPEN(FILE='decoupe.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model:   visu/decoupe.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per element:   process   visu/decoupe.indices'
    CLOSE(12)
    RETURN
  END SUBROUTINE write_partitioning_paraview



!>Initializes the file of the partitioning for Paraview visualization
!!@details This method extracts details on partitioning from the
!!<em>decoupe.x</em> files and writes them in <em>visu/affectation.*</em>
!!files. Extra files are written in a case of interface partitioning.
!!<em>visu/affectation-interface.*</em> files are the assignment for
!!the interface.
!!@note <em>x</em> is the identifier (number) of a process.
!!@see module_calcul::write_partitioning(), partition_with_interface()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_assignment_paraview(params)
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
    INTEGER :: nb_points
    INTEGER :: nb_slaves
    INTEGER :: offset
    INTEGER, DIMENSION(:), POINTER :: ids
    INTEGER, DIMENSION(:), POINTER :: proc_ids
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    ! Reading of files
    IF (params%nb_proc==1) THEN
       offset=1
       nb_slaves=1
    ELSE
       offset=0
       nb_slaves=params%nb_proc-1
    ENDIF
    DO i=offset,nb_slaves
       IF ((i==0).AND.(params%is_interfacing==1).AND.(params%nb_proc>1)) THEN
          ! File aside for interfacing 
          PRINT *,'-> visu/affectation-interface.geo'
          PRINT *,'-> visu/affectation-interface.indices'
          OPEN(FILE='visu/affectation-interface.geo',UNIT=10)
          OPEN(FILE='visu/affectation-interface.indices',UNIT=11)
          WRITE(10,*) '** Output of visuclusters **'
          WRITE(10,*) '** Points on the interface **'
          WRITE(10,'(a)') 'node id assign'
          WRITE(10,'(a)') 'element id assign'
          WRITE(11,*) '** Indexes of processes **'
       ELSEIF (((i==1).AND.(params%is_interfacing==1)).OR. &
            ((i==0).AND.(params%is_interfacing==0))) THEN
          ! General file for others subdomains
          PRINT *,'-> visu/affectation.geo'
          PRINT *,'-> visu/affectation.indices'
          OPEN(FILE='visu/affectation.geo',UNIT=10)
          OPEN(FILE='visu/affectation.indices',UNIT=11)
          WRITE(10,*) '** Output of visuclusters **'
          WRITE(10,*) '** Partitioning of subclusters **'
          WRITE(10,'(a)') 'node id assign'
          WRITE(10,'(a)') 'element id assign'
          WRITE(11,*) '** Indexes of processes **'
       ENDIF
       ! File name
       WRITE(num,*) i
       files='decoupe.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=20)
       READ(20,*) nb_points
       PRINT *,'  > ',i,' :',nb_points
        ALLOCATE(coords(nb_points,params%dim))
       coords(:,:)=0.
       ALLOCATE(ids(nb_points))
       ids(:)=0
       ALLOCATE(proc_ids(nb_points))
       proc_ids(:)=0
       IF (nb_points>0) THEN
          DO j=1,nb_points
             IF (params%coords==1) THEN
                ! Partitioning by coordinates
                READ(20,*) coords(j,:)
                ids(j)=i
             ELSE
                ! Partitioning 1D picture
                READ (20,*) proc_ids(j)
                ids(j)=i
             ENDIF
          ENDDO
          ! Writing
          IF (params%coords==1) THEN
             ! Partitioning by coordinates
             CALL write_points_coord_format(params%dim, 1, nb_points, 10, 11, ids, coords)
          ELSE
             ! Partitioning 1D picture
             CALL write_points_picture_format(params, nb_points, 10, 11, ids, proc_ids)
          ENDIF
       ENDIF
       DEALLOCATE(coords)
       DEALLOCATE(ids)
       CLOSE(20)
       IF ((i==0).AND.(params%is_interfacing==1)) THEN
          ! Closes interfacing files
          CLOSE(10)
          CLOSE(11)
       ENDIF
    ENDDO
    CLOSE(10)
    CLOSE(11)
    ! Writes interfacing data
    IF ((params%is_interfacing==1).AND.(params%nb_proc>1)) THEN
       PRINT *,'-> affectation-interface.CASE'
       OPEN(FILE='affectation-interface.CASE',UNIT=12)
       WRITE(12,'(a)') 'FORMAT'
       WRITE(12,'(a)') 'type: ensight gold'
       WRITE(12,*) 
       WRITE(12,'(a)') 'GEOMETRY'
       WRITE(12,'(a)') 'model:   visu/affectation-interface.geo'
       WRITE(12,*) 
       WRITE(12,'(a)') 'VARIABLE'
       WRITE(12,'(a)') 'scalar per node:   process   visu/affectation-interface.indices'
       CLOSE(12)
    ENDIF
    ! writes others subdomains
    PRINT *,'-> affectation.CASE'
    OPEN(FILE='affectation.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model:   visu/affectation.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per node:   process   visu/affectation.indices'
    CLOSE(12)
    RETURN
  END SUBROUTINE write_assignment_paraview



!>Writes the clusters before grouping for Paraview visualization
!!@details This methods extracts details on computed clusters
!!on each domain from <em>cluster.partiel.x</em> files and writes
!!them in corresponding <em>visu/cluster.partiel.x.*</em> files.
!!@note <em>x</em> is the identifier (number) of a process.
!!@see module_calcul::write_partial_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_partial_clusters_paraview(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    CHARACTER (LEN=30) :: extension
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: length
    INTEGER :: nb_points
    INTEGER :: nb_zeros
    INTEGER, DIMENSION(:), POINTER :: ids
    INTEGER, DIMENSION(:), POINTER :: matchings
    INTEGER, DIMENSION(:), POINTER :: proc_ids
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    ! Number of generic characters to be used
    nb_zeros=floor(log(real(params%nb_proc-1))/log(real(10)))+1
    DO i=1,nb_zeros
       extension(i:i)='0'
    ENDDO
    DO i=0,params%nb_proc-1
       ! File name
       WRITE(num,*) i
       num=adjustl(num)
       length=len(trim(num))
       files='cluster.partiel.'//trim(num)
       OPEN(FILE=files,UNIT=20)
       ! Output
       ! General file for others subdomains
       files='cluster.partiel.'//extension(1:nb_zeros-length)//num(1:length)
       PRINT *,'-> visu/'//trim(files)//'.geo'
       PRINT *,'-> visu/'//trim(files)//'.indices'
       OPEN(FILE='visu/'//trim(files)//'.geo',UNIT=10)
       OPEN(FILE='visu/'//trim(files)//'.indices',UNIT=11)
       WRITE(10,*) '** Output of visuclusters **'
       WRITE(10,*) '** Partitioning subcluster '//trim(num)//' **'
       WRITE(10,'(a)') 'node id assign'
       WRITE(10,'(a)') 'element id assign'
       WRITE(11,*) '** Indexes of clusterized elements **'
       ! Reading
       READ(20,*) nb_points,k
       PRINT *,'  > ',i,' :',nb_points,' -> ',files
       ALLOCATE(coords(nb_points,k))
       ALLOCATE(ids(max(1,nb_points)))
       ids(:)=0
       ALLOCATE(proc_ids(max(1,nb_points)))
       proc_ids(:)=0
       IF ((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) THEN
          ! Reading the matchings
          files='decoupe.'//trim(num)
          OPEN(FILE=files,UNIT=21)
          READ(21,*)
          ALLOCATE(matchings(nb_points))
          DO j=1,nb_points
             READ(21,*) matchings(j)
          ENDDO
          CLOSE(21)
       ENDIF
       DO j=1,nb_points
          IF (params%coords==1) THEN
             READ(20,*) coords(j,:),ids(j)
          ELSE
             READ(20,*) proc_ids(j),ids(j)
             proc_ids(j)=matchings(proc_ids(j))
          ENDIF
       ENDDO
       CLOSE(20) 
       ! Writing
       IF (params%coords==1) THEN
          ! Partitioning by coordinates
          CALL write_points_coord_format(params%dim, 1, nb_points, 10, 11, ids, coords)
       ELSE
          ! Partitioning 1D picture
          CALL write_points_picture_format(params, nb_points, 10, 11, ids, proc_ids)
       ENDIF
       DEALLOCATE(coords)
       DEALLOCATE(ids)
       DEALLOCATE(proc_ids)
       CLOSE(10)
       CLOSE(11)
       IF ((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) THEN
          DEALLOCATE(matchings)
       ENDIF
    ENDDO
    PRINT *,'-> cluster.partiel.CASE'
    DO i=1,nb_zeros
       extension(i:i)='*'
    ENDDO
    OPEN(FILE='cluster.partiel.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model: 1 visu/cluster.partiel.'//extension(1:nb_zeros)//'.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per node: 1 cluster visu/cluster.partiel.'//&
         extension(1:nb_zeros)//'.indices'
    WRITE(12,*) 
    WRITE(12,'(a)') 'TIME'
    WRITE(12,'(a)') 'time set:         1'
    WRITE(num,*) params%nb_proc
    WRITE(12,'(a)') 'number of steps: '//trim(adjustl(num))
    WRITE(12,'(a)') 'filename start number: 0'
    WRITE(12,'(a)') 'filename increment: 1'
    WRITE(12,'(a)') 'time values:'
    DO i=0,params%nb_proc-1
       WRITE(12,*) i
    ENDDO
    CLOSE(12)
    RETURN
  END SUBROUTINE write_partial_clusters_paraview



!>Writes the clusters after grouping for Paraview visualization
!!@details This methods extracts details on computed clusters
!!from <em>cluster.final.x</em> files and writes them in 
!!<em>cluster.final.visu</em> file.
!!@note <em>x</em> is the identifier (number) of a cluster.
!!@see module_calcul::write_final_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE write_final_clusters_paraview(params)
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
    INTEGER :: nb_points
    INTEGER :: nb_points_temp
    INTEGER, DIMENSION(:), POINTER :: ids
    INTEGER, DIMENSION(:), POINTER :: proc_ids
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF (params%coords==1) THEN
       OPEN(FILE=params%input_file,UNIT=1)
       READ(1,*) j,k
       PRINT *, 'Reading mesh file...', j, k
       ALLOCATE(coords(j,k))
       nb_points_temp=0
       DO i=1,j
          READ(1,*,END=100) coords(i,:)
          nb_points_temp=nb_points_temp+1
       ENDDO
100    PRINT *, 'Number of points ', nb_points_temp
       j=nb_points_temp
       CLOSE(1)
    ENDIF

    nb_points_temp=k ! Points dimension
    ! Output
    PRINT *,'-> visu/cluster.final.geo'
    PRINT *,'-> visu/cluster.final.indices'
    OPEN(FILE='visu/cluster.final.geo',UNIT=10)
    OPEN(FILE='visu/cluster.final.indices',UNIT=11)
    WRITE(10,*) '** Output of visuclusters **'
    WRITE(10,*) '** Partitioning final cluster **'
    WRITE(10,'(a)') 'node id assign'
    WRITE(10,'(a)') 'element id assign'
    WRITE(11,*) '** clusters **'
    ALLOCATE(ids(max(1,params%nb_points)))
    ids(:)=0
    ALLOCATE(proc_ids(max(1,params%nb_points)))
    proc_ids(:)=0
    ! Reading files
    DO i=1,params%nb_clusters
       ! File name
       WRITE(num,*) i
       files='cluster.final.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=20)
       READ(20,*) nb_points
       PRINT *,'  > ',i,' :',nb_points
       DO j=1,nb_points
          READ(20,*) k
          ids(k)=i
          proc_ids(k)=k
       ENDDO
       CLOSE(20)
    ENDDO
    IF (params%coords==1) THEN
       ! Classical coordinates
       CALL write_points_coord_format(params%dim, 1, params%nb_points, 10, 11, ids, coords)
    ELSE
       ! Pictures reassembly
       CALL write_points_picture_format(params, params%nb_points, 10, 11, ids, proc_ids)
    ENDIF
    CLOSE(10)
    CLOSE(11)
    DEALLOCATE(ids)
    DEALLOCATE(proc_ids)
    IF (params%coords==1) DEALLOCATE(coords)
    ! Master file
    PRINT *,'-> cluster.final.CASE'
    OPEN(FILE='cluster.final.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model:   visu/cluster.final.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per element:   cluster   visu/cluster.final.indices'
    CLOSE(12)
    RETURN
  END SUBROUTINE write_final_clusters_paraview




!>Writes point coordinates from coordinates format for Paraview visualization
  SUBROUTINE write_points_coord_format(dim, k, nb_points, unit_geo, unit_ind, ids, coords)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: dim
    INTEGER :: k
    INTEGER :: nb_points
    INTEGER :: unit_geo
    INTEGER :: unit_ind
    INTEGER, DIMENSION(:), POINTER :: ids
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    
    !#### Variables  ####
    INTEGER :: i

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    WRITE(unit_geo,'(a)') 'part'
    WRITE(unit_geo,*) ids(k)
    WRITE(unit_geo,*) '** Partitionings **'
    WRITE(unit_geo,'(a)') 'coordinates'
    WRITE(unit_geo,*) nb_points
    WRITE(unit_ind,'(a)') 'part'
    WRITE(unit_ind,*) ids(k)
    WRITE(unit_ind,'(a)') 'point'
    DO i=1,nb_points
       WRITE(unit_geo,*) coords(i,1)
       WRITE(unit_ind,*) ids(i)
    ENDDO
    DO i=1,nb_points
       WRITE(unit_geo,*) coords(i,2)
    ENDDO
    DO i=1,nb_points
       IF (dim==2) THEN
          WRITE(unit_geo,*) 0.
       ELSE
          WRITE(unit_geo,*) coords(i,3)
       ENDIF
    ENDDO
    WRITE(unit_geo,'(a)') 'point'
    WRITE(unit_geo,*)nb_points
    DO i=1,nb_points
       WRITE(unit_geo,*) i
    ENDDO
    RETURN
  END SUBROUTINE write_points_coord_format



!>Writes point coordinates from picture format for Paraview visualization
!!@details The following is written in the geometric file for each point :
!!<ul>
!!<li><b> For 2D case</b>: <em> x y </em></li>
!!<li><b> For 3D case</b>: <em> x y z</em></li>
!!</ul>
!!And the following is written in the process partitioning file for each point :</br>
!!<em> ind </em>
!!@note x,y and z are the coordinates and ind is the identifier of the process.
!!@note This method writes headers at first in the files.

!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] nb_pixels the number of points
!! @param[in] unit_geo the written geometric file unit
!! @param[in] unit_ind the written processus identifiers file unit
!! @param[in] ids the process identifier for each point
!! @param[in] proc_ids the pixel identifiers
  SUBROUTINE write_points_picture_format(params, nb_pixels, unit_geo, unit_ind, ids, proc_ids)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    INTEGER :: nb_pixels
    INTEGER :: unit_geo
    INTEGER :: unit_ind
    INTEGER, DIMENSION(:), POINTER :: ids
    INTEGER, DIMENSION(:), POINTER :: proc_ids
    !=== IN/OUT ===
    !====  OUT ====
    
    !#### Variables  ####
    INTEGER :: i
    INTEGER :: ix
    INTEGER :: iy
    INTEGER :: k
    DOUBLE PRECISION, DIMENSION(:), POINTER :: data
    DOUBLE PRECISION, DIMENSION(:), POINTER :: kx
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ky
    DOUBLE PRECISION, DIMENSION(:), POINTER :: kz
    
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ALLOCATE(kx(nb_pixels))
    kx(:)=0
    ALLOCATE(ky(nb_pixels))
    ky(:)=0
    ALLOCATE(kz(nb_pixels))
    kz(:)=0
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
    ! Search the points
    DO i=1,nb_pixels
       k=proc_ids(i)
       ix=params%image_ref(k,1)
       iy=params%image_ref(k,2)
       IF (params%image_dim==2) THEN
          ! 2D points
          IF (params%is_geom==1) THEN
             kx(i)=iy*params%step(2)
             ky(i)=-ix*params%step(1)
          ELSE
             kx(i)=float(iy)
             ky(i)=-float(ix)
          ENDIF
          IF (((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) &
               .AND.(params%image_dim==2)) THEN
             kz(i)=data(k)
          ELSE
             kz(i)=0.
          ENDIF
       ELSEIF (params%image_dim==3) THEN
          ! 3D points
          IF (params%is_geom==1) THEN
             kx(i)=iy*params%step(2)
             ky(i)=-ix*params%step(1)
             kz(i)=float(params%image_ref(k,3))*params%step(3)
          ELSE

             kx(i)=float(iy)
             ky(i)=-float(ix)
             kz(i)=float(params%image_ref(k,3))
          ENDIF
       ENDIF
    ENDDO
    ! Writing
    WRITE(unit_geo,'(a)') 'part'
    WRITE(unit_geo,*) ids(1)
    WRITE(unit_geo,*) '** Partitionings **'
    WRITE(unit_geo,'(a)') 'coordinates'
    WRITE(unit_geo,*) nb_pixels
    WRITE(unit_ind,'(a)') 'part'
    WRITE(unit_ind,*) ids(1)
    WRITE(unit_ind,'(a)') 'point'
    DO i=1,nb_pixels
       WRITE(unit_geo,*) kx(i)
       WRITE(unit_ind,*) ids(i)
    ENDDO
    DO i=1,nb_pixels
       WRITE(unit_geo,*) ky(i)
    ENDDO
    DO i=1,nb_pixels
       WRITE(unit_geo,*) kz(i)
    ENDDO
    WRITE(unit_geo,'(a)') 'point'
    WRITE(unit_geo,*)nb_pixels
    DO i=1,nb_pixels
       WRITE(unit_geo,*) i
    ENDDO
    DEALLOCATE(kx)
    DEALLOCATE(ky)
    DEALLOCATE(kz)
    IF (((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) &
         .AND.(params%image_dim==2)) DEALLOCATE(data)
    RETURN
  END SUBROUTINE write_points_picture_format



!>Lists the commands related to Paraview visualization
  SUBROUTINE list_commands_paraview
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *,'paraview --data=decoupe.CASE'
    PRINT *,'paraview --data=affectation.CASE'
    PRINT *,'paraview --data=affectation-interface.CASE'
    PRINT *,'paraview --data=cluster.partiel.CASE'
    PRINT *,'paraview --data=cluster.final.CASE'
    RETURN
  END SUBROUTINE list_commands_paraview
  
END MODULE module_visuclusters_paraview
