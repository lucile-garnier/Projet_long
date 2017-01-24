!>Contains methods enabling writing results in a selected data file format (for now : Paraview or GMSH)
MODULE module_visuclusters
  USE module_visuclusters_structure
  USE module_visuclusters_gmsh
  USE module_visuclusters_paraview
CONTAINS


!>Reads a file containing metadata on the data and the computed clusters
!!@details This function extracts the following information from the input 
!!<em>fort.3</em> file :
!!<ol>
!!<li> The data file name </li>
!!<li> The number of the points in the entire data set </li>
!!<li> The number of processes used </li>
!!<li> The partitioning mode (by interface or overlapping) </li>
!!<li> The number of clusters found </li>
!!<li> The data file format </li>
!!</ol>
!!In the case of a picture format: this extra information is
!!written :
!!<ol>
!!<li> The image dimension </li>
!!<li> The image partitioning </li>
!!<li> The number of attributes </li>
!!<li> The number of steps (only in geometric format) </li>
!!</ol> 
!!@see write_metadata()
!! @param[in,out] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE read_metadata(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !=== IN/OUT ===
    TYPE(type_params) :: params
    
    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: n
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################     
    params%is_geom=0
    params%is_threshold=0
    PRINT *
    PRINT *, 'Reading info...'
    READ (3,*)
    READ(3,*) params%input_file
    PRINT *, '> Mesh file name : ', params%input_file
    READ (3,*)
    READ(3,*) params%nb_points
    PRINT *, '> Number of points : ', params%nb_points
    READ (3,*)
    READ (3,*) params%dim
    PRINT *, '> Dimension : ', params%dim
    READ (3,*)
    READ(3,*) params%nb_proc
    PRINT *, '> Number of process : ', params%nb_proc
    READ (3,*)
    READ(3,*) params%is_interfacing
    PRINT *, '> Partitioning by interfacing ? : ', params%is_interfacing
    READ (3,*)
    READ(3,*) params%is_overlapping
    PRINT *, '> Partitioning by overlapping ? : ', params%is_overlapping
    READ (3,*)
    READ(3,*) params%nb_clusters  
    PRINT *, '> Number of clusters got : ', params%nb_clusters
    READ (3,*)
    READ(3,*) params%coords
    PRINT *, '> Coordinates format ? : ', params%coords
    READ (3,*)
    READ(3,*) params%is_image
    PRINT *, '> Image format ? : ', params%is_image
    READ (3,*)
    READ(3,*) params%is_geom
    PRINT *, '> Geometric format ? : ', params%is_geom
    READ (3,*)
    READ(3,*) params%is_threshold
    PRINT *, '> Threshold format ? : ', params%is_threshold
    IF ((params%is_image==1).OR.(params%is_geom==1).OR.(params%is_threshold==1)) THEN
       READ (3,*)
       READ(3,*) params%image_dim
       PRINT *, '> Image dimension : ', params%image_dim
       ALLOCATE(params%partitioning(params%image_dim))
       READ (3,*)
       READ(3,*) params%partitioning(:)
       PRINT *, '> Image partitioning : ', params%partitioning
       READ (3,*)
       READ(3,*) params%image_times
       PRINT *, '> Number of time : ', params%image_times
       ! Points referencing
       ALLOCATE(params%image_ref(params%nb_points,params%image_dim))
       params%image_ref(:,:)=0
       n=0
       DO i=1,params%partitioning(1)
          DO j=1,params%partitioning(2)
             IF (params%image_dim==2) THEN
                n=n+1
                params%image_ref(n,1)=i
                params%image_ref(n,2)=j
             ELSE
                DO k=1,params%partitioning(3)
                   n=n+1
                   params%image_ref(n,1)=i
                   params%image_ref(n,2)=j
                   params%image_ref(n,3)=k
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       IF (params%is_geom==1) THEN
          ALLOCATE(params%step(params%image_dim))
          READ (3,*)
          READ (3,*) params%step(:)
          PRINT *, '> Mesh step : ', params%step(:)
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE read_metadata


!>Writes the geometry of the partitioning (Gmsh or Paraview) and calls the eponym function
!!@details This methods extracts the domain definitions
!!from the <em>fort.2</em> file.
!!@see module_calcul::write_partitioning()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
  SUBROUTINE write_partitioning(params, format_output)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    PRINT *
    PRINT *, 'Writing partitioning geometry...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_partitioning_gmsh(params)
    CASE('paraview')
       CALL write_partitioning_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_partitioning




!>Initializes the file of the partitionning
!!@details This method extracts details on partitioning from the
!!<em>decoupe.x</em> files.
!!@see module_calcul::write_partial_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
  SUBROUTINE write_assignment(params, format_output)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    PRINT *
    PRINT *, 'Writing partitioning allocations...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_assignment_gmsh(params)
    CASE('paraview')
       CALL write_assignment_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_assignment




!>Writes the clusters before grouping by calling the corresponding method (Gmsh or Paraview)
!!@details This methods extracts details on computed clusters
!!on each domain from <em>cluster.partiel.x</em> files.
!!@see module_calcul::write_partial_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
  SUBROUTINE write_partial_clusters(params, format_output)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output

    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    PRINT *
    PRINT *, 'Reading clusters before grouping...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_partial_clusters_gmsh(params)
    CASE('paraview')
       CALL write_partial_clusters_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_partial_clusters




!>Writes the clusters after grouping by calling the corresponding method (Gmsh or Paraview)
!!@details This methods extracts details on computed clusters
!!from <em>cluster.final.x</em> files.
!!@see module_calcul::write_final_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
  SUBROUTINE write_final_clusters(params, format_output)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    PRINT *
    PRINT *, 'Reading clusters after grouping...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_final_clusters_gmsh(params)
    CASE('paraview')
       CALL write_final_clusters_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_final_clusters




!>Lists the commands related to Gmsh or Paraview depending on the input parameter
!! @param[in] format_output the file format for visualization
  SUBROUTINE list_commands(format_output)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: format_output
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    PRINT *
    PRINT *, '-------------------------------------'
    PRINT *, 'Command list of visualisation : '
    PRINT *
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL list_commands_gmsh
    CASE('paraview')
       CALL list_commands_paraview
    END SELECT
    RETURN
  END SUBROUTINE list_commands


END MODULE module_visuclusters
