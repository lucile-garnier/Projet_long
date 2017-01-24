!>Contains methods enabling data and parameter file reading
MODULE module_entree
 USE module_structure
CONTAINS

!>Displays help on used formats and keywords for <em>param.in</em> file
  SUBROUTINE help
    IMPLICIT NONE
    include 'mpif.h'

    INTEGER :: ierr
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *
    PRINT *,'Calling syntax : clusters input_file'
    PRINT *
    PRINT *,'Input file keywords : '
    PRINT *
    PRINT *,'DATA'
    PRINT *,'COORD (if data with coordinates)'
    PRINT *,'IMAGE (if mesh file = image + partitioning by pixel)'
    PRINT *,'GEOM  (if mesh file = image + geom partitioning)'
    PRINT *,'SEUIL (if mesh file = image + threshold partitioning)'
    PRINT *,'  mesh_file'
    PRINT *
    PRINT *,'EPAISSEUR'
    PRINT *,'  thickness_of_the_slice'
    PRINT *
    PRINT *,'NBLIMIT'
    PRINT *,'  max_nb_of_clusters'
    PRINT *
    PRINT *,'NBCLUST'
    PRINT *,'  (facultatif)'
    PRINT *,'  nb_of_clusters_by_subdomain'
    PRINT *
    PRINT *, 'CLUSTMETHID'
    PRINT *,'  (optional)'
    PRINT *,'  clustering_method_id : 1 Spectral, 2 MeanShift, 3 K-KMeans'
    PRINT *
    PRINT *, 'NBFINALCLUST'
    PRINT *,'  (optional)'
    PRINT *,'  nb_final_clust desired cluster number for KKMeans'
    PRINT *
    PRINT *, 'KERNEL'
    PRINT *,'  (optional)'
    PRINT *,'  ker_table 0 polynomial, 1 gaussian'
    PRINT *
    PRINT *,'SIGMA'
    PRINT *,'  (optional)'
    PRINT *,'  imposed_value_of_sigma'
    PRINT *
    PRINT *, 'GAMMA'
    PRINT *,'  (optional)'
    PRINT *,'  ker_table'
    PRINT *
    PRINT *, 'DELTA'
    PRINT *,'  (optional)'
    PRINT *,'  ker_table'
    PRINT *
    PRINT *, 'BANDWIDTH'
    PRINT *,'  (optional)'
    PRINT *,'  bandwidth_for_mean_shift'
    PRINT *
    PRINT *,'DECOUPAGE'
    PRINT *,'INTERFACE (partitioning with interface) '
    PRINT *,'RECOUVREMENT (partitioning with overlapping)'
    PRINT *,'  nb_of_subdomains_by_DIMENSION'
    PRINT *
    PRINT *,'END'
    PRINT *,'  (end of input file)'
    CALL MPI_ABORT(ierr)
    STOP
    RETURN
  END SUBROUTINE help


!>Reads a file in which there is the whole information on the data
!! @param[in] nb_proc the number of processors used
!! @param[in,out] data the entire data for computing
!! @param[out] input_file the name of the text file where input data is written
!! @param[out] nb_sd number of sub-domains
!! @param[out] nb_slaves number of slaves on each sub-domains
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
!! @param[out] epsilon the slice thickness
!! @param[out] sigma the affinity parameter
!! @param[out] nb_clusters_max the maximum number of clusters
!! @param[out] list_nb_clusters the imposed numbers of clusters list for testing purpose (useless)
!! @param[out] partitioning the partitionning (number of processors along each dimension)
  SUBROUTINE read_params(nb_proc, data, clust_param, input_file, nb_sd, nb_slaves, &
               nb_clusters_max, list_nb_clusters, partitioning, epsilon, &
               sigma, coord_max, coord_min)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nb_proc

    !=== IN/OUT ===
    TYPE(type_data) :: data
    TYPE(type_clustering_param) :: clust_param

    !====  OUT ====
    CHARACTER (LEN=30) :: input_file
    INTEGER :: nb_sd
    INTEGER :: nb_slaves
    INTEGER :: nb_clusters_max
    INTEGER, DIMENSION(:), POINTER :: list_nb_clusters
    INTEGER, DIMENSION(:), POINTER :: partitioning
    DOUBLE PRECISION :: epsilon
    DOUBLE PRECISION :: sigma
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min
    

    !#### Variables  ####
    CHARACTER (LEN=30) :: word
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: tot
    LOGICAL :: ok
    LOGICAL :: partitioning_bool

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    epsilon=0.0
    sigma=-1.0
    data%coords=0
    data%is_image=0
    data%is_geom=0
    data%is_threshold=0
    data%is_interfacing=0
    data%is_overlapping=0


    nb_sd = -1
    nb_slaves = 1
    nb_clusters_max=4
    partitioning_bool=.FALSE.

    ! Reading
    ok=.FALSE.
    DO WHILE (.NOT. ok)
       ok=.TRUE.
       READ(1,*) word
       PRINT *, word
       SELECT CASE(word)
       CASE('DATA')
          ok=.FALSE.
          READ(1,*) input_file
          IF (input_file=='IMAGE') THEN
             data%is_image=1
             READ (1,*) input_file
             PRINT *, '> Input image format + partitioning by pixel'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_picture_data(input_file, data, coord_max, coord_min)
          ELSEIF (input_file=='GEOM') THEN
             data%is_geom=1
             READ (1,*) input_file
             PRINT *, '> Input image format + geometric partitioning'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_geometric_data(input_file, data, coord_max, coord_min)
          ELSEIF (input_file=='SEUIL') THEN
             data%is_threshold=1
             READ (1,*) input_file
             PRINT *, '> Input image format + partitioning by threshold'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_threshold_data(input_file, data, coord_max, coord_min)
          ELSEIF (input_file=='COORD') THEN
             data%coords=1
             READ (1,*) input_file
             PRINT *, '> Input image format + partitioning by threshold'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_coordinates_data(input_file, data, coord_max, coord_min)
          ELSE
             PRINT *
             PRINT *, 'Non-recognized data format !!!'
             CALL help
          ENDIF
          IF ((data%is_image==1).OR.(data%is_geom==1).OR.(data%is_threshold==1)) THEN
             ! Creation of array pixels/coordinates
             PRINT *, '> Decoding image format...'
             CALL assign_picture_array(data)
          ENDIF
       CASE('EPAISSEUR')
          ok=.FALSE.
          READ(1,*) epsilon
          PRINT *, '> Thickness of the slice :', epsilon
       CASE('NBLIMIT')
          ok=.FALSE.
          READ(1,*) nb_clusters_max
          PRINT *, '> Maximal number of searched clusters : ', nb_clusters_max
          clust_param%nbLimitClust = nb_clusters_max
       CASE('NBCLUST')
          ok=.FALSE.
          IF(nb_sd < 0) THEN
            PRINT *, 'Problem NBCLUST'
            CALL help
          ELSE
            READ(1,*) list_nb_clusters(:)
            PRINT *, '> Test for number of clusters=', list_nb_clusters
          END IF
       CASE('KERNEL')
          ok=.FALSE.
          READ(1,*) clust_param%kernelfunindex
          PRINT *, '> Kernel fun index =', clust_param%kernelfunindex
       CASE('CLUSTMETHID')
          ok=.FALSE.
          READ(1,*) clust_param%clustering_method_id
          PRINT *, '> clustering method id =', clust_param%clustering_method_id
       CASE('GAMMA')
          ok=.FALSE.
          READ(1,*) clust_param%gamma
          PRINT *, '> 	gamma =', clust_param%gamma
       CASE('DELTA')
          ok=.FALSE.
          READ(1,*) clust_param%delta
          PRINT *, '> 	delta =', clust_param%delta
       CASE('SIGMA')
          ok=.FALSE.
          READ(1,*) sigma
          PRINT *, '> Imposed value of sigma : ', sigma
          clust_param%sigma=sigma
          IF (data%is_image==1) THEN
             IF (sigma<1.0) THEN
                PRINT *, 'Too small thickness for image mode !!!!'
                STOP
             ENDIF
          ENDIF
       CASE('BANDWIDTH')
          ok=.FALSE.
          READ(1,*) clust_param%bandwidth
          PRINT *, '> 	BANDWIDTH =', clust_param%bandwidth
       CASE('DECOUPAGE')
          partitioning_bool=.TRUE.
          ok=.FALSE.
          READ (1,*) word
          SELECT CASE(word)
          CASE('INTERFACE')
             data%is_interfacing=1
             PRINT *, '> Partitioning by interface activated.'
          CASE('RECOUVREMENT')
             data%is_overlapping=1
             PRINT *, '> Partitioning by overlapping activated.'
          CASE DEFAULT
             PRINT *
             PRINT *, 'Bad partitioning format !!!'
             PRINT *
             CALL help
          END SELECT

          PRINT *, 'dim : ', data%dim

          IF ((data%coords==1).OR.(data%is_geom==1).OR.(data%is_threshold==1)) THEN
             ALLOCATE(partitioning(data%dim))
          ELSEIF (data%is_image==1) THEN
             ! Partitioning per pixel
             ALLOCATE(partitioning(data%image_dim))
          ENDIF
          READ(1,*) partitioning(:)
          PRINT *, 'partitioning', partitioning
          IF (nb_proc>1) THEN
             nb_sd = 1
             IF ((data%coords==1).OR.(data%is_geom==1).OR.(data%is_threshold==1)) THEN
                DO i=1,data%dim
                    nb_sd =nb_sd*partitioning(i)
                ENDDO
             ELSEIF (data%is_image==1) THEN
                ! Partitioning per pixel
                DO i=1,data%image_dim
                   nb_sd = nb_sd*partitioning(i)
                ENDDO
             ENDIF
             !IF (tot/=nb_proc-data%is_interfacing) THEN
             !   PRINT *, 'Invalidated partitioning !'
             !   PRINT *, 'Number of process must be equal to ', tot+data%is_interfacing
             !   CALL MPI_ABORT(ierr)
             !   STOP
             !ENDIF
          ELSE
             ! 1 proc
             IF ((data%coords==1).OR.(data%is_geom==1).OR.(data%is_threshold==1)) THEN
                DO i=1,data%dim
                   partitioning(i)=1
                ENDDO
                nb_sd = 1
             ELSEIF (data%is_image==1) THEN
                ! Partitioning per pixel
                DO i=1,data%image_dim
                   partitioning(i)=1
                ENDDO
                nb_sd = 1
             ENDIF
          ENDIF
          IF (nb_sd>1) THEN
            ALLOCATE(list_nb_clusters(0:nb_sd-1))
          ELSE
            ALLOCATE(list_nb_clusters(1))
          ENDIF
          list_nb_clusters(:)=0
          IF(mod(nb_proc, nb_sd) /= 0) THEN
            PRINT *, 'Problem mod(nb_proc, nb_sd) /= 0'
            STOP
          END IF
          nb_slaves = nb_proc / nb_sd
          PRINT *, '> partitioning :', partitioning
          PRINT *, '> nb_slaves :', nb_slaves
       CASE('END')
          ok=.TRUE.
       CASE DEFAULT
          ok=.FALSE.
          PRINT *, 'Unknown keyword : ', word
       END SELECT
    ENDDO

    ! Partitioning parameter
    IF ((nb_sd > 1).AND.(.NOT.partitioning_bool)) THEN
       PRINT *
       PRINT *, 'The keyword <<DECOUPAGE>> has not been found !'
       CALL help 
    ENDIF

    ! 1 proc
    IF (nb_proc==1) THEN
       ! Initialization to 1 by default of all the partitioning parameters
       IF (partitioning_bool) DEALLOCATE(partitioning)
       IF ((data%coords==1).OR.(data%is_geom==1).OR.(data%is_threshold==1)) THEN
          ALLOCATE(partitioning(data%dim) )
       ELSEIF (data%is_image==1) THEN
          ALLOCATE(partitioning(data%image_dim))
       ENDIF
       partitioning(:)=1
       epsilon=1.0
    ENDIF   
    ! Validation of the combinations of input parameters
    tot = data%is_geom+data%is_threshold+data%coords+data%is_image
    IF (tot /= 1) THEN
       PRINT *
       PRINT *, 'Problem with data format !'
       CALL help
    ENDIF
    tot = data%is_interfacing+data%is_overlapping
    IF (tot /= 1) THEN
       PRINT *
       PRINT *, 'Problem with data format !'
       CALL help
    ENDIF
    RETURN
  END SUBROUTINE read_params


!>Reads data written in coordinates format
  SUBROUTINE read_coordinates_data(input_file, data, coord_max, coord_min)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Reading classic data
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%nb_points,data%dim
    data%nb_clusters=0
    PRINT *, '> Number of points : ', data%nb_points
    PRINT *, '> Dimension : ', data%dim
    ALLOCATE(data%points(data%nb_points))
    ALLOCATE(coord_max(data%dim))
    ALLOCATE(coord_min(data%dim))
    nb=0
    DO i=1,data%nb_points
       ALLOCATE(data%points(i)%coords(data%dim))
       READ(2,*,END=100) data%points(i)%coords(:)
       nb=nb+1
       data%points(i)%cluster=-1
       IF (i==1) THEN
          coord_max(:)=data%points(1)%coords(:)
          coord_min(:)=data%points(1)%coords(:)
       ELSE
          DO j=1,data%dim
             coord_min(j)=min(coord_min(j),data%points(i)%coords(j))
             coord_max(j)=max(coord_max(j),data%points(i)%coords(j))
          ENDDO
       ENDIF
    ENDDO
100 PRINT *, 'Number of points : ',nb
    data%nb_points=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates : '
    DO j=1,data%dim
       PRINT *, '> ', j, ' : ', coord_min(j), coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_coordinates_data


!>Reads data written in picture format
!!@details This type of file starts with two <em>integers</em>
!!separated by a blank space: the first one corresponds to the 
!!dimension of the image and the second one the number of attributes
!!(typically, the attributes could be the color channel intensities).
!!The following line is composed of two <em>integers</em> separated by 
!!a blank space: the first one corresponds to the number of pixels in
!!a column and the second one to the number of pixels in a row. 
!!Then, next lines are composed of <em>double</em> numbers separated 
!!by blank spaces. Each of these lines corresponds to the value of
!!the attributes for each pixel. The coordinates of the pixels are 
!!implicit as they are ordered by rows.
!!@note The attributes are stored in the field type_data::coord and so ONLY the color will be taken into account.
!!@note The outputs <em>coordmax</em> and <em>coordmin</em> are computed according to the position of the pixels.
!!@see assign_picture_array(), read_coordinates_data(), read_threshold_data(), read_geometric_data
!! @param[in] input_file the name of the text file where input data is written
!! @param[in,out] data the entire data for computing
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
  SUBROUTINE read_picture_data(input_file, data, coord_max, coord_min)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%image_dim,data%image_times
    PRINT *, '> Image dimension : ', data%image_dim
    PRINT *, '> Number of time : ', data%image_times
    ALLOCATE(data%partitioning(data%image_dim))
    READ(2,*) data%partitioning(:)
    PRINT *, '> Spatial partitioning : ', data%partitioning
    data%nb_points=1
    DO i=1,data%image_dim
       data%nb_points=data%nb_points*data%partitioning(i)
    ENDDO
    data%dim=data%image_times
    data%nb_clusters=0
    PRINT *, '> Number of points to read : ', data%nb_points
    ALLOCATE(data%points(data%nb_points))
    ALLOCATE(coord_max(data%image_dim))
    ALLOCATE(coord_min(data%image_dim))
    coord_min(:)=0.9
    DO i=1,data%image_dim
       coord_max(i)=data%partitioning(i)+0.1
    ENDDO
    nb=0
    DO i=1,data%nb_points
       ALLOCATE(data%points(i)%coords(data%dim))
       READ(2,*,END=200) data%points(i)%coords(:)
       nb=nb+1
       data%points(i)%cluster=-1
    ENDDO
200 PRINT *, '> Number of points read : ', nb       
    data%nb_points=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates :'
    DO j=1,data%dim
       PRINT *, '> ', j, ' : ', coord_min(j), coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_picture_data
  

!>Reads data written in geometric format
!!@details This type of file starts with two <em>integers</em>
!!separated by a blank space: the first one corresponds to the 
!!dimension of the image and the second one the number of attributes
!!(typically, the attributes could be the color channel intensities).
!!The following line is composed of two <em>integers</em> separated by 
!!a blank space: the first one corresponds to the number of pixels in
!!a column and the second one to the number of pixels in a row. 
!!Then, next lines are composed of <em>double</em> numbers separated 
!!by blank spaces. Each of these lines corresponds to the value of
!!the attributes for each pixel. The coordinates of the pixels are 
!!implicit as they are ordered by rows.
!!@note The attributes and the coordinates are stored in the field type_data::coord and so the position AND the color will be taken into account.
!!@note The partitioning is made according to the coordinates
!!@see assign_picture_array(), read_coordinates_data(), read_threshold_data(), read_picture_data
!! @param[in] input_file the name of the text file where input data is written
!! @param[in,out] data the entire data for computing
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
  SUBROUTINE read_geometric_data(input_file, data, coord_max, coord_min)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb
    DOUBLE PRECISION :: max_step

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%image_dim,data%image_times
    PRINT *, '> Image dimension : ', data%image_dim
    PRINT *, '> Number of time : ', data%image_times
    ALLOCATE(data%step(data%image_dim))
    data%step(:)=0.0
    ALLOCATE(data%partitioning(data%image_dim))
    READ(2,*) data%partitioning(:)
    PRINT *, '> Spatial partitioning : ', data%partitioning
    data%nb_points=1
    DO i=1,data%image_dim
       data%nb_points=data%nb_points*data%partitioning(i)
    ENDDO
    data%dim=data%image_dim+data%image_times
    data%nb_clusters=0
    PRINT *,'> Number of points to read : ', data%nb_points
    ALLOCATE(data%points(data%nb_points))
    ALLOCATE(coord_max(data%dim))
    ALLOCATE(coord_min(data%dim))
    nb=0
    DO i=1,data%nb_points
       ALLOCATE(data%points(i)%coords(data%dim))
       data%points(i)%coords(:)=0.0
       READ(2,*,END=300) data%points(i)%coords(data%image_dim+1:data%image_dim+data%image_times)
       nb=nb+1
       data%points(i)%cluster=-1
       IF (i==1) THEN
          coord_max(:)=data%points(1)%coords(:)
          coord_min(:)=data%points(1)%coords(:)
       ELSE
          DO j=1,data%dim
             coord_min(j)=min(coord_min(j),data%points(i)%coords(j))
             coord_max(j)=max(coord_max(j),data%points(i)%coords(j))
          ENDDO
       ENDIF
    ENDDO
300 PRINT *, '> Number of points read : ', nb
    data%nb_points=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates : '
    max_step=1.e-13
    DO j=data%image_dim+1,data%image_dim+data%image_times
       max_step=max(max_step,coord_max(j)-coord_min(j))
       PRINT *, '> ', j, ' : ', coord_min(j), coord_max(j)
    ENDDO
    PRINT *,'> Maximal step : ', max_step
    ! Searching steps by picture dimension
    DO j=1,data%image_dim
       data%step(j)=max_step/data%partitioning(j)
       PRINT *, '> Step : ', j, data%step(j)
       coord_min(j)=0.9*data%step(j)
       coord_max(j)=(data%partitioning(j)+1)*data%step(j)
    ENDDO
    RETURN
  END SUBROUTINE read_geometric_data



!>Reads data written in threshold format
!!@details This type of file starts with two <em>integers</em>
!!separated by a blank space: the first one corresponds to the 
!!dimension of the image and the second one the number of attributes
!!(typically, the attributes could be the color channel intensities).
!!The following line is composed of two <em>integers</em> separated by 
!!a blank space: the first one corresponds to the number of pixels in
!!a column and the second one to the number of pixels in a row. 
!!Then, next lines are composed of <em>double</em> numbers separated 
!!by blank spaces. Each of these lines corresponds to the value of
!!the attributes for each pixel. The coordinates of the pixels are 
!!implicit as they are ordered by rows.
!!@note The attributes are stored in the field type_data::coord and so ONLY the color will be taken into account
!!@note The outputs <em>coordmax</em> and <em>coordmin</em> are computed according to the color of the pixels.
!!@see assign_picture_array(), read_picture_data(), read_coordinates_data(), read_geometric_data
!! @param[in] input_file the name of the text file where input data is written
!! @param[in,out] data the entire data for computing
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
  SUBROUTINE read_threshold_data(input_file, data, coord_max, coord_min)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Reading classic data
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%image_dim,data%image_times
    PRINT *, '> Image dimension : ', data%image_dim
    PRINT *, '> Number of time : ', data%image_times
    ALLOCATE(data%partitioning(data%image_dim))
    READ(2,*) data%partitioning(:)
    PRINT *, '> Spatial dimension : ', data%partitioning
    data%nb_points=1
    DO i=1,data%image_dim
       data%nb_points=data%nb_points*data%partitioning(i)
    ENDDO
    data%dim=data%image_times
    data%nb_clusters=0
    PRINT *, '> Number of points to read : ', data%nb_points
    ALLOCATE(data%points(data%nb_points))
    ALLOCATE(coord_max(data%dim))
    ALLOCATE(coord_min(data%dim))
    nb=0
    DO i=1,data%nb_points
       ALLOCATE(data%points(i)%coords(data%dim))
       READ(2,*,END=400) data%points(i)%coords(:)
       nb=nb+1
       data%points(i)%cluster=-1
       IF (i==1) THEN
          coord_max(:)=data%points(1)%coords(:)
          coord_min(:)=data%points(1)%coords(:)
       ELSE
          DO j=1,data%dim
             coord_min(j)=min(coord_min(j),data%points(i)%coords(j))
             coord_max(j)=max(coord_max(j),data%points(i)%coords(j))
          ENDDO
       ENDIF
    ENDDO
400 PRINT *, 'Number of points : ', nb
    data%nb_points=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates : '
    DO j=1,data%dim
       PRINT *, '> ', j, ' : ' , coord_min(j), coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_threshold_data



!>Puts the index number of the pixels into an array
!!@details Each point of the data set is mapped
!!to an array corresponding to an image. Thus, instead
!!of one index for accessing one point, two will be used
!!(row and column).
!! @param[in,out] data the entire data for computing
  SUBROUTINE assign_picture_array(data)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !=== IN/OUT ===
    TYPE(type_data) :: data

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER, DIMENSION(:), POINTER :: plane
    LOGICAL :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Creation of array points/image_coordinates
    ALLOCATE(data%image_ref(data%nb_points,data%image_dim))
    ALLOCATE(plane(data%image_dim))
    plane(:)=1
    DO i=1,data%nb_points
       DO j=1,data%image_dim
          ! Index in the array points/pixel
          data%image_ref(i,j)=plane(j)
          IF (data%is_geom==1) THEN
             ! Input of coordinates 1:imgdim for the geometric cluster
             data%points(i)%coords(j)=plane(j)*data%step(j)
          ENDIF
       ENDDO
       ok=.FALSE.
       k=data%image_dim
       DO WHILE(.NOT. ok)
          IF (plane(k)<data%partitioning(k)) THEN
             plane(k)=plane(k)+1
             ok=.TRUE.
          ELSE
             plane(k)=1
             k=k-1
          ENDIF
          IF (k==0) ok=.TRUE.
       ENDDO
    ENDDO
    DEALLOCATE(plane)
    RETURN
  END SUBROUTINE assign_picture_array


END MODULE module_entree
