PROGRAM clusters
  USE module_structure
  USE module_entree
  USE module_decoupe
  USE module_MPI
  USE module_calcul
  USE module_sortie

  IMPLICIT NONE
 
  ! MPI library
  INCLUDE 'mpif.h'
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  TYPE(type_clustering_param) :: clust_param
  TYPE(type_clusters), DIMENSION(:), POINTER :: array_clusters
  TYPE(type_data) :: data
  TYPE(type_data) :: partitioned_data
  CHARACTER (LEN=80) :: proc_name ! MPI variable
  CHARACTER (LEN=30) :: input
  CHARACTER (LEN=30) :: input_file
  INTEGER :: i
  INTEGER :: ierr ! MPI variable
  INTEGER :: j
  INTEGER :: length ! MPI variable
  INTEGER :: n_max
  INTEGER :: nb_clusters
  INTEGER :: nb_clusters_max
  INTEGER :: nb_clusters_opt
  INTEGER :: nb_proc ! MPI variable
  INTEGER :: proc_id ! MPI variable
  INTEGER :: status(MPI_STATUS_SIZE) ! MPI variable
  INTEGER :: tag ! MPI variable
  INTEGER :: nb_sd
  INTEGER,DIMENSION(:), POINTER :: list_nb_clusters
  INTEGER,DIMENSION(:), POINTER :: partitioning
  INTEGER,DIMENSION(:), POINTER :: points_by_cluster
  INTEGER,DIMENSION(:), POINTER :: points_by_domain
  INTEGER,DIMENSION(:,:), POINTER :: assignments
  INTEGER,DIMENSION(:,:), POINTER :: cluster_map
  DOUBLE PRECISION :: end_time
  DOUBLE PRECISION :: epsilon
  DOUBLE PRECISION :: sigma, delta
  DOUBLE PRECISION :: start_time
  DOUBLE PRECISION :: t1
  DOUBLE PRECISION :: t2
  DOUBLE PRECISION :: t_parall
  DOUBLE PRECISION :: t_parallg
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min
  DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
  LOGICAL :: exist_bool

  CHARACTER (LEN=20) :: filename, sproc_id

  INTEGER :: nb_slaves ! number of slaves by sub-domain
  INTEGER :: sd_com ! sub-domain communicator
  INTEGER :: master_com ! master communicator
  INTEGER :: nb_masters ! number of masters (should be equal to number of sub-domains)
  INTEGER :: lproc_id ! id of a MPI process in sub-domain communicator
  INTEGER :: mproc_id ! id of a MPI process in master communicator
  INTEGER, ALLOCATABLE :: masters(:) ! master array
  INTEGER :: world_group ! MPI_COMM_WORLD group
  INTEGER :: master_group ! com_master group
  INTEGER :: nb_points
  LOGICAL :: tab_nb_slaves(20)
 
  !###########################################
  ! INSTRUCTIONS
  !###########################################
  ! Timers init
  t1 = 0.0
  t2 = 0.0
  start_time = 0.0
  end_time = 0.0
  ! MPI init
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nb_proc, ierr)
  CALL MPI_GET_PROCESSOR_NAME(proc_name, length, ierr)
  PRINT *, 'Process ID : ', proc_id, '. Process name : ', proc_name, nb_proc

  ! available values for nb_slaves
  tab_nb_slaves(:) = .false.
  tab_nb_slaves(1) = .true.
  tab_nb_slaves(2) = .true.
  tab_nb_slaves(4) = .true.
  tab_nb_slaves(6) = .true.
  tab_nb_slaves(8) = .true.
  tab_nb_slaves(9) = .true.
  tab_nb_slaves(12) = .true.
  tab_nb_slaves(16) = .true.
  tab_nb_slaves(20) = .true.

  !=================== READ PARAM (0) =========================
  IF (proc_id==0) THEN

    t1 = MPI_WTIME()
    IF (iargc()>0) THEN
      ! Gets the input file name
      CALL getarg(1,input)
      INQUIRE(FILE=input,EXIST=exist_bool)
      IF (.NOT.exist_bool) CALL help
    ELSE
      ! If no input file, default=fort.1
      CALL help
    ENDIF
    ! File reading
    OPEN(FILE=input,UNIT=1)

    CALL read_params(nb_proc, data, clust_param, input_file, nb_sd, &
      nb_slaves, nb_clusters_max, list_nb_clusters, partitioning,  &
      epsilon, sigma, coord_max, coord_min)


    t2 = MPI_WTIME()
    PRINT *, 'Time for reading data : ', t2-t1
    CLOSE(1)
  ENDIF

  !=================== MPI COMMUNICATOR (All) =========================
  CALL MPI_BCAST(nb_slaves, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(nb_sd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  nb_masters = nb_proc/nb_slaves 
  IF(nb_sd /= nb_masters) THEN
    PRINT *, 'Problem : NB_SD /= NB_MASTERS'
    CALL MPI_FINALIZE(ierr)
    STOP
  END IF

  IF(nb_slaves > 20 .or. .not. tab_nb_slaves(nb_slaves)) THEN
    PRINT *, 'Problem : incompatible nb_slaves', nb_slaves
    CALL MPI_FINALIZE(ierr)
    STOP
  END IF

  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, proc_id/nb_slaves, mod(proc_id, nb_slaves), sd_com, ierr)
  CALL MPI_COMM_RANK(sd_com, lproc_id, ierr)

  ALLOCATE(masters(0:nb_masters-1))

  masters(0) = 0
  DO i = 1, nb_masters-1
    masters(i) = masters(i-1) + nb_slaves
  ENDDO

  IF(proc_id == 0) PRINT *, 'MASTERS', masters

  call MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierr)

  CALL MPI_GROUP_INCL(world_group, nb_masters, masters, master_group, ierr)

  CALL MPI_COMM_CREATE(MPI_COMM_WORLD, master_group, master_com, ierr)

  IF(proc_id==0) THEN
    start_time = MPI_WTIME()
  ENDIF

  !===================== MASTER's ZONE ======================================
  IF(mod(proc_id, nb_slaves) == 0) THEN

    call MPI_COMM_RANK(master_com, mproc_id, ierr)
#if aff
    PRINT *,'DEBUG : Launching process ',proc_id,' of ', nb_proc
#endif
    PRINT *,'----------------------------------------'


    ! Partitioning and sigma computing
    IF (proc_id==0) THEN
      ! Partitioning
#if aff
      PRINT *
      PRINT *,'DEBUG : Partitioning data...'
#endif
      PRINT *, 'NB_SD', nb_sd, nb_proc, nb_masters
      t1 = MPI_WTIME()
      CALL partition_data(data, nb_sd, partitioning, epsilon, coord_max, coord_min, points_by_domain, assignments, bounds)

      !    open(file="debug1.txt", unit = 22)
      !    write(22,*) 'coord'
      !    write(22,*) size(coord_min, 1), size(coord_max, 1)
      !    write(22,*) coord_min
      !    write(22,*) coord_max
      !
      !    write(22,*) 'partitioning'
      !    write(22,*) size(partitioning, 1)
      !    write(22,*) partitioning
      !
      !    write(22,*) 'points_by_domain'
      !    write(22,*) size(points_by_domain, 1)
      !    write(22,*) points_by_domain
      !
      !    write(22,*) 'assignments'
      !    write(22,*) size(assignments, 1)
      !    write(22,*) size(assignments, 2)
      !    do i = 1, size(assignments, 2)
      !      write(22,*) i, assignments(:, i)
      !    enddo
      !
      !    write(22,*) 'bounds'
      !    write(22,*) size(bounds, 1)
      !    write(22,*) size(bounds, 2)
      !    write(22,*) size(bounds, 3)
      !
      !    do i = 1, size(bounds, 1)
      !      do j = 1, size(bounds, 2)
      !        write(22,*) i,j, bounds(i,j, :)
      !      enddo
      !    enddo
      !
      !    close(22)


      t2 = MPI_WTIME()
      PRINT *,'Time for partitioning data : ', t2-t1
    ENDIF

    !  CALL MPI_FINALIZE(ierr)
    !  STOP

    IF(proc_id==0) THEN
      t1 = MPI_WTIME()
    ENDIF

    !  write(sproc_id, "(I2)") proc_id
    !  open(file="bug1_"//trim(adjustl(sproc_id)), unit = 22)

    ! Exchanges
    IF (nb_proc > 1) THEN
      ! Case of several proc 

      ! Nblimit sending
      CALL MPI_BCAST(nb_clusters_max, 1, MPI_INTEGER, 0, master_com, ierr)

      ! Nbideal sending
      IF (proc_id == 0) THEN
        DO i = 1, nb_masters-1
          tag = i
          CALL MPI_SEND(list_nb_clusters(i), 1, MPI_INTEGER, i, tag, master_com, ierr)
        ENDDO
        nb_clusters_opt = list_nb_clusters(0)
      ELSE
        tag = mproc_id
        CALL MPI_RECV(nb_clusters_opt, 1, MPI_INTEGER, 0, tag, master_com, status, ierr)
      ENDIF

      ! Data transferring
      IF (proc_id==0) THEN
        ! Data sending
#if aff
        PRINT *
        PRINT *,'DEBUG : Transferring partitioned data...'
#endif
        CALL send_partitioning(nb_masters, points_by_domain, assignments, master_com, data, partitioned_data)
#if aff
        PRINT *
        PRINT *,'DEBUG : Computing clusters...'
#endif
      ELSE
        ! Data receiving
        CALL receive_partitioning(mproc_id, master_com, partitioned_data)
      ENDIF

    ELSE
      ! Case of 1 proc alone
      partitioned_data%nb_points=data%nb_points
      partitioned_data%dim=data%dim
      partitioned_data%nb_clusters=0
      ALLOCATE(partitioned_data%points(data%nb_points))
      DO i=1,data%nb_points
      ALLOCATE(partitioned_data%points(i)%coords(data%dim))
      partitioned_data%points(i)%coords=data%points(i)%coords
      partitioned_data%points(i)%cluster=0
      ENDDO
      nb_clusters_opt=list_nb_clusters(1)
    ENDIF

    ! Sigma computing if auto individual mode
    IF(proc_id == 0) THEN
      IF (sigma < 0.0) THEN
        CALL get_delta(data, delta)
        CALL get_sigma(data, sigma)
        print *, '**** delta - sigma', delta, sigma
        sigma = delta / 2.D0 
      END IF
    ENDIF
    CALL MPI_BCAST(sigma, 1, MPI_DOUBLE_PRECISION, 0, master_com, ierr)
    clust_param%sigma = sigma

    IF(proc_id==0) THEN
      t2 = MPI_WTIME()
      PRINT *,'Time for sending data and computing sigma : ', t2-t1
    ENDIF

    !  write(22,*) 'sigma = ', sigma

    !  close(22)
    !  CALL MPI_FINALIZE(ierr)
    !  STOP

    ! Clusters computing
    ! Parallel part
    t1 = MPI_WTIME()
    IF (partitioned_data%nb_points>0) THEN
#if aff
      PRINT *,'DEBUG : ', proc_id, ' : computing clusters...'
#endif

      CALL send_clustering_param(clust_param, master_com)

      SELECT CASE (clust_param%clustering_method_id)
      CASE (1)
        PRINT *, 'Calling Spectral clustering...'
        CALL apply_spectral_clustering(clust_param, nb_clusters_max, nb_clusters_opt, proc_id, &
          sigma, partitioned_data, sd_com, nb_masters, nb_slaves)
      CASE (2)
        PRINT *, 'Calling Mean shift...'
        CALL mean_shift(proc_id,nb_clusters_max,nb_clusters_opt,partitioned_data,clust_param)
      CASE (3)
        PRINT *, 'Calling Kernel k means...'
        CALL apply_kernel_k_means(proc_id,nb_clusters_max,nb_clusters_opt,partitioned_data,clust_param)
      END SELECT

    ENDIF

    t2 = MPI_WTIME()
    t_parall = t2 - t1
    PRINT *, proc_id, ' : computing parallel cluster : ', t_parall

    CALL MPI_REDUCE(t_parall, t_parallg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, master_com, ierr)
    IF(proc_id==0) PRINT *, 'Time for computing global parallel cluster :', t_parallg

    IF(proc_id==0) THEN
      t1 = MPI_WTIME()
    ENDIF

    ! Saves the partial clusters
    CALL write_partial_clusters(proc_id, partitioned_data)

    ! Exchanges part
    IF (nb_proc>1) THEN
      ! Clusters grouping
      IF (proc_id == 0) THEN
#if aff
        PRINT *
        PRINT *,'DEBUG : Grouping clusters...'
#endif
        ! Receiving of the number of clusters with duplications
        CALL receive_number_clusters(partitioned_data, nb_masters, points_by_domain, master_com, &
          nb_clusters, array_clusters)
#if aff
        PRINT *,'DEBUG : number of clusters with duplications found : ', nb_clusters
#endif
        ! Receiving of clusters info
        PRINT *, nb_clusters, data%nb_points
        ALLOCATE(cluster_map(nb_clusters,data%nb_points))
        print *, 'after allocate'
        CALL receive_clusters(array_clusters, partitioned_data, nb_clusters, nb_sd, points_by_domain,  &
          assignments, master_com, points_by_cluster, cluster_map)
      ELSE
        ! Sends the number of clusters
        CALL send_number_clusters(partitioned_data, mproc_id, master_com)
        ! Sends the clusters
        CALL send_clusters(partitioned_data, mproc_id, master_com)
      ENDIF

      ! End of post-process
      IF (proc_id==0) THEN
        ! Groups the clusters and removes duplicates from the set of found clusters
        CALL group_clusters(nb_clusters, points_by_cluster, cluster_map, data)
      ENDIF

    ELSE
      ! Case of 1 proc alone
      nb_clusters=partitioned_data%nb_clusters
      ALLOCATE(points_by_cluster(nb_clusters))
      points_by_cluster(:)=0
      n_max=0
      DO i=1,partitioned_data%nb_points
        j=partitioned_data%points(i)%cluster
        points_by_cluster(j)=points_by_cluster(j)+1
        n_max=max(n_max,points_by_cluster(j))
      ENDDO
      ALLOCATE(cluster_map(nb_clusters,n_max))
      cluster_map(:,:)=0
      points_by_cluster(:)=0
      DO i=1,partitioned_data%nb_points
        j=partitioned_data%points(i)%cluster
        points_by_cluster(j)=points_by_cluster(j)+1
        cluster_map(j,points_by_cluster(j))=i
      ENDDO
    ENDIF

    IF(proc_id==0) THEN
      t2 = MPI_WTIME()
      PRINT *,'Time for grouping clusters : ', t2-t1
    ENDIF

    IF(proc_id==0) THEN
      t1 = MPI_WTIME()
    ENDIF
    ! Outputs
    IF (proc_id==0) THEN
      ! Writing of the cluster.final files
      PRINT*, 'Calling write_final_clusters'
      CALL write_final_clusters(points_by_cluster, cluster_map, nb_clusters)

      ! Information writing
      CALL write_metadata(data, input_file, nb_clusters, nb_sd)
    ENDIF
    IF(proc_id==0) THEN
      t2 = MPI_WTIME()
      PRINT *,'Time for writing clusters : ', t2-t1
    ENDIF

    ! Ending of MPI
    IF (proc_id==0) THEN
      end_time = MPI_WTIME()
      PRINT *,'End of computing'
      PRINT *,'Overall time : ', end_time-start_time
    ENDIF

  !========================= END of MASTER CODE ===========================
  ELSE
  !========================= LES AUTRES ===========================
    !print *,'*/*/*/*/*/*/*/*/ je suis un autre'
    !CALL slave_solve_pdsyev(sd_com, nb_sd, nb_slaves)

  END IF
  
  CALL MPI_FINALIZE(ierr)
  STOP

END PROGRAM clusters
