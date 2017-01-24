!>Contains the spectral clustering method and methods that computes affinity parameters for kernels and overlapping
MODULE module_calcul
  USE module_structure
  USE module_solve
  USE module_embed
CONTAINS

!>Computes the distance delta
!! @param[in] partitioned_data the partitioned data for computing
!! @param[out] delta, the distance
  SUBROUTINE get_delta(partitioned_data, delta)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data

    !====  OUT ====
    DOUBLE PRECISION :: delta

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    DOUBLE PRECISION :: produit
    DOUBLE PRECISION :: minx, maxx

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    produit = 1.D0
    DO i = 1, partitioned_data%dim

       minx = partitioned_data%points(i)%coords(1)
       maxx = partitioned_data%points(i)%coords(1)
       DO j = 2, partitioned_data%nb_points
         minx = min(minx, partitioned_data%points(j)%coords(i))
         maxx = max(maxx, partitioned_data%points(j)%coords(i))
       ENDDO
       produit = produit * (maxx-minx)
    ENDDO
    delta = (produit / partitioned_data%nb_points)**(1.D0 / partitioned_data%dim)
    ! Safety
    IF (delta == 0.0) delta = 1.0
    RETURN
  END SUBROUTINE get_delta


!>Computes the affinity parameter @f$sigma@f$
!!@details The Gaussian affinity matrix is widely used and depends on a free parameter
!!@f$\sigma@f$. It is known that this parameter affects the results in spectral clustering
!!and spectral embedding. With an assumption that the @f$p@f$ dimensionnal data set S composed
!!of @f$n@f$ points is isotropic enough, this data set is included in a @f$p@f$ dimensional 
!!box bounded by @f$D_{max}@f$ the largest distance between pairs of points in @f$S@f$ : 
!!@f{equation}{D_{max} = \max_{1\leq i,j\leq n} \Vert x_i - x_j \Vert @f}
!!So a reference distance noted @f$\sigma@f$ could be defined. This distance represents the
!!case of an uniform distribution in the sense that all pair of points are seprated by the
!!same distance @f$\sigma@f$ in the box of edge size @f$D_{max}@f$ :
!!@f{equation}{\sigma = \frac{D_{max}}{n^{1/p}} @f}
!!<br>
!!This function is used to compute this parameter on one domain. The algorithm is simple :
!!<ol> 
!!<li> <b> Find the maximum distance using two nested loops over the points in the domain </b> </li>
!!<li> <b> Divide it by the @f$p^{th}@f$ root of @f$n@f$ </b> </li>
!!</ol>
!! @param[in] partitioned_data the partitioned data for computing
!! @param[out] sigma the affinity parameter
  SUBROUTINE get_sigma(partitioned_data, sigma)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data

    !====  OUT ====
    DOUBLE PRECISION :: sigma

    !#### Variables  ####
    INTEGER :: i1
    INTEGER :: j1
    INTEGER :: k1
    DOUBLE PRECISION :: norm

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    sigma=0.0
    DO i1=1,partitioned_data%nb_points
       DO j1=i1+1,partitioned_data%nb_points
          norm=0.0
          DO k1=1,partitioned_data%dim
             norm=norm+&
                  (partitioned_data%points(i1)%coords(k1)-partitioned_data%points(j1)%coords(k1))**2
          ENDDO
          sigma=max(sigma,sqrt(norm))
       ENDDO
    ENDDO
    sigma=sigma/(2*exp(log(float(partitioned_data%nb_points))*(1.0/float(partitioned_data%dim))))
    ! Safety
    IF (sigma==0.0) sigma=1.0
    RETURN
  END SUBROUTINE get_sigma


!>Computes the affinity parameter @f$sigma@f$ for the interface
!!@details This method is useful when the partitioning is
!!made by interface. Because, the domain defining the interface
!!has a volume whose topology changes drastically, a specific
!!computation of @f$\sigma@f$ has to be made.
!!@deprecated Use get_sigma() instead
!! @param[in] partitioned_data the partitioned data for computing
!! @param[in] bounds the intervals along each dimension representing the bounds of each partition
!! @param[in] epsilon the slice thickness
!! @param[in] proc_id the processus identifier
!! @param[in] partitioning the partitionning (number of processors along each dimension)
!! @param[out] sigma the affinity parameter
  SUBROUTINE get_sigma_interface(partitioned_data, proc_id, partitioning, epsilon, bounds, sigma)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data
    INTEGER :: proc_id
    INTEGER, DIMENSION(:), POINTER :: partitioning
    DOUBLE PRECISION :: epsilon
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds

    !====  OUT ====
    DOUBLE PRECISION :: sigma

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: nb
    INTEGER, DIMENSION(:,:), POINTER :: array
    INTEGER, DIMENSION(:), POINTER :: partitioning_tmp
    DOUBLE PRECISION :: ext_volume
    DOUBLE PRECISION :: int_volume
    DOUBLE PRECISION :: long
    DOUBLE PRECISION :: sigma0

    !###########################################
    ! INSTRUCTIONS
    !###########################################
#if aff
    ! Number of partitionings
    nb=1
    DO i=1,partitioned_data%dim
       nb=nb*partitioning(i)
    ENDDO
    ! Creation of partitioning
    ALLOCATE(array(nb,0:partitioned_data%dim))
    ALLOCATE(partitioning_tmp(partitioned_data%dim))
    partitioning_tmp(:)=1
    DO i=1,nb
       DO j=1,partitioned_data%dim
          array(i,j)=partitioning_tmp(j)
       ENDDO
       partitioning_tmp(1)=partitioning_tmp(1)+1
       k=1
       DO WHILE(partitioning_tmp(k)>partitioning(k))
          partitioning_tmp(k)=1
          IF (k<partitioned_data%dim) partitioning_tmp(k+1)=partitioning_tmp(k+1)+1
       ENDDO
    ENDDO
    DEALLOCATE(partitioning_tmp)
    ! Value of sigma
    sigma0=0.0
    DO i=1,nb
       ext_volume=1.0
       int_volume=1.0
       DO j=1,partitioned_data%dim
          k=array(i,j)
          long=bounds(j,k,2)-bounds(j,k,1)
          ext_volume=ext_volume*long
          int_volume=int_volume*max(0.0D1,long-2.0*epsilon)
       ENDDO
       sigma0=sigma0+ext_volume-int_volume
    ENDDO
    DEALLOCATE(array)
    ! Computing of scale length
    sigma0=exp(1.0/float(partitioned_data%dim)*log(sigma0))
    ! Sigma computing
    sigma0=sigma0/(2.0*exp(log(float(partitioned_data%nb_points))*(1.0/float(partitioned_data%dim))))
    PRINT *, 'DEBUG : ', proc_id, ' : value of computed sigma for interfacing : ', sigma0
#endif
    ! Sigma computing, global formula
    CALL get_sigma(partitioned_data, sigma)
#if aff
    PRINT *, 'DEBUG : ', proc_id,' : value of sigma for interfacing', sigma
#endif
    RETURN
  END SUBROUTINE get_sigma_interface


  FUNCTION poly_kernel( partitioned_data, gam, delta )
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data
    DOUBLE PRECISION :: delta
    DOUBLE PRECISION :: gam 

    !====  OUT  ====
    DOUBLE PRECISION, DIMENSION(partitioned_data%nb_points,partitioned_data%nb_points) :: poly_kernel

    !#### Variables  ####
    INTEGER :: k
    INTEGER :: i
    INTEGER :: j
    INTEGER :: n
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Ker
    DOUBLE PRECISION :: aux

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    n=partitioned_data%nb_points
    ALLOCATE(Ker(n,n))
    Ker(:,:) = 0.D0

    DO i = 1, n-1
      DO j = 1, n-1
        aux = 0.D0
        DO k = 1, partitioned_data%dim
          aux = aux + partitioned_data%points(i)%coords(k) * partitioned_data%points(j)%coords(k)
        ENDDO 
        Ker(i,j) = (aux + gam)**delta
      ENDDO
    ENDDO
    poly_kernel = Ker
    RETURN
  END FUNCTION


  FUNCTION gaussian_kernel( partitioned_data, sigma )
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data
    DOUBLE PRECISION sigma

    !====  OUT  ====
    DOUBLE PRECISION, DIMENSION(partitioned_data%nb_points,partitioned_data%nb_points) :: gaussian_kernel

    !#### Variables  ####
    INTEGER k
    INTEGER i
    INTEGER j
    INTEGER n
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Ker
    DOUBLE PRECISION :: aux

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    n = partitioned_data%nb_points
    ALLOCATE(Ker(n,n))
    !ALLOCATE(gaussian_kernel)
    Ker(:,:) = 0.D0

    DO i = 1, n-1
      DO j = i+1, n
        aux = 0.D0
        DO k = 1, partitioned_data%dim
          aux = aux + (partitioned_data%points(i)%coords(k) - partitioned_data%points(j)%coords(k))**2
        ENDDO
        !Ker(i,j) = exp(- aux / (2*sigma**2))
        Ker(i,j) = exp(- aux / ((sigma/2.D0)**2))
        ! Symetry
        Ker(j,i) = Ker(i,j)
      ENDDO
    ENDDO
    gaussian_kernel = Ker
    RETURN
  END FUNCTION
!Stop when converged compute E = sum_N(sum_M( Indicatrice (xi E Ck)*||phi(xi)-mk||ÃÂÃÂ²))



SUBROUTINE apply_kernel_k_means(proc_id,nb_clusters_max,nb_clusters_opt, &
                  partitioned_data,clust_param)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_clustering_param) :: clust_param
    INTEGER :: nb_clusters_max
    INTEGER :: nb_clusters_opt
    INTEGER :: proc_id

    !=== IN/OUT ===
    TYPE(type_data) :: partitioned_data

    !#### Variables  ####
    INTEGER :: cluster(partitioned_data%nb_points) ! indicates which cluster each point belongs to
    INTEGER :: cluster_id (partitioned_data%nb_clusters)
    INTEGER :: cluster_population (partitioned_data%nb_clusters) ! the number of points in each cluster
    INTEGER :: i
    INTEGER :: it_max ! the maximum number of iterations
    INTEGER :: it_num ! the number of iterations taken
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l
    INTEGER :: p
    INTEGER :: stockpopulation (clust_param%nbLimitClust)
    INTEGER :: swap
    DOUBLE PRECISION :: cluster_center (partitioned_data%dim, clust_param%nbLimitClust) ! the cluster centers
    DOUBLE PRECISION :: den1
    DOUBLE PRECISION :: den2
    DOUBLE PRECISION :: listnorm (partitioned_data%nb_points, clust_param%nbLimitClust)
    DOUBLE PRECISION :: norm
    DOUBLE PRECISION :: num1
    DOUBLE PRECISION :: num2
    DOUBLE PRECISION :: seuil
    DOUBLE PRECISION :: stockcenter (partitioned_data%dim, clust_param%nbLimitClust)
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: valmax
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Ker
    LOGICAL :: ok 
    LOGICAL :: ok2

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  

    ALLOCATE(Ker(partitioned_data%nb_points,partitioned_data%nb_points))


    Ker(:,:)=0.0


    PRINT*,'clust_param%nbLimitClust : ',clust_param%nbLimitClust
    it_num = 0
 
    !###########################################      
    ! Kernel Calculus
    !########################################### 
    IF (clust_param%kernelfunindex==0) THEN
        Ker=poly_kernel( partitioned_data, clust_param%gamma, clust_param%delta)
    ELSEIF (clust_param%kernelfunindex==1) THEN
        Ker=gaussian_kernel(partitioned_data, clust_param%sigma)
    ENDIF

    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.

    !  Assign one point to each cluster center.
    cluster_center(:,1) = partitioned_data%points(1)%coords(:) !point(:,1) %%
    cluster_id(:)=0
    cluster_id(1)=1
    p=2
    seuil=0.4
    !###########################################      
    PRINT *, 'Center research'
    !###########################################      
    DO i = 2, clust_param%nbLimitClust
       ok=.FALSE.
       DO WHILE(.NOT.ok)
          valmax=2.0*seuil
          !Search if the point is already used as a center
          ok2=.FALSE.
          DO j=1,i-1
             IF (partitioned_data%points(j)%cluster==p) ok2=.TRUE.
          ENDDO
          !If not a center check if it is out of a threshold
          IF (.NOT.ok2) THEN
             DO j=1,i-1
                val=0.0
                norm=0.0
                DO k=1,partitioned_data%dim
                   val=max(val,abs(cluster_center(k,j)-partitioned_data%points(p)%coords(k))) 
                  !VOIR SI CELA DOIT ETRE MODIFIE EN FONCTION DES KERNEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ENDDO
                valmax=min(val,valmax)
             ENDDO
             IF (valmax>=seuil) THEN
             ok=.TRUE.
             ENDIF
          ENDIF
         p=p+1

         !Lower the threshold if not enough center are found after going through all the dataset
         IF ((p>partitioned_data%nb_points).AND.(.NOT.ok)) THEN 
            ! Reinitialise the process with a lower threshold
            seuil=0.9*seuil
            p=1
            PRINT *,'Lower threshold :',seuil
         ENDIF
       ENDDO
       ! We are out of the loop => the previous point was out previous center threshold=> it s a new cluster center
       p=p-1
       cluster_center(:,i)= partitioned_data%points(p)%coords(:) !point(:,p) 
       cluster_id(i)=p
    ENDDO
    !###########################################      
    PRINT *,'Points checked to find the nb cluster Max : ',p
    !###########################################      

!!! boucle            
    it_num = 0
    swap=1
    !Initialise the cluster points are part of
    partitioned_data%points(:)%cluster=1 

    DO WHILE ((it_num<partitioned_data%nb_points**2).AND.(swap/=0))
   
       it_num = it_num + 1
       swap=0


       !! Computing the distances
       cluster_population(1:clust_param%nbLimitClust) = 1
       listnorm(:,:)=0.0
       num1=0.0
       den1=0.0
       num2=0.0
       den2=0.0
       DO k=1,clust_param%nbLimitClust
           DO i=1,partitioned_data%nb_points
               DO j=1,partitioned_data%nb_points
                   IF ( partitioned_data%points(j)%cluster.EQ.k) THEN
                   num1=num1 + 2*(Ker(i,j))
                   den1=den1+1 
                   ENDIF
                   DO l=1,partitioned_data%nb_points
                       IF ( (partitioned_data%points(j)%cluster.EQ.k).AND.(partitioned_data%points(l)%cluster.EQ.k)) THEN
                       num2=num2 + Ker(j,l)
                       den2=den2+1
                       ENDIF
                   ENDDO
               ENDDO
           listnorm(i,k)= Ker(i,i)
               IF (den1.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k) - (num1/den1)
               ENDIF
               IF (den2.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k)  + (num2/den2)
               ENDIF
           ENDDO
       ENDDO

       !!Min distance assignation 
       cluster_population(:)=0
       DO i=1,partitioned_data%nb_points

          DO j=1,clust_param%nbLimitClust
             IF (listnorm(i,j)<listnorm(i,partitioned_data%points(i)%cluster).AND.( partitioned_data%points(i)%cluster.NE.j)) THEN
                partitioned_data%points(i)%cluster=j
                swap=swap+1
                ! A point change its cluster
             ENDIF
          ENDDO

          cluster_population(partitioned_data%points(i)%cluster)=cluster_population(partitioned_data%points(i)%cluster)+1
       ENDDO
       
    ENDDO
    partitioned_data%nb_clusters=clust_param%nbLimitClust
    DEALLOCATE(Ker)
    PRINT*,'OUT apply_kernel_k_means   clust_param%nbLimitClust : ',clust_param%nbLimitClust
    !RETURN

  END SUBROUTINE apply_kernel_k_means


  SUBROUTINE apply_spectral_clustering(clust_param, nb_clusters_max, nb_clusters_opt, proc_id, sigma, &
          partitioned_data, sd_com, nb_masters, nb_slaves)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_clustering_param) :: clust_param
    INTEGER :: nb_clusters_max
    INTEGER :: nb_clusters_opt
    INTEGER :: proc_id
    DOUBLE PRECISION :: sigma
    INTEGER :: sd_com
    INTEGER :: nb_masters
    INTEGER :: nb_slaves

    !=== IN/OUT ===
    TYPE(type_data) :: partitioned_data

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: n
    INTEGER :: nb
    INTEGER :: nb_clusters
    INTEGER :: nbvp
    INTEGER :: solver ! solveur au valeur propre => parametre de controle
    INTEGER, DIMENSION(:), POINTER :: cluster
    INTEGER, DIMENSION(:), POINTER :: nb_info
    INTEGER, DIMENSION(:), POINTER :: points_by_clusters
    DOUBLE PRECISION :: norm
    DOUBLE PRECISION :: ratio
    DOUBLE PRECISION :: ratio1
    DOUBLE PRECISION :: ratio2
    DOUBLE PRECISION :: seuilrij
    DOUBLE PRECISION :: t1
    DOUBLE PRECISION :: t2
    DOUBLE PRECISION :: t_cons_vp
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: value
    DOUBLE PRECISION, DIMENSION(:), POINTER :: clusters_energies
    DOUBLE PRECISION, DIMENSION(:), POINTER :: D
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomin
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomoy
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiorii
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiorij
    DOUBLE PRECISION, DIMENSION(:), POINTER :: W
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A2
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: clusters_centers
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z

    CHARACTER (LEN=20) :: sproc_id
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Matrix creation
    PRINT *, proc_id, ' : value of sigma : ', sigma
    n=partitioned_data%nb_points
    ! Forall i, A(i,i) = 0
    ALLOCATE(A(n,n))
    A(:,:)=0.0

    IF (clust_param%kernelfunindex==0) THEN
        A=poly_kernel( partitioned_data, clust_param%gamma, clust_param%delta)
    ELSEIF (clust_param%kernelfunindex==1) THEN
        A=gaussian_kernel(partitioned_data, clust_param%sigma)
    ENDIF

    ! Normalizing of affinity matrix
    ALLOCATE(D(n))
    D(:)=0.0
    DO i=1,n
       DO j=1,i-1
          D(i)=D(i)+A(i,j)
       ENDDO
       DO j=i+1, n
          D(i)=D(i)+A(j,i)
       ENDDO
    ENDDO

    ! D-1 * A
    DO i=1,n
       DO j=1,n
          ! Matrix A is not symmetric anymore
          !A(i,j)=A(i,j)/D(i)
       ENDDO
    ENDDO

    ! D -1/2 * A * D -1/2 
    DO i=1,n
       DO j=1,n
          A(i,j)=A(i,j)/sqrt(D(i))
       ENDDO
    ENDDO
    DO j=1,n
       DO i=1,n
          A(i,j)=A(i,j)/sqrt(D(j))
       ENDDO
    ENDDO


    DEALLOCATE(D)

    solver = 0


    !write(sproc_id, "(I2)") proc_id
    !open(file="A_"//trim(adjustl(sproc_id)), unit = 22)

    !DO i=1,n
    !   DO j=1,n
    !      write(22,*) i, j, A(i,j)
    !   ENDDO
    !ENDDO

    !close(22)

    IF(solver == 0) THEN

      nbvp = n

      ALLOCATE(A2(n,n))
      A2(:,:)=0.0
      DO i=1,n
        DO j=1,n
          A2(i,j)=A(i,j)
        ENDDO
      ENDDO

      t1 = MPI_WTIME()

      if(nb_slaves == 1) then
        print *, proc_id, 'Lapack solver'
        CALL solve_dgeev(n, A2, W, Z)
      else
        print *, proc_id, 'Scalapack solver'
        !CALL master_solve_pdsyev(n, A2, sd_com, nb_masters, nb_slaves, W, Z)
      endif

!      write(sproc_id, "(I2)") proc_id
!      open(file="Z_"//trim(adjustl(sproc_id)), unit = 22)
!
!      DO i=1,n
!        DO j=1,n
!          write(22,*) i, j, Z(i,j)
!        ENDDO
!      ENDDO

    ELSE

      PRINT *, proc_id, ' : Arpack solver'

      nb = 2*nb_clusters_max
      nbvp = nb
      CALL solve_arpack_full(A, n, nb, W, Z)

    ENDIF

    t2 = MPI_WTIME()

    t_cons_vp = t2 - t1
    PRINT *, proc_id, ' : Time for eigen values construction : ', t_cons_vp

!    maintenant solve_dgeev et master_solve renvoient les valeurs propres rangées
!    ATTENTION à ARPACK
!    DO i=1,nbvp-1
!       DO j=i+1,nbvp
!          IF (W(i)<W(j)) THEN
!             val=W(i)
!             W(i)=W(j)
!             W(j)=val
!             DO k=1,n
!                val=Z(k,i)
!                Z(k,i)=Z(k,j)
!                Z(k,j)=val
!             ENDDO
!          ENDIF
!       ENDDO
!    ENDDO

    ! Spectral embedding
    IF ((nb_clusters_opt==0).AND.(n>2)) THEN
       ! Search of the best partitioning
       ALLOCATE(ratiomax(nb_clusters_max))
       ratiomax(:)=0
       ALLOCATE(ratiomin(nb_clusters_max))
       ratiomin(:)=0
       ALLOCATE(ratiomoy(nb_clusters_max))
       ratiomoy(:)=0
       ALLOCATE(ratiorii(nb_clusters_max))
       ratiorii(:)=0
       ALLOCATE(ratiorij(nb_clusters_max))
       ratiorij(:)=0

       ALLOCATE(nb_info(nb_clusters_max))
       nb_info(:)=0
       DO nb_clusters=2,min(n,nb_clusters_max)

          ALLOCATE(cluster(n))
          cluster(:)=0
          ALLOCATE(clusters_centers(nb_clusters,nb_clusters))
          clusters_centers(:,:)=0.0
          ALLOCATE(points_by_clusters(nb_clusters))
          points_by_clusters(:)=0
          ALLOCATE(clusters_energies(nb_clusters))
          clusters_energies(:)=0.0

          CALL apply_spectral_embedding(n, nb_clusters, proc_id, A, Z,  &
                  nb_info(nb_clusters), cluster, points_by_clusters,  &
                  ratiomax(nb_clusters), ratiomoy(nb_clusters), ratiorii(nb_clusters), &
                   ratiorij(nb_clusters), clusters_centers, clusters_energies)


          DEALLOCATE(cluster)
          DEALLOCATE(clusters_centers)
          DEALLOCATE(clusters_energies)
          DEALLOCATE(points_by_clusters)
       ENDDO

#if aff
PRINT *, 'DEBUG : Frobenius ratio'
#endif
       ! Norm of frobenius ratio
       ratio=ratiomax(nb_clusters_max)
       partitioned_data%nb_clusters=nb_clusters_max
       ratio1=0.0
       ratio2=1e+10
       DO i=2,nb_clusters_max
          seuilrij=1e-4
          IF ((ratiorii(i)>=0.95*ratio1).AND.(ratiorij(i)-ratio2<=seuilrij)) THEN  
             partitioned_data%nb_clusters=i
             ratio1=ratiorii(i)
             ratio2=ratiorij(i)
          ENDIF
       ENDDO

    ELSEIF ((nb_clusters_opt==1).AND.(n>nb_clusters_opt)) THEN
       ! Test with an imposed cluster
       ALLOCATE(nb_info(nb_clusters_opt))
       nb_info(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       partitioned_data%nb_clusters=nb_clusters_opt
    ELSE
       ! Case of a domain with less points than nb_clusters_opt or only one point
       ALLOCATE(nb_info(n))
       nb_info(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       partitioned_data%nb_clusters=n
       ALLOCATE(ratiomax(n))
       ratiomax(:)=0
       ALLOCATE(ratiomoy(n))
       ratiomoy(:)=0
       ALLOCATE(ratiomin(n))
       ratiomin(:)=0
       ALLOCATE(ratiorii(n))
       ratiorii(:)=0
       ALLOCATE(ratiorij(n))
       ratiorij(:)=0
    ENDIF
    ! Case of nb_clusters==1
    IF (partitioned_data%nb_clusters==2) THEN
       PRINT *, 'Ratio difference : ', ratiorij(2)/ratiorii(2)
       IF (ratiomax(2)>=0.6) THEN 
          partitioned_data%nb_clusters=1
       ELSE 
          partitioned_data%nb_clusters=2
       ENDIF
    ENDIF
#if aff
    PRINT *, 'DEBUG : ', proc_id,' : final cluster got : ', partitioned_data%nb_clusters
#endif

    ! Final clustering computing
    IF (partitioned_data%nb_clusters>1) THEN
       CALL apply_spectral_embedding(n, partitioned_data%nb_clusters, proc_id, A, Z,  &
                  nb_info(partitioned_data%nb_clusters), cluster, points_by_clusters,  &
                  ratio, ratiomin(1), ratiorii(1), ratiorij(1), clusters_centers,  &
                  clusters_energies)
       DO i=1,partitioned_data%nb_points
          partitioned_data%points(i)%cluster=cluster(i)
       ENDDO
       DEALLOCATE(cluster)
       DEALLOCATE(points_by_clusters)
       DEALLOCATE(ratiomax)
       DEALLOCATE(clusters_energies)
       DEALLOCATE(ratiomin)
       DEALLOCATE(ratiomoy)
       DEALLOCATE(ratiorii)
       DEALLOCATE(ratiorij)
       DEALLOCATE(A)
       DEALLOCATE(Z)
       IF(solver == 0) DEALLOCATE(A2)
       DEALLOCATE(clusters_centers)
       DEALLOCATE(W)
    ELSE 
#if aff
       PRINT *, proc_id, ' : OK'
#endif
       DO i=1,partitioned_data%nb_points
          partitioned_data%points(i)%cluster=1
       ENDDO
#if aff
       PRINT *, proc_id,' : Cluster'
#endif
    ENDIF

    RETURN
  END SUBROUTINE apply_spectral_clustering


  SUBROUTINE mean_shift(proc_id,nb_clusters_max,nb_clusters_opt,partitioned_data,clust_param)
    INCLUDE 'mpif.h'
    !IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nb_clusters_max
    INTEGER :: nb_clusters_opt
    INTEGER :: proc_id

   TYPE(type_clustering_param) :: clust_param

    !=== IN/OUT ===
    TYPE(type_data) :: partitioned_data
    
    !#### Variables  ####
    INTEGER ::point_num !number of points
    INTEGER ::dim_num !number of dimensions
    INTEGER :: beenVisitedFlag(partitioned_data%nb_points) !track if a point has been seen already
    INTEGER :: clusterVotes(clust_param%nbLimitClust,partitioned_data%nb_points) !number of votes for each point for each cluster
    INTEGER :: i
    INTEGER :: j
    INTEGER :: mergeWith !used to merge clusters
    INTEGER :: myMembers(partitioned_data%nb_points) !1 if the point belongs to the cluster, else 0
    INTEGER :: num
    INTEGER :: cN
    INTEGER :: numClust !the cluster number  
    INTEGER :: numInitPts !number of points to possibly use as initialization points
    INTEGER :: stInd !start point of mean
    INTEGER :: thisClusterVotes(partitioned_data%nb_points) !used to resolve conflicts on cluster membership
    DOUBLE PRECISION :: bandSqG, bandSqC !square of bandwidth
    DOUBLE PRECISION :: clustCent(partitioned_data%dim,clust_param%nbLimitClust) !centers of each cluster
    DOUBLE PRECISION :: myMean(partitioned_data%dim) !mean of this cluster
    DOUBLE PRECISION :: myOldMean(partitioned_data%dim) !old mean computed for this cluster
    DOUBLE PRECISION :: sqDist, sqDistG, sqDistI
    DOUBLE PRECISION :: stopThresh !when mean has converged

    INTEGER :: nbMembers
    INTEGER :: nbIt

    DOUBLE PRECISION :: bandWidthI

!PRINT*, 'MEAN SHIFT IN'
    bandWidthI =  30
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    !INITIALIZE STUFF    	

    point_num = partitioned_data%nb_points
    dim_num = partitioned_data%dim
    print *, '#######################', dim_num
    numClust = 0
    bandSqG = clust_param%bandwidth**2
    bandSqI = bandWidthI
    stopThresh = 1e-1*clust_param%bandwidth
    beenVisitedFlag(:) = 0
    numInitPts = point_num
    clusterVotes(:,:) = 0
    stInd = 0 ! init a verifier
    nbMerg = 0
    
    DO WHILE (numInitPts > 0)

      !take the first point as start of mean
      i = 1;
      do while((i <= point_num).and.(beenVisitedFlag(i) /= 0))
        i = i + 1;
      ENDDO
      stInd = i
      !print *, i

    !  myMean = partitioned_data%point(stInd)%coords !initialize mean to this points location
      DO j = 1, dim_num
        myMean(j) = partitioned_data%points(stInd)%coords(j)
        !print *, 'myMean', j, myMean(j)
      ENDDO

      thisClusterVotes(:) = 0 !used to resolve conflicts on cluster membership
      nbIt = 0
      DO

        nbIt = nbIt + 1
        !print *, 'Construction du Cluster', numClust, nbIt
        myMembers(:) = 0
        nbMembers = 0

        DO i = 1, point_num
          !dist squared from mean to all points still active
          !IF (beenVisitedFlag(i) == 0) THEN
          if(dim_num <= 2) then
            sqDistG = 0.D0
            DO j=1, 2
              sqDistG = sqDistG + (partitioned_data%points(i)%coords(j) - myMean(j))**2
            ENDDO
            IF (sqDistG < bandSqG) THEN
              thisClusterVotes(i) = thisClusterVotes(i) + 1 !add a vote for all the in points belonging to this cluster
              !print *, i, thisClusterVotes(i)
              myMembers(i) = 1 !add any point within bandwidth to the cluster
              beenVisitedFlag(i) = 1 !mark that these points have been visited
              ! print *, 'Members', i
              nbMembers = nbMembers + 1;
            ENDIF
            ! Image 2+3
          else
            sqDistG = 0.D0
            DO j=1, 2
              sqDistG = sqDistG + (partitioned_data%points(i)%coords(j) - myMean(j))**2
            ENDDO
            sqDistI = 0.D0
            DO j=2, dim_num
            sqDistI = sqDistG + (partitioned_data%points(i)%coords(j) - myMean(j))**2
            ENDDO
            IF ((sqDistG < bandSqG).AND.(sqDistI < bandSqI)) THEN
              thisClusterVotes(i) = thisClusterVotes(i) + 1 !add a vote for all the in points belonging to this cluster
              !print *, i, thisClusterVotes(i)
              myMembers(i) = 1 !add any point within bandwidth to the cluster
              beenVisitedFlag(i) = 1 !mark that these points have been visited
              !print *, 'Members', i
              nbMembers = nbMembers + 1;
            ENDIF

            
          end if
          !ENDIF
        ENDDO 
        !print *, 'nbMembers', nbMembers
    
        myOldMean = myMean      
        !compute the new mean
        myMean = 0.D0
        num = 0

        DO i=1, point_num
          IF (myMembers(i) == 1) THEN
            !print *, 'Members', i
            DO j=1, dim_num
              myMean(j) = myMean(j) + partitioned_data%points(i)%coords(j)
            ENDDO
            num = num + 1
          ENDIF
        ENDDO

!PRINT*, 'MEAN SHIFT num  ', num 
        myMean = myMean/num
        !print *, 'myMean', myMean, num

        !compute the distance from myMean to myOldMean
        sqDist = 0

        DO j=1, dim_num
          sqDist = sqDist + (myOldMean(j) - myMean(j))**2
        ENDDO

        !if mean doesn't move much stop this cluster
        !print *, 'MEAN SHIFT sqDist ', sqDist
        
        IF (sqDist < stopThresh) THEN   

          !DO i=1, point_num
          !  IF (myMembers(i) == 1) THEN
          !    print *, 'Members', i
          !  ENDIF
          !ENDDO
          !check for merge posibilities
          mergeWith = 0
          DO cN = 1, numclust
            !compute the distance from possible new clust max to old clust max
            !print *, 'dans boucle cn'
            sqDist = 0
            DO j=1, dim_num
              sqDist = sqDist + (clustCent(j,cN) - myMean(j))**2
            ENDDO
            IF (sqDist < (clust_param%bandwidth/2)**2) THEN
              mergeWith = cN
              EXIT
            ENDIF
          ENDDO

          !print *, 'mergeWith', mergeWith

          IF (mergeWith > 0) THEN !something to merge
            !print *, '####################################################', mergeWith

            clustCent(:,mergeWith) = (myMean+clustCent(:,mergeWith))/2 !mean of centers
            !print*,  'clustCent', clustCent(:,numClust)
            clusterVotes(mergeWith,:) = clusterVotes(mergeWith,:) + thisClusterVotes !add these votes to the merged cluster
            !do i = 1, point_num
            !  print *, i, clusterVotes(mergeWith, i)
            !enddo
          !PRINT*, 'MEAN SHIFT MERGING CLUSTERS'


          ELSE
            numClust = numClust + 1
            clustCent(:,numClust) = myMean
            clusterVotes(numClust,:) = thisClusterVotes
          ENDIF
          EXIT
        ENDIF
      ENDDO
      !do i = 1, point_num
      !  print *, i, clusterVotes(:, i)
      !enddo
      numInitPts = 0
      DO i=1, point_num
        IF (beenVisitedFlag(i)==0) THEN
          numInitPts = numInitPts + 1
        ENDIF
      ENDDO
      !print *, "========================== fin construction ======================="
      !print *, numClust, nbIt, num, numInitPts
      !print*,  'clustCent', clustCent(:,numClust)
      !print *, "==================================================================="
      !if (numClust >= 6) stop
    ENDDO

    DO i=1, point_num

    !PRINT*, 'MEAN SHIFT Affectation a un cluster, point, nb point tot',i, point_num
    !DO j=1, 5
    !	PRINT*, 'MEAN SHIFT clusterVotes', j, clusterVotes(j,i)
    !ENDDO

    !PRINT*, 'MEAN SHIFT MAX LOC', MAXLOC(clusterVotes(:,i), DIM=1)
      partitioned_data%points(i)%cluster = MAXLOC(clusterVotes(:,i), DIM=1)

    ENDDO

    !compute number of clusters
    clust_param%nbLimitClust = numClust
    partitioned_data%nb_clusters=numClust !!!!!!!!!!!!!!!!
    !PRINT*, 'MEAN SHIFT OUT', proc_id
END SUBROUTINE mean_shift


END MODULE module_calcul
