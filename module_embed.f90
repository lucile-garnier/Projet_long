!>Contains K-means and spectral embedding algorithms
MODULE module_embed
  USE module_structure
CONTAINS



!>Computes the clusters using eigen vector matrix
!!@details The first part of the method performs the following:
!!<ol>
!!<li> Extract the <em>nb_clusters</em> first columns of the eigen vector matrix (corresponding to the highest eigen values)</li>
!!<li> Normalize the matrix and transposes it </li>
!!<li> Apply K-Means on it to find the clusters </li>
!!</ol>
!!Then it operates a quality measurement on the found clusters.
!!It computes the sum of the ratios that indicate if the number of
!!clusters is optimal. The lower this sum is, the best is the number
!!of clusters. This method also compute the number of clusters that
!!have at least one point and that have an internal affinity greater
!!than zero.
!!@note We refer you to the article <em>"On a strategy for Spectral Clustering with parallel computation"</em> for a better understanding.
!! @param[in] A the affinity matrix
!! @param[in] Z the matrix of eigen vectors
!! @param[in] n 
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] proc_id the processus identifier
!! @param[out] clusters_centers the cluster centers
!! @param[out] clusters_energies the cluster energies
!! @param[out] ratio the sum of the ratios between Frobenius norm of the off-diagonal and the diagonal blocks of the normalized affinity matrix
!! @param[out] ratio_moy 
!! @param[out] ratio_rii the sum of the denominator of each ratio
!! @param[out] ratio_rij the sum of the numerators of each ratio
!! @param[out] nb_info the reduced number of clusters (?)
!! @param[out] clusters indicates which cluster each point belongs to
!! @param[out] points_by_clusters the number of points in each cluster
  SUBROUTINE apply_spectral_embedding(n, nb_clusters, proc_id, A, Z, nb_info,  &
                  clusters, points_by_clusters, ratio, ratio_moy, ratio_rii, ratio_rij,  &
                  clusters_centers, clusters_energies)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: n
    INTEGER :: nb_clusters ! nbre de clusters
    INTEGER :: proc_id
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z ! matrice des vecteurs propres
    
    !====  OUT ====
    INTEGER :: nb_info
    INTEGER, DIMENSION(:), POINTER :: clusters ! appartenance des clusters
    INTEGER, DIMENSION(:), POINTER :: points_by_clusters ! nbre de points par cluster
    DOUBLE PRECISION :: ratio ! max des ration de frob sur matrice aff reordonnancee suivant
    DOUBLE PRECISION :: ratio_moy
    DOUBLE PRECISION :: ratio_rii
    DOUBLE PRECISION :: ratio_rij
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: clusters_centers ! centre des clusters
    DOUBLE PRECISION, DIMENSION(:), POINTER :: clusters_energies ! somme des energies par cluster
    
    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: ki
    INTEGER :: kj
    INTEGER :: nb_iter
    INTEGER :: nb_iter_max
    INTEGER :: nb_max
    INTEGER :: ni
    INTEGER :: nj
    INTEGER, DIMENSION(:,:), POINTER :: corresp_cluster
    DOUBLE PRECISION :: ratio_min
    DOUBLE PRECISION :: test
    DOUBLE PRECISION, DIMENSION(:), POINTER :: Z3
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Frob
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z1
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z2
    LOGICAL :: ok
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    ALLOCATE(clusters(n))
    ALLOCATE(clusters_centers(nb_clusters,nb_clusters))
    ALLOCATE(points_by_clusters(nb_clusters))
    ALLOCATE(clusters_energies(nb_clusters))
    ALLOCATE(Z1(n,nb_clusters))
    ALLOCATE(Z2(nb_clusters,n))
    ALLOCATE(Z3(n))
    Z3(:)=0.0

    DO i=1,n
       DO j=1,nb_clusters
          Z1(i,j)=Z(i,j)
          Z3(i)=Z3(i)+Z1(i,j)**2
       ENDDO
    ENDDO

    DO i=1,n
       test=0.0
       DO j=1,nb_clusters
          Z2(j,i)=Z1(i,j)/(sqrt(Z3(i)))
          test=test+Z2(j,i)**2
       ENDDO
    ENDDO

    nb_iter_max=n*n

    CALL apply_kmeans(nb_clusters, nb_clusters, nb_iter_max, n, proc_id,  &
                  Z2, clusters_centers, clusters, nb_iter, points_by_clusters,  &
                  clusters_energies)

    ! Quality measure
    nb_max=0
    DO i=1,nb_clusters
       nb_max=max(nb_max,points_by_clusters(i))
    ENDDO
    ALLOCATE(corresp_cluster(nb_clusters,nb_max))
    corresp_cluster(:,:)=0
    DO i=1,n
       j=clusters(i)
       ok=.FALSE.
       k=1
       DO WHILE(.NOT. ok)
          IF (corresp_cluster(j,k)==0) THEN
             ok=.TRUE.
          ELSE
             k=k+1
          ENDIF
       ENDDO
       corresp_cluster(j,k)=i
    ENDDO
    ALLOCATE(Frob(nb_clusters,nb_clusters))
    Frob(:,:)=0.0 
    DO i=1,nb_clusters
       DO j=1,nb_clusters
          DO ki=1,points_by_clusters(i)
             ni=corresp_cluster(i,ki)
             DO kj=1,points_by_clusters(j)
                nj=corresp_cluster(j,kj)
                Frob(i,j)=Frob(i,j)+A(ni,nj)**2
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(corresp_cluster)
    ratio=0.0
    ratio_min=1.D+16
    ratio_rii=0.0
    ratio_rij=0.0
    nb_info=nb_clusters
    DO i=1,nb_clusters
       IF ((points_by_clusters(i)/=0).AND.(Frob(i,i)/=0)) THEN
          DO j=1,nb_clusters
             IF (i/=j) THEN
                ratio=ratio+Frob(i,j)/Frob(i,i)
                ratio_moy=ratio_moy+Frob(i,j)/Frob(i,i)
                ratio_rij=ratio_rij+Frob(i,j)
                ratio_rii=ratio_rii+Frob(i,i)
                ratio_min=min(ratio_min,Frob(i,j)/Frob(i,i))
             ENDIF
          ENDDO
       ELSE
          nb_info=nb_info-1
       ENDIF
       ratio_rij=ratio_rij*2/(nb_clusters*(nb_clusters-1))
    ENDDO
    DEALLOCATE(Frob)

#if aff
    PRINT *, 'DEBUG : ', proc_id,' : nb_info=', nb_info, ' nb_clusters=', nb_clusters
#endif

    RETURN 
  END SUBROUTINE apply_spectral_embedding



!>Implements K-Means algorithm (required by spectral clustering and Kernel K-Means methods)
!!The algorithm works as follows:
!!<ol>
!!<li>Choose <em>nb_clusters</em> starting points as cluster centers randomly</li>
!!<li>For each point in the data set, find the minimum distance from a cluster center and attach it to the corresponding cluster</li>
!!<li>Compute the density center of each cluster</li>
!!<li>Stop if the density centers are similar to the cluster centers</li>
!!<li>Select the density centers as new cluster centers and start from the beginning</li>
!!</ol>
!! @param[in] points the points
!! @param[in] dim the number of spatial dimensions
!! @param[in] dim the number of spatial dimensions
!! @param[in] dim 
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_iter_max the maximum number of iterations
!! @param[in] nb_points the number of points
!! @param[in] nb_points the number of points
!! @param[in,out] clusters_centers the cluster centers
!! @param proc_id the processus identifier
!! @param[out] clusters_energies the cluster energies
!! @param[out] nb_iter the number of iterations taken
!! @param[out] clusters indicates which cluster each point belongs to
!! @param[out] points_by_clusters the number of points in each cluster
  SUBROUTINE apply_kmeans(dim, nb_clusters, nb_iter_max, nb_points,  &
                  proc_id, points, clusters_centers, clusters, nb_iter,  &
                  points_by_clusters, clusters_energies)

    !*****************************************************************************80
    !
    !! KMEANS_01 applies the K-Means algorithm.
    !
    !  Discussion:
    !
    !    Given a matrix of NB_POINTS observations on DIMENSION variables, the
    !    observations are to be ALLOCATEd to NB_CLUSTERS clusters in such 
    !    a way that the within-cluster sum of squares is minimized.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    18 November 2004
    !
    !  Author:
    !
    !    FORTRAN77 original version by David Sparks
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    David Sparks,
    !    Algorithm AS 58: 
    !    Euclidean Cluster Analysis,
    !    Applied Statistics,
    !    Volume 22, Number 1, 1973, pages 126-130.
    !

    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: dim ! the number of spatial dimensions
    INTEGER :: nb_clusters ! the number of clusters
    INTEGER :: nb_iter_max ! the maximum number of iterations
    INTEGER :: nb_points ! the number of points + TODO reorganisation
    INTEGER :: proc_id ! UNUSED
    DOUBLE PRECISION :: points (dim, nb_points) ! the points

    !=== IN/OUT ===
    DOUBLE PRECISION :: clusters_centers (dim, nb_clusters) ! the cluster centers

    !====  OUT ====
    INTEGER :: clusters (nb_points) ! indicates which cluster each point belongs to
    INTEGER :: nb_iter ! the number of iterations taken
    INTEGER :: points_by_clusters (nb_clusters) ! the number of points in each cluster
    DOUBLE PRECISION :: clusters_energies (nb_clusters) ! the cluster energies
    
    !#### Variables  ####
    INTEGER :: cluster_id (nb_clusters)
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: p
    INTEGER :: stock_population (nb_clusters)
    INTEGER :: swap
    DOUBLE PRECISION :: list_norm (nb_points, nb_clusters)
    DOUBLE PRECISION :: max_value
    DOUBLE PRECISION :: stock_center (dim, nb_clusters)
    DOUBLE PRECISION :: stock_energy (nb_clusters)
    DOUBLE PRECISION :: threshold
    DOUBLE PRECISION :: value
    LOGICAL :: ok
    LOGICAL :: ok2
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################   
    nb_iter = 0
    !
    !  Idiot checks.
    !
    IF ( nb_clusters < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  NB_CLUSTERS < 1.0'
       STOP
    ENDIF

    IF ( dim < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  DIMENSION < 1.0'
       STOP
    ENDIF

    IF ( nb_points < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  NB_POINTS < 1.0'
       STOP
    ENDIF
    !
    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.
    !

    !
    !  Assign one point to each cluster center.
    !
    clusters_centers(:,1) = points(:,1)
    cluster_id(:)=0
    cluster_id(1)=1
    p=2
    threshold=0.4
#if aff
PRINT *, 'DEBUG : searching centers'
#endif
    DO i = 2, nb_clusters
       ok=.FALSE.
       DO WHILE(.NOT. ok)
          max_value=2.0*threshold
          ! Test if the point is already used as center
          ok2=.FALSE.
          DO j=1,i-1
             IF (cluster_id(j)==p) ok2=.TRUE.
          ENDDO
          ! If the point is not a center, test against the threshold
          IF (.NOT. ok2) THEN
             DO j=1,i-1
                value=0.0
                DO k=1,dim
                   value=max(value,abs(clusters_centers(k,j)-points(k,p)))
                ENDDO
                max_value=min(value,max_value)
             ENDDO
             IF (max_value>=threshold) ok=.TRUE.
          ENDIF
         p=p+1

         ! Lower the threshold if not enough centers found
         IF ((p>nb_points).AND.(.NOT. ok)) THEN 
            threshold=0.9*threshold
#if aff
            PRINT *, 'DEBUG : Lower threshold : ', threshold
#endif
            p=1
          ENDIF
       ENDDO
       p=p-1
       clusters_centers(:,i)=points(:,p)
       cluster_id(i)=p
    ENDDO
#if aff
   PRINT *, 'DEBUG : initial centers : ', p
#endif
    nb_iter = 0
    swap=1
    clusters(:)=1
    DO WHILE ((nb_iter<nb_iter_max).AND.(swap/=0))
       nb_iter = nb_iter + 1
       swap=0
       DO i=1,nb_clusters
          stock_energy(i)=clusters_energies(i)
          stock_population(i)=points_by_clusters(i)
          DO j=1,dim
             stock_center(j,i)=clusters_centers(j,i)
          ENDDO
       ENDDO

       ! Computing of the distances
       points_by_clusters(1:nb_clusters) = 1
       list_norm(:,:)=0.0
       DO i=1,nb_points
          DO j=1,nb_clusters
             DO k=1,dim
                list_norm(i,j)=list_norm(i,j)+(points(k,i)-clusters_centers(k,j))**2
             ENDDO
          ENDDO
       ENDDO

       ! Allocation related to the minimum of the distances
       points_by_clusters(:)=0
       DO i=1,nb_points
          DO j=1,nb_clusters
             IF (list_norm(i,j)<list_norm(i,clusters(i))) THEN
                clusters(i)=j
                swap=swap+1
             ENDIF
          ENDDO
          clusters_energies(clusters(i))=clusters_energies(clusters(i))&
               +list_norm(i,clusters(i))
          points_by_clusters(clusters(i))=points_by_clusters(clusters(i))+1
       ENDDO

       ! Update of centers
       clusters_centers(:,:)=0.0
       DO j=1,nb_points
          i=clusters(j) 
          DO k=1,dim
             clusters_centers(k,i)=clusters_centers(k,i)+points(k,j)
          ENDDO
       ENDDO
       DO i=1,nb_clusters
          clusters_centers(:,i)=clusters_centers(:,i)/points_by_clusters(i)
       ENDDO



    ENDDO

    RETURN
  END SUBROUTINE apply_kmeans



END MODULE module_embed
