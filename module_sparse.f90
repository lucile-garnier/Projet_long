!>
MODULE module_sparse
  USE module_structure
  USE module_solve
  USE module_embed
CONTAINS

!>Computes the clusters using spectral clustering algorithm using sparsity
!!@see apply_spectral_clustering()
!! @param[in] sigma the affinity parameter
!! @param[in] nb_clusters_max the maximum number of clusters
!! @param[in] nb_clusters_opt the optimal number of clusters
!! @param[in] proc_id the processus identifier
!! @param[in,out] partitioned_data the partitioned data for computing
  SUBROUTINE apply_spectral_clustering_sparse(nb_clusters_max, nb_clusters_opt, proc_id, sigma, partitioned_data)
    IMPLICIT INTEGER(i, j, q)
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nb_clusters_max
    INTEGER :: nb_clusters_opt
    INTEGER :: proc_id
    DOUBLE PRECISION :: sigma

    !=== IN/OUT ===
    TYPE(type_data) :: partitioned_data

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    DOUBLE PRECISION :: factor
    DOUBLE PRECISION :: norm
    DOUBLE PRECISION :: ratio
    DOUBLE PRECISION :: ratio1
    DOUBLE PRECISION :: ratio2
    DOUBLE PRECISION :: t1
    DOUBLE PRECISION :: t2
    DOUBLE PRECISION :: t_cons_a
    DOUBLE PRECISION :: t_cons_vp
    DOUBLE PRECISION :: threshold
    DOUBLE PRECISION :: threshold_rij
    DOUBLE PRECISION, DIMENSION(:), POINTER :: AS
    DOUBLE PRECISION, DIMENSION(:), POINTER :: clusters_energies
    DOUBLE PRECISION, DIMENSION(:), POINTER :: D
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratio_moy
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratio_rii
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratio_rij
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomin
    DOUBLE PRECISION, DIMENSION(:), POINTER :: W
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: clusters_centers
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z
    INTEGER :: k
    INTEGER :: l
    INTEGER :: n
    INTEGER :: nb
    INTEGER :: nb_clusters
    INTEGER :: nbproc
    INTEGER :: nnz
    INTEGER :: nnz2
    INTEGER, DIMENSION(:), POINTER :: clusters
    INTEGER, DIMENSION(:), POINTER :: IAS
    INTEGER, DIMENSION(:), POINTER :: JAS
    INTEGER, DIMENSION(:), POINTER :: nb_info
    INTEGER, DIMENSION(:), POINTER :: points_by_clusters

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    ! Matrix creation
#if aff
    PRINT *, 'DEBUG : ', proc_id, ' : value of sigma : ', sigma
#endif
    n=partitioned_data%nb_points

    ! Beginning of sparsification
    nnz = 0
    ! Arbitrary threshold value 
    ! TODO : mettre la valeur du facteur dans le fichier param
    factor = 3.0
    threshold = factor*sigma

    t1 = MPI_WTIME()
    DO i=1,n-1  ! bound ?
       DO j=i+1,n ! bound ?

          norm=0.0

          DO k=1,partitioned_data%dim
             norm=norm+(partitioned_data%points(i)%coords(k)-partitioned_data%points(j)%coords(k))**2
          ENDDO

          IF(sqrt(norm) <= threshold) THEN
            nnz = nnz + 1
          ENDIF

       ENDDO
    ENDDO

    t2 = MPI_WTIME()
    t_cons_a = t2 - t1
    PRINT *, proc_id, ' : t_cons A : ', t_cons_a

    t1 = MPI_WTIME()
    nnz2 = nnz*2

    ALLOCATE(AS(nnz2))
    ALLOCATE(IAS(nnz2))
    ALLOCATE(JAS(nnz2))
    l = 1
    DO i=1,n-1
       DO j=i+1,n
          norm=0.0
          DO k=1,partitioned_data%dim
             norm=norm+(partitioned_data%points(i)%coords(k)-partitioned_data%points(j)%coords(k))**2
          ENDDO
          value=exp(-norm/sigma)
          ! kepp if value <= threshold
          ! (if we want to keep it all, do comment line IF, ENDIF)
          IF(sqrt(norm) <= threshold) THEN
            AS(l) = value
            IAS(l) = i
            JAS(l) = j
            l = l+1
            !------
            AS(l) = value
            IAS(l) = j
            JAS(l) = i
            l = l+1
          ENDIF
       ENDDO
    ENDDO
    WRITE(*,*) '========== factor, n*n nnz2 = ', factor, n*n, nnz2

  ALLOCATE(D(n))
  D(:)=0.0
  DO l=1, nnz2
    D(IAS(l)) = D(IAS(l)) + AS(l)
  ENDDO

  DO l=1, nnz2
    AS(l)=AS(l)/D(IAS(l))
  ENDDO

  DEALLOCATE(D)

    ! nb and nb_clusters_max same value ?
    nb = 2*nb_clusters_max

    t1 = MPI_WTIME()
    CALL solve_arpack(n, nb, nnz2, IAS, JAS, AS, W, Z)
    PRINT *, "---------- W -------------"
    DO i=1,nb
       PRINT *, 'Pure eigen values Arpack : ', i, W(i)
    ENDDO
    

    ! Reorder eigen values... QUESTION: essential with arpack ?
    DO i=1,nb-1
       DO j=i+1,nb
          IF (W(i)<W(j)) THEN
             value=W(i)
             W(i)=W(j)
             W(j)=value
             DO k=1,n
                value=Z(k,i)
                Z(k,i)=Z(k,j)
                Z(k,j)=value
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    DO i=1,nb
       PRINT *, 'Reordered eigen values Arpack : ', i, W(i)
    ENDDO

    ! Test spectral embedding with different nb_clusters   
    ! Spectral embedding

    IF ((nb_clusters_opt==0).AND.(n>2)) THEN
       ! Searching the best partitioning
       ALLOCATE(ratiomax(nb_clusters_max))
       ratiomax(:)=0
       ALLOCATE(ratiomin(nb_clusters_max))
       ratiomin(:)=0
       ALLOCATE(ratio_moy(nb_clusters_max))
       ratio_moy(:)=0
       ALLOCATE(ratio_rii(nb_clusters_max))
       ratio_rii(:)=0
       ALLOCATE(ratio_rij(nb_clusters_max))
       ratio_rij(:)=0

       ALLOCATE(nb_info(nb_clusters_max))
       nb_info(:)=0

       DO nb_clusters = 2 ,min(n,nb_clusters_max)

          ALLOCATE(clusters(n))
          clusters(:)=0.0
          ALLOCATE(clusters_centers(nb_clusters,nb_clusters))
          clusters_centers(:,:)=0.0
          ALLOCATE(points_by_clusters(nb_clusters))
          points_by_clusters(:)=0.0
          ALLOCATE(clusters_energies(nb_clusters))
          clusters_energies(:)=0.0

          CALL apply_spectral_embedding_sparse(n, nb_clusters, nnz2, proc_id, IAS, JAS, AS, Z, clusters, points_by_clusters, ratiomax(nb_clusters), ratio_moy(nb_clusters), ratio_rii(nb_clusters), ratio_rij(nb_clusters), clusters_energies, clusters_centers, nb_info(nb_clusters))

          DEALLOCATE(clusters)
          DEALLOCATE(clusters_centers)
          DEALLOCATE(clusters_energies)
          DEALLOCATE(points_by_clusters)
       ENDDO


#if aff
PRINT *, 'DEBUG : Frobenius ratio'
#endif
       ! Ratio of Frobenius norm
       ratio=ratiomax(nb_clusters_max)
       partitioned_data%nb_clusters=nb_clusters_max
       ratio1=0.0
       ratio2=1e+10

       DO i=2,nb_clusters_max
          IF ((proc_id==0).AND.(nbproc>1)) THEN 
             threshold_rij=1e-1
          ELSE
             threshold_rij=1e-4
          ENDIF

          IF ((ratio_rii(i)>=0.95*ratio1).AND.(ratio_rij(i)-ratio2<=threshold_rij)) THEN  
             partitioned_data%nb_clusters=i
             ratio1=ratio_rii(i)
             ratio2=ratio_rij(i)
          ENDIF
       ENDDO

    ELSEIF ((nb_clusters_opt==1).AND.(n>nb_clusters_opt)) THEN
       ! Test with an imposed cluster
       ALLOCATE(nb_info(nb_clusters_opt))
       nb_info(:) = 0
       ALLOCATE(ratiomin(1))
       ratiomin(:) = 0.0
       partitioned_data%nb_clusters = nb_clusters_opt
    ELSE
       ! Case of a domain with less points than nb_clusters_opt or only one point
       ALLOCATE(nb_info(n))
       nb_info(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       partitioned_data%nb_clusters=n
       ALLOCATE(ratiomax(n))
       ratiomax(:)=0
       ALLOCATE(ratio_moy(n))
       ratio_moy(:)=0
       ALLOCATE(ratiomin(n))
       ratiomin(:)=0
       ALLOCATE(ratio_rii(n))
       ratio_rii(:)=0
       ALLOCATE(ratio_rij(n))
       ratio_rij(:)=0
    ENDIF
    ! Case with nb_clusters==1
    IF (partitioned_data%nb_clusters==2) THEN
       PRINT *, 'Ratio difference : ', ratio_rij(2)/ratio_rii(2)
       IF (ratiomax(2)>=0.6) THEN 
          partitioned_data%nb_clusters=1
       ELSE 
          partitioned_data%nb_clusters=2
       ENDIF
    ENDIF
#if aff
    PRINT *, 'DEBUG : ', proc_id,' : final clusters got : ', partitioned_data%nb_clusters
#endif

    ! Computing final clustering
    IF (partitioned_data%nb_clusters>1) THEN

       CALL apply_spectral_embedding_sparse(n, partitioned_data%nb_clusters, nnz2, proc_id, IAS, JAS, AS, Z, clusters, points_by_clusters, ratio, ratiomin(1), ratio_rii(1), ratio_rij(1), clusters_energies, clusters_centers, nb_info(partitioned_data%nb_clusters))

       DO i=1,partitioned_data%nb_points
          partitioned_data%points(i)%clusters=clusters(i)
       ENDDO

       DEALLOCATE(clusters)
       DEALLOCATE(points_by_clusters)
       DEALLOCATE(ratiomax)
       DEALLOCATE(clusters_energies)
       DEALLOCATE(ratiomin)
       DEALLOCATE(ratio_moy)
       DEALLOCATE(ratio_rii)
       DEALLOCATE(ratio_rij)
       DEALLOCATE(clusters_centers)

    ELSE 
#if aff
       PRINT *, 'DEBUG : ', proc_id, ' : OK'
#endif
       DO i=1,partitioned_data%nb_points
          partitioned_data%points(i)%clusters=1
       ENDDO
#if aff
       PRINT *, 'DEBUG : ', proc_id, ' : cluster'
#endif
    ENDIF

    !deallocations
    DEALLOCATE(AS)
    DEALLOCATE(IAS)
    DEALLOCATE(JAS)
    DEALLOCATE(W)
    DEALLOCATE(Z)

    RETURN
  END SUBROUTINE apply_spectral_clustering_sparse

!>Computes the ideal number of clusters using sparsity
!!@see apply_spectral_embedding()
!! @param[in] Z the matrix of eigen vectors
!! @param[in] AS the affinity sparse matrix
!! @param[in] n 
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nnz the number of non-zero coefficients
!! @param[in] proc_id the processus identifier
!! @param[in] IAS the row indices of the affinity matrix coefficients
!! @param[in] JAS the column indices of the affinity matrix coefficients
!! @param[out] clusters_centers the cluster centers
!! @param[out] clusters_energies the cluster energies
!! @param[out] ratio the sum of the ratios between Frobenius norm of the off-diagonal and the diagonal blocks of the normalized affinity matrix
!! @param[out] ratio_moy 
!! @param[out] ratio_rii the sum of the denominator of each ratio
!! @param[out] ratio_rij the sum of the numerators of each ratio
!! @param[out] nb_info the reduced number of clusters (?)
!! @param[out] clusters indicates which cluster each point belongs to
!! @param[out] points_by_clusters the number of points in each cluster
    SUBROUTINE apply_spectral_embedding_sparse(n, nb_clusters, nnz, proc_id, IAS, JAS, AS, Z, clusters, points_by_clusters, ratio, ratio_moy, ratio_rii, ratio_rij, clusters_energies, clusters_centers, nb_info)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: n
    INTEGER :: nb_clusters      
    INTEGER :: nnz
    INTEGER :: proc_id
    INTEGER, DIMENSION(:), POINTER :: IAS
    INTEGER, DIMENSION(:), POINTER :: JAS
    DOUBLE PRECISION, DIMENSION(:), POINTER:: AS
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z

    !====  OUT ====
    INTEGER, DIMENSION(:), POINTER :: clusters
    INTEGER, DIMENSION(:), POINTER :: points_by_clusters
    DOUBLE PRECISION :: ratio
    DOUBLE PRECISION :: ratio_moy
    DOUBLE PRECISION :: ratio_rii
    DOUBLE PRECISION :: ratio_rij
    DOUBLE PRECISION, DIMENSION(:), POINTER :: clusters_energies
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: clusters_centers
    INTEGER :: nb_info

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: it_max
    INTEGER :: it_num
    INTEGER :: j
    INTEGER :: k
    INTEGER :: ki
    INTEGER :: kj
    INTEGER :: l
    INTEGER :: nb_max
    INTEGER :: ni
    INTEGER :: nj
    INTEGER :: num1
    INTEGER :: num2
    INTEGER :: ok
    INTEGER, DIMENSION(:,:), POINTER :: matchings
    DOUBLE PRECISION :: ratiomin
    DOUBLE PRECISION :: test
    DOUBLE PRECISION, DIMENSION(:), POINTER :: Z3
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Frob
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z1
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z2
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

    PRINT *, '************ sp_spectral_embedding *************'
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

    PRINT *, proc_id,' : kmeans method'

    it_max=n*n !1000.0

    CALL apply_kmeans(nb_clusters, nb_clusters, it_max, n, proc_id, Z2, clusters_centers, clusters, it_num, points_by_clusters, clusters_energies)

    ! Quality measure

    nb_max=0
    DO i=1,nb_clusters
       nb_max=max(nb_max,points_by_clusters(i))
    ENDDO
    PRINT *, 'points_by_clusters : ', points_by_clusters
    ALLOCATE(matchings(nb_clusters,nb_max))
    matchings(:,:)=0
    DO i=1,n
       j=clusters(i)
       ok=0
       k=1
       DO WHILE(ok==0)
          IF (matchings(j,k)==0) THEN
             ok=1
          ELSE
             k=k+1
          ENDIF
       ENDDO
       matchings(j,k)=i
    ENDDO


    ! Beginning of sparsification
    ALLOCATE(Frob(nb_clusters,nb_clusters))
    Frob(:,:)=0.0
    DO i=1, nnz
      num1 = clusters(IAS(i))
      num2 = clusters(JAS(i))
      Frob(num1, num2) = Frob(num1, num2) + AS(i)**2
    ENDDO
    ! End of sparsification


    ! Beginning of sparsification
    ratio=0.0
    ratiomin=1.D+16
    ratio_rii=0.0
    ratio_rij=0.0
    ratio_moy = 0.0
    nb_info=nb_clusters
    DO i=1,nb_clusters
       IF ((points_by_clusters(i)/=0).AND.(Frob(i,i)/=0)) THEN
          DO j=1,nb_clusters
             IF (i/=j) THEN
                ratio=ratio+Frob(i,j)/Frob(i,i)
                ratio_moy=ratio_moy+Frob(i,j)/Frob(i,i)
                ratio_rij=ratio_rij+Frob(i,j)
                ratio_rii=ratio_rii+Frob(i,i)
                ratiomin=min(ratiomin,Frob(i,j)/Frob(i,i))
             ENDIF
          ENDDO
       ELSE
          nb_info=nb_info-1
       ENDIF
       ratio_rij=ratio_rij*2/(nb_clusters*(nb_clusters-1))
       ratio_moy=ratio_moy*2/(nb_clusters*(nb_clusters-1))
       ratio_rii=ratio_rii!/nb_clusters
    ENDDO

    PRINT *, "============= ratio ================", ratio_moy, ratio_rij

    DEALLOCATE(Frob)
    ! End of sparsification

#if aff
    PRINT *, proc_id,' : nb_info=', nb_info, ' nb_clusters=', nb_clusters
#endif

    RETURN 
    END SUBROUTINE apply_spectral_embedding_sparse

!>Computes the matrix vector product using sparsity
!! @param[in] A the sparse matrix
!! @param[in] X the input vector
!! @param[in] n 
!! @param[in] nnz the number of non-zero coefficients
!! @param[in] IA the row indices of the matrix coefficients
!! @param[in] JA the column indices of the matrix coefficients
!! @param[out] Y the resulting vector
    SUBROUTINE compute_matvec_prod(IA, JA, n, nnz, X, A, Y)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER, DIMENSION(nnz) :: IA
    INTEGER, DIMENSION(nnz) :: JA
    INTEGER, :: n
    INTEGER, :: nnz
    DOUBLE PRECISION, DIMENSION(n) :: X
    DOUBLE PRECISION, DIMENSION(nnz) :: A

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(n) :: Y

    !#### Variables  ####
    INTEGER :: l

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    Y(:) = dfloat(0)
    DO l = 1, nnz
        Y(IA(l)) = Y(IA(l)) + A(l)*X(JA(l))
    ENDDO
    RETURN
    END SUBROUTINE compute_matvec_prod



!>
!! @param[in] A the affinity sparse matrix
!! @param[in] dim the number of spatial dimensions
!! @param[in] dim the number of spatial dimensions
!! @param[in] dim 
!! @param[in] nb_clusters_max the maximum number of clusters
!! @param[in] nnz the number of non-zero coefficients
!! @param[in] IA the row indices of the affinity matrix coefficients
!! @param[in] JA the column indices of the affinity matrix coefficients
!! @param[out] Z the matrix of eigen vectors
!! @param[out] W 
  SUBROUTINE solve_arpack(dim, nb_clusters_max, nnz, IA, JA, A, W(:), Z(:,:))
!     %-------------------------------------------------%
!     | The following INCLUDE statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting mnaupd = 1.                    |
!     %-------------------------------------------------%
    INCLUDE 'debug.h'
    INTRINSIC abs
    EXTERNAL dlapy2
    EXTERNAL dnrm2
    EXTERNAL daxpy 

    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: dim
    INTEGER :: nb_clusters_max
    INTEGER :: nnz
    INTEGER, DIMENSION(:) :: IA
    INTEGER, DIMENSION(:) :: JA
    DOUBLE PRECISION, DIMENSION(:) :: A

    !====  OUT ====
    DOUBLE PRECISION, POINTER :: W(:)
    DOUBLE PRECISION, POINTER :: Z(:,:)

    !#### Variables  ####
    CHARACTER :: bmat*1
    CHARACTER :: which*2
    INTEGER :: i
    INTEGER :: ido
    INTEGER :: ierr
    INTEGER :: info
    INTEGER :: iparam(11)
    INTEGER :: ipntr(14)
    INTEGER :: ishfts
    INTEGER :: j
    INTEGER :: ldv
    INTEGER :: lworkl
    INTEGER :: maxitr
    INTEGER :: maxn
    INTEGER :: maxncv
    INTEGER :: maxnev
    INTEGER :: mode1
    INTEGER :: n
    INTEGER :: nbite
    INTEGER :: nconv
    INTEGER :: ncv
    INTEGER :: nev
    INTEGER :: nx
    DOUBLE PRECISION :: dlapy2
    DOUBLE PRECISION :: dnrm2
    DOUBLE PRECISION :: sigmai
    DOUBLE PRECISION :: sigmar
    DOUBLE PRECISION :: zero
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ax
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: resid
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workd
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workl
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: d
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: v
    LOGICAL :: first
    LOGICAL :: rvec
    LOGICAL, DIMENSION(:), ALLOCATABLE :: array_select

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    PARAMETER (zero = 0.0D+0)

! link between lengths
      maxn = dim
      ldv = maxn
      maxnev = nb_clusters_max
      maxncv = 2*maxnev + 1

! memory allocation
      ALLOCATE(SELECT(maxn))
      ALLOCATE(ax(maxn), resid(maxn), workd(3*maxn), &
               workev(3*maxncv), workl(3*maxncv*maxncv+6*maxncv))
      ALLOCATE(d(maxncv, 3), v(ldv, maxncv))

      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 1
      mnaup2 = 0
      mneigh = 0
      mneupd = 0
!
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
      n     = dim

!     %-----------------------------------------------%
!     |                                               |
!     | Specifications for ARPACK usage are set       |
!     | below:                                        |
!     |                                               |
!     |    1) NEV = 4  asks for 4 eigenvalues to be   |
!     |       computed.                               |
!     |                                               |
!     |    2) NCV = 20 sets the length of the Arnoldi |
!     |       factorization.                          |
!     |                                               |
!     |    3) This is a standard problem.             |
!     |         (indicated by bmat  = 'I')            |
!     |                                               |
!     |    4) Ask for the NEV eigenvalues of          |
!     |       largest magnitude.                      |
!     |         (indicated by which = 'LM')           |
!     |       See DOcumentation in DNAUPD for the     |
!     |       other options SM, LR, SR, LI, SI.       |
!     |                                               |
!     | Note: NEV and NCV must satisfy the following  |
!     | conditions:                                   |
!     |              NEV <= MAXNEV                    |
!     |          NEV + 2 <= NCV <= MAXNCV             |
!     |                                               |
!     %-----------------------------------------------%
!
      nev   = maxnev
      ncv   = maxncv
      bmat  = 'I'
      which = 'LR'
!
      IF ( n .GT. maxn ) THEN
         PRINT *, ' ERROR with _NSIMP: N is greater than MAXN '
         GOTO 9000
      ELSEIF ( nev .GT. maxnev ) THEN
         PRINT *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         GOTO 9000
      ELSEIF ( ncv .GT. maxncv ) THEN
         PRINT *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         GOTO 9000
      ENDIF
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling DNAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .LE. 0,  THEN TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION PARAMETER         |
!     |      used to specify actions to be taken on return  |
!     |      from DNAUPD. (see usage below)                 |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to DNAUPD.                                |
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     |
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID).  |
!     |                                                     |
!     | The work array WORKL is used in DNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+6*ncv 
      tol    = 1.D-6
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting IPARAM(1) = 1).             |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the DOcumentation in |
!     | DNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode1 = 1
!
      iparam(1) = ishfts
!
      iparam(3) = maxitr
!
      iparam(7) = mode1
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
      nbite = 1
 10   CONTINUE
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         CALL dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                       v, ldv, iparam, ipntr, workd, workl, lworkl, & 
                       info )
!
         IF (ido .EQ. -1 .OR. ido .EQ. 1) THEN
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- Op*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and RETURN the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
            CALL compute_matvec_prod(IA, JA, dim, nnz, workd(ipntr(1)), A, workd(ipntr(2)))

            nbite = nbite + 1
!
!           %-----------------------------------------%
!           | L O O P   B A C K to CALL DNAUPD again. |
!           %-----------------------------------------%
!
            GOTO 10
!
         ENDIF
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      IF ( info .LT. 0 ) THEN
!
!        %--------------------------%
!        | Error message, check the |
!        | DOcumentation in DNAUPD. |
!        %--------------------------%
!
         PRINT *, ' '
         PRINT *, ' Error with _naupd, info = ',info
         PRINT *, ' Check the DOcumentation of _naupd'
         PRINT *, ' '
!
      ELSE 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .TRUE.)    |
!        |                                           |
!        | The routine DNEUPD now CALLed to DO this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1,)                                   |
!        |                                           |
!        %-------------------------------------------%
!
         rvec = .TRUE.
!
         CALL dneupd ( rvec, 'A', array_select, d, d(1,2), v, ldv, &
              sigmar, sigmai, workev, bmat, n, which, nev, tol, & 
              resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
              lworkl, ierr )
!
!        %------------------------------------------------%
!        | The REAL parts of the eigenvalues are returned |
!        | in the first column of the two dimensional     |
!        | array D, and the IMAGINARY part are returned   |
!        | in the second column of D.  The corresponding  |
!        | eigenvectors are returned in the first         |
!        | NCONV (= IPARAM(5)) columns of the two         |
!        | dimensional array V if requested.  Otherwise,  |
!        | an orthogonal basis for the invariant subspace |
!        | corresponding to the eigenvalues in D is       |
!        | returned in V.                                 |
!        %------------------------------------------------%
!
         IF ( ierr .NE. 0) THEN
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
!
            PRINT *, ' '
            PRINT *, ' Error with _neupd, info = ', ierr
            PRINT *, ' Check the documentation of _neupd. '
            PRINT *, ' '
!
         ELSE
!
            first = .TRUE.
            nconv =  iparam(5)
            DO 20 j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (IPARAM(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               IF (d(j,2) .EQ. zero)  THEN
!
!                 %--------------------%
!                 | Ritz value is REAL |
!                 %--------------------%
!
                  !CALL av(nx, v(1,j), ax)
                  CALL compute_matvec_prod(IA, JA, dim, nnz, v(1, j), A, ax)
                  CALL daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
!
               ELSEIF (first) THEN
!
!                 %------------------------%
!                 | Ritz value is COMPLEX. |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
                  !CALL av(nx, v(1,j), ax)
                  CALL compute_matvec_prod(IA, JA, dim, nnz, v(1, j), A, ax)
                  CALL daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  CALL daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  !CALL av(nx, v(1,j+1), ax)
                  CALL compute_matvec_prod(IA, JA, dim, nnz, v(1, j+1), A, ax)
                  CALL daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  CALL daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .FALSE.
               ELSE
                  first = .TRUE.
               ENDIF
!
 20         CONTINUE
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            CALL dmout(6, nconv, 3, d, maxncv, -6, &
                 'Ritz values (Real, Imag) and residual residuals')
         ENDIF
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         IF ( info .EQ. 1) THEN
             PRINT *, ' '
             PRINT *, ' Maximum number of iterations reached.'
             PRINT *, ' '
         ELSEIF ( info .EQ. 3) THEN
             PRINT *, ' ' 
             PRINT *, ' No shifts could be applied during IMPLICIT &
                        Arnoldi update, try increasing NCV.'
             PRINT *, ' '
         ENDIF      
!
         PRINT *, ' '
         PRINT *, ' _NSIMP '
         PRINT *, ' ====== '
         PRINT *, ' '
         PRINT *, ' Size of the matrix is ', n
         PRINT *, ' The number of Ritz values requested is ', nev
         PRINT *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         PRINT *, ' What portion of the spectrum: ', which
         PRINT *, ' The number of converged Ritz values is ', &
                    nconv 
         PRINT *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         PRINT *, ' The number of OP*x is ', iparam(9)
         PRINT *, ' The convergence criterion is ', tol
         PRINT *, ' '
!
      ENDIF
!
!     %---------------------------%
!     | Done with program dnsimp. |
!     %---------------------------%
!
 9000 CONTINUE

      ALLOCATE(W(nb_clusters_max))
      ALLOCATE(Z(n, nb_clusters_max))

      W(1:nb_clusters_max) = d(1:nb_clusters_max, 1)

      DO i = 1, nb_clusters_max
        Z(:,i) = v(:,i)
      ENDDO

      DEALLOCATE(array_select)
      DEALLOCATE(ax, resid, workd, workev, workl)
      DEALLOCATE(d, v)

  END SUBROUTINE solve_arpack

  END MODULE module_sparse
