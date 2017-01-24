!>Contains structure types required by the different modules
MODULE module_structure

  !#### Datatype for computing ####
  TYPE type_data
     ! Data parameters
     TYPE(type_points), DIMENSION(:), POINTER :: points
     INTEGER :: dim ! Dimension of points
     INTEGER :: nb_clusters ! Number of clusters
     INTEGER :: nb_points ! Number of points

     ! Input parameters : format + processing
     INTEGER :: coords ! Data format classic points
     INTEGER :: is_geom ! Image in mode geom ?
     INTEGER :: is_image ! If image mode activated
     INTEGER :: is_interfacing ! Partitioning + interfacing ?
     INTEGER :: is_overlapping ! Partitioning + overlapping ?
     INTEGER :: is_threshold ! Image in threshold mode ?

     ! Image processing parameters
     INTEGER :: image_dim ! Dimension of images
     INTEGER :: image_times ! Number of "times"
     INTEGER, DIMENSION(:), POINTER :: partitioning   ! Pixel partitionings of picture
     INTEGER, DIMENSION(:,:), POINTER :: image_ref ! Point reference in pixel coordinates
     REAL, DIMENSION(:), POINTER :: step ! Step for geom mode

  END TYPE type_data


  !#### Points description ####
  TYPE type_points
     INTEGER :: cluster
     DOUBLE PRECISION, DIMENSION(:), POINTER :: coords
  END TYPE type_points


  !#### Sub-clusters description ####
  TYPE type_clusters
     INTEGER :: nb
     INTEGER, DIMENSION(:), POINTER :: nb_elements
  END TYPE type_clusters


  !### Kernel parameter ####
  TYPE type_clustering_param
     ! Clustering method id
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: gamma
      DOUBLE PRECISION :: sigma
      INTEGER :: clustering_method_id
      INTEGER :: kernelfunindex
      INTEGER :: nbLimitClust
     !mean shift
     DOUBLE PRECISION :: bandwidth !bandwidth for mean shift

  END TYPE type_clustering_param

END MODULE module_structure
