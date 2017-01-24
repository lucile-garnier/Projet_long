!>Contains useful data structures
MODULE module_visuclusters_structure

  TYPE type_params
     CHARACTER (LEN=30) :: input_file
     INTEGER :: coords
     INTEGER :: dim
     INTEGER :: image_dim
     INTEGER :: image_times
     INTEGER :: is_geom
     INTEGER :: is_image
     INTEGER :: is_interfacing
     INTEGER :: is_overlapping
     INTEGER :: is_threshold
     INTEGER :: nb_clusters
     INTEGER :: nb_points
     INTEGER :: nb_proc
     INTEGER, DIMENSION(:,:), POINTER :: image_ref
     INTEGER, DIMENSION(:), POINTER :: partitioning
     DOUBLE PRECISION, DIMENSION(:), POINTER :: step
  END TYPE type_params

END MODULE module_visuclusters_structure
