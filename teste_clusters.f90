PROGRAM teste_clusters
  USE module_teste_clusters

  IMPLICIT NONE
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  TYPE(type_test), DIMENSION(:), POINTER :: test
  CHARACTER (LEN=30) :: files
  INTEGER :: i
  INTEGER :: nb_tests
  
  !###########################################
  ! INSTRUCTIONS
  !###########################################
  PRINT *
  PRINT *,'-------------------------------------------------------'
  PRINT *,'Programme de test de CLUSTERS et VISUCLUSTERS '
  PRINT *,'-------------------------------------------------------'
  PRINT *
  PRINT *,'teste pour plusieurs configurations differentes :'
  PRINT *,'  - le programme CLUSTERS'
  PRINT *,'  - le programme VISUCLUSTERS pour des sorties paraview'
  PRINT *,'  - le programme VISUCLUSTERS pour des sorties gmsh'
  PRINT *
  PRINT *,'le programme verifie l''existence d''un fichier genere'
  PRINT *,'  en fin d''execution de ces programmes.'
  PRINT *,'il retourne "T" ou "F" selon que ce programme s''est'
  PRINT *,'  bien execute ou non.'
  PRINT *
  PRINT *,'-------------------------------------------------------'
  
  ! Data creation
  PRINT *
  PRINT *,'creation des datas...'
  CALL create_data

  ! Tests declaration
  nb_tests=12
  ALLOCATE(test(nb_tests))

  ! Testing: coord_mono
  test(1)%directory='test_coord_mono'
  test(1)%nb_proc=1
  test(1)%output='verif.coord_mono.txt'
  test(1)%visu_paraview='paraview.coord_mono.txt'
  test(1)%visu_gmsh='gmsh.coord_mono.txt'
  test(1)%file='cible'
  test(1)%data_type='COORD'
  test(1)%partition_type='INTERFACE'
  ALLOCATE(test(1)%partitioning(2))
  test(1)%partitioning(:)=1
  test(1)%thickness=0.

  ! Testing: coord multi interface
  test(2)%directory='test_coord_multi_interface'
  test(2)%nb_proc=5
  test(2)%output='verif.coord_multi_interface.txt'
  test(2)%visu_paraview='paraview.coord_multi_interface.txt'
  test(2)%visu_gmsh='gmsh.coord_multi_interface.txt'
  test(2)%file='cible'
  test(2)%data_type='COORD'
  test(2)%partition_type='INTERFACE'
  ALLOCATE(test(2)%partitioning(2))
  test(2)%partitioning(1)=2
  test(2)%partitioning(2)=2
  test(2)%thickness=0.3

  ! Testing: coord multi recouvre
  test(3)%directory='test_coord_multi_recouvre'
  test(3)%nb_proc=4
  test(3)%output='verif.coord_multi_recouvre.txt'
  test(3)%visu_paraview='paraview.coord_multi_recouvre.txt'
  test(3)%visu_gmsh='gmsh.coord_multi_recouvre.txt'
  test(3)%file='cible'
  test(3)%data_type='COORD'
  test(3)%partition_type='RECOUVREMENT'
  ALLOCATE(test(3)%partitioning(2))
  test(3)%partitioning(1)=2
  test(3)%partitioning(2)=2
  test(3)%thickness=0.3

  ! Testing: image mono
  test(4)%directory='test_image_mono'
  test(4)%nb_proc=1
  test(4)%output='verif.image_mono.txt'
  test(4)%visu_paraview='paraview.image_mono.txt'
  test(4)%visu_gmsh='gmsh.image_mono.txt'
  test(4)%file='image1d'
  test(4)%data_type='IMAGE'
  test(4)%partition_type='INTERFACE'
  ALLOCATE(test(4)%partitioning(2))
  test(4)%partitioning(:)=1
  test(4)%thickness=0.

  ! Testing: image multi interface
  test(5)%directory='test_image_multi_interface'
  test(5)%nb_proc=7
  test(5)%output='verif.image_multi_interface.txt'
  test(5)%visu_paraview='paraview.image_multi_interface.txt'
  test(5)%visu_gmsh='gmsh.image_multi_interface.txt'
  test(5)%file='image1d'
  test(5)%data_type='IMAGE'
  test(5)%partition_type='INTERFACE'
  ALLOCATE(test(5)%partitioning(2))
  test(5)%partitioning(1)=3
  test(5)%partitioning(2)=2
  test(5)%thickness=1.01

  ! Testing: image multi recouvrem
  test(6)%directory='test_image_multi_recouvre'
  test(6)%nb_proc=6
  test(6)%output='verif.image_multi_recouvre.txt'
  test(6)%visu_paraview='paraview.image_multi_recouvre.txt'
  test(6)%visu_gmsh='gmsh.image_multi_recouvre.txt'
  test(6)%file='image1d'
  test(6)%data_type='IMAGE'
  test(6)%partition_type='RECOUVREMENT'
  ALLOCATE(test(6)%partitioning(2))
  test(6)%partitioning(1)=3
  test(6)%partitioning(2)=2
  test(6)%thickness=1.01

  ! Testing: seuil mono
  test(7)%directory='test_seuil_mono'
  test(7)%nb_proc=1
  test(7)%output='verif.seuil_mono.txt'
  test(7)%visu_paraview='paraview.seuil_mono.txt'
  test(7)%visu_gmsh='gmsh.seuil_mono.txt'
  test(7)%file='image1d'
  test(7)%data_type='SEUIL'
  test(7)%partition_type='INTERFACE'
  ALLOCATE(test(7)%partitioning(1))
  test(7)%partitioning(:)=1
  test(7)%thickness=0.

  ! Testing: seuil multi interface
  test(8)%directory='test_seuil_multi_interface'
  test(8)%nb_proc=9
  test(8)%output='verif.seuil_multi_interface.txt'
  test(8)%visu_paraview='paraview.seuil_multi_interface.txt'
  test(8)%visu_gmsh='gmsh.seuil_multi_interface.txt'
  test(8)%file='image1d'
  test(8)%data_type='SEUIL'
  test(8)%partition_type='INTERFACE'
  ALLOCATE(test(8)%partitioning(1))
  test(8)%partitioning(:)=8
  test(8)%thickness=0.01

  ! Testing: seuil multi recouvre
  test(9)%directory='test_seuil_multi_recouvre'
  test(9)%nb_proc=8
  test(9)%output='verif.seuil_multi_recouvre.txt'
  test(9)%visu_paraview='paraview.seuil_multi_recouvre.txt'
  test(9)%visu_gmsh='gmsh.seuil_multi_recouvre.txt'
  test(9)%file='image1d'
  test(9)%data_type='SEUIL'
  test(9)%partition_type='RECOUVREMENT'
  ALLOCATE(test(9)%partitioning(1))
  test(9)%partitioning(:)=8
  test(9)%thickness=0.01

  ! Testing: geom mono
  test(10)%directory='test_geom_mono'
  test(10)%nb_proc=1
  test(10)%output='verif.geom_mono.txt'
  test(10)%visu_paraview='paraview.geom_mono.txt'
  test(10)%visu_gmsh='gmsh.geom_mono.txt'
  test(10)%file='image1d'
  test(10)%data_type='GEOM'
  test(10)%partition_type='INTERFACE'
  ALLOCATE(test(10)%partitioning(3))
  test(10)%partitioning(:)=1
  test(10)%thickness=0.

  ! Testing: geom multi interface
  test(11)%directory='test_geom_multi_interface'
  test(11)%nb_proc=28
  test(11)%output='verif.geom_multi_interface.txt'
  test(11)%visu_paraview='paraview.geom_multi_interface.txt'
  test(11)%visu_gmsh='gmsh.geom_multi_interface.txt'
  test(11)%file='image1d'
  test(11)%data_type='GEOM'
  test(11)%partition_type='INTERFACE'
  ALLOCATE(test(11)%partitioning(3))
  test(11)%partitioning(:)=3
  test(11)%thickness=0.01

  ! Testing: geom multi recouvre
  test(12)%directory='test_geom_multi_recouvre'
  test(12)%nb_proc=27
  test(12)%output='verif.geom_multi_recouvre.txt'
  test(12)%visu_paraview='paraview.geom_multi_recouvre.txt'
  test(12)%visu_gmsh='gmsh.geom_multi_recouvre.txt'
  test(12)%file='image1d'
  test(12)%data_type='GEOM'
  test(12)%partition_type='RECOUVREMENT'
  ALLOCATE(test(12)%partitioning(3))
  test(12)%partitioning(:)=3
  test(12)%thickness=0.01

  ! Launching
  files=''
  DO i=1,nb_tests
     PRINT *
     PRINT *,'TEST : '//test(i)%directory
     IF (files/='t') THEN
        PRINT *,'  > launch test ? [o,n,t]'
        READ *,files
     ENDIF
     IF (files/='n') THEN
        CALL create_executable(test(i))
        CALL create_test(test(i))  
        CALL execute_test(test(i))
     ENDIF
  ENDDO
  PRINT *
  PRINT *,'-------------------------------------------------------'

END PROGRAM teste_clusters
