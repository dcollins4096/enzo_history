;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/16/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + read data by dataset_name, not sds_num

pro find_max_value, option,  u_g, data_dir, sds_num, xyz_vector
; using all grid in the array of grid structure: grid_info
; find the spatial coord[0..1,0..1,0..1]  of the highest value of the
; quantitiy specified by sds_num (see hdf files)
; and return it in xyz_vector

  dataset_name=list_str[sds_num]

  dim_grid_info = size(u_g)
  num_of_grids  = dim_grid_info(1)
  IF (num_of_grids lt 1) THEN BEGIN
      print, 'find_max_value: no grid information in grid_info:', u_g
      STOP
  ENDIF

  IF (option eq 'FIND MAX.from the finest levels') THEN BEGIN
      sec_finest_level = max(u_g.level)   - 1 
      grid_info    = u_g(where(u_g.level ge sec_finest_level))
      END ELSE grid_info = u_g

  dim_grid_info = size(grid_info)
  num_of_grids  = dim_grid_info(1)


  max_value = -1.e30 ; that's small 
  For i=0,num_of_grids-1 DO BEGIN
      read_hdf, data_dir+grid_info(i).baryon_file,temp,byname=dataset_name
      
      data_points = grid_info(i).End_index-grid_info(i).Start_index+1
      delta_dist  = DOUBLE(grid_info(i).Right_edge- $
                     grid_info(i).Left_edge)/DOUBLE(data_points)
      max_value_tg =  max(temp,full_index)
      if (max_value_tg ge max_value) THEN BEGIN
          max_value  =  max_value_tg
          index = indgen(3)
          index(0)   =  full_index mod data_points(0)
          index(1)   =  full_index/data_points(0) mod data_points(1)
          index(2)   =  full_index/ $
            (data_points(0) * data_points(1)) mod data_points(2)
          xyz_vector = DOUBLE(grid_info(i).Left_edge + (index)*delta_dist)
          print, 'xyz_vec, max, maxcheck:', xyz_vector, max_value,$
            temp(index(0),index(1),index(2))
      ENDIF
  ENDFOR

  Print, 'Maximum value is:', max_value
END
;.compile find_max_value.pro
