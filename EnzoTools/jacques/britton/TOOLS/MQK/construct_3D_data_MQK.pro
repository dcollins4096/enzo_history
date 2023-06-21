;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/17/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + takes dataset_name as an argument, instead of sds_num.
;; + load data by dataset_name instead of sds_num (see read_grid_hdf.pro).
;;

pro construct_3D_data, grid_info, cube_center, cube_length,$
             data_dir, dataset_name, cube_size, cube_data, smooth=smooth
; using grid_info read in data (found in data_dir directory) of all
; necessary grids for specified cube ;
; cube_length: side length [0..1]
; cube_center: [x,y,z]

  help,dataset_name

  num_dim = 3
; rough check of input data
  cube_data = dblarr(cube_size, cube_size, cube_size)
  cube_data(*,*,*) = 0.
  
  if (n_elements(grid_info) lt 1) then begin
      print, 'construct_3d_data: no grid information in grid_info:',$
        grid_info
      stop
  endif
  dim_grid_info = size(grid_info)
  num_of_grids  = dim_grid_info[1]

; min left edge and scaling
  min_left  = $
    [cube_center(0),cube_center(1),cube_center(2)]-cube_length/2. > 0.
  max_right = $
    [cube_center(0),cube_center(1),cube_center(2)]+cube_length/2. < 1.

  scale_up    = min(1./(max_right-min_left))

; determine grids that have sufficient res and are at least partially
; in the cube 
  u_g = grid_in_cube(grid_info, cube_center, cube_length, cube_size)
  
  index_range = u_g.End_index-u_g.Start_index+1
  delta_distance = (u_g.Right_edge - $
                    u_g.Left_edge)/DOUBLE(index_range)
  dim_grid_info = size(u_g)
  num_of_grids  = dim_grid_info(1)
;   print, size(u_g)
;  print,  'max:', min_left, max_right, scale_up
   this_level = 0
  for i=0,num_of_grids-1 do begin
; check this
      delta_dist   = delta_distance(*,i)
      data_points  = index_range(*,i)  

      i_l = CEIL((min_left - u_g(i).Left_edge)/delta_dist >0 $
                  < (data_points-1))
      i_r   = FLOOR(((max_right -u_g(i).Left_edge)/ $
                    delta_dist - 1) < (data_points-1) ) > 0

;      print, 'index', i, i_r, i_l

      Left_point = FIX(DOUBLE(cube_size) *scale_up * $
                         ((u_g(i).Left_edge+ $
                           (i_l)*delta_dist)-min_left) > 0.)
      Right_point= Left_point + $
        FIX(DOUBLE(cube_size) *scale_up * $
              double(i_r-i_l+1)*delta_dist) < (cube_size-1)
      sub_cube_size = Right_point-Left_point
;      print, 'left and right',  Left_point  , Right_point
      if ((min(i_r-i_l) ge 3) and (min(sub_cube_size) ge 3)) THEN BEGIN
          read_hdf, data_dir+u_g(i).baryon_file,temp,byname=dataset_name
          temp_s = size(temp) 
; chedck whether data is really what promised by the hierarchy file !          
          if (min(i_r-i_l - [temp_s(1),temp_s(2),temp_s(3)]) lt 1) then begin
              temp = REFORM(temp(i_l(0):i_r(0),i_l(1):i_r(1),i_l(2):i_r(2)))
          
;          print, 'sub_cube_size', sub_cube_size
              sub_cube= CONGRID(temp, $
                sub_cube_size(0)+1,sub_cube_size(1)+1,sub_cube_size(2)+1)
              
              cube_data(Left_point(0):Right_point(0),      $
                        Left_point(1):Right_point(1),      $
                        Left_point(2):Right_point(2)       ) = sub_cube
          endif
      endif
      next_level = u_g((i+1)<num_of_grids-1).level 
      if (keyword_set(smooth) and (this_level lt next_level)) THEN BEGIN 
          this_level = next_level
          smooth_radius = ROUND(DOUBLE(smooth) * $
                                min(sub_cube_size/(i_r-i_l))) < $
            ROUND(cube_size/2)
          if smooth_radius gt 2 then BEGIN
              print, 'smoothing radius:', smooth_radius,$
                ' next level',this_level
              cube_data = smooth(cube_data,smooth_radius ,$
                                  /edge_truncate)
          ENDIF
      ENDIF
  ENDFOR


END
