;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/17/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + read data by dataset_name, not sds_num

pro construct_slice_data, grid_info, slice_ori, slice_coord,$
             data_dir,sds_num,slice_size, slice_data, SMOOTH=smooth, SECMIN=SECMIN
; using grid_info read in data (found in data_dir directory) of all
; necessary grids for a slice ; specified by the 3 component vector
; slice_ori. slice_coord has 4 components specifying the left and right
; edge of the slice
; slice_size: is a two component vector giving the x and y size in pixels
; note that slice_ori has to have two elements set to zero since we
; only do orthogonal slices.

  dataset_name=list_str[sds_num]

  num_dim = 3
; rough check of input data
  slice_data = fltarr(slice_size(0), slice_size(1))
  slice_data(*,*) = -1.e30
  axes = (slice_ori gt 1.e-30)
  const_sub =  where(axes eq max(axes))
  const_sub = const_sub(0)
  other_subs = where(axes eq 0)
  
  IF (Fix(TOTAL(axes)) gt 1) then BEGIN
      print, 'construct_slice_data: slice_ori:', $
        slice_ori,' not allowed !'
      STOP
  ENDIF
  IF (N_ELEMENTS(grid_info) lt 1) THEN BEGIN
      print, 'construct_slice_data: no grid information in grid_info:',$
        grid_info
      STOP
  ENDIF
  dim_grid_info = size(grid_info)
  num_of_grids  = dim_grid_info(1)

; min left edge and scaling
  min_left  = [slice_coord(0),slice_coord(1)]
  max_right = [slice_coord(2),slice_coord(3)]
  scale_up    = min(1./(max_right-min_left))

; use only grids that will show up on the image:
  b = 0.7*(max_right-min_left)/FLOAT(slice_size)
; determine grids that are at least partially in the slice
  grid_flag      = grid_in_slice(grid_info, slice_ori, slice_coord)
  index_range    = grid_info.End_index-grid_info.Start_index+1
  delta_distance = (grid_info.Right_edge - $
                    grid_info.Left_edge)/(index_range)

  u_g = grid_info(FIX(where((delta_distance(other_subs,*) ge max(b)) $
                       and (grid_flag gt 0)  ) ))
  dim_grid_info = size(u_g)
  num_of_grids  = dim_grid_info(1)

; initialize minimum of data 
  min_data = 1.e30

  FOR i =0, num_of_grids-1 DO BEGIN
      if ((u_g(i).Left_edge(other_subs(0)) le min_left(0))$
        and  (u_g(i).Left_edge(other_subs(1)) le min_left(1)) $
        and (u_g(i).Right_edge(other_subs(0)) ge max_right(0)) $
        and (u_g(i).Right_edge(other_subs(1)) ge max_right(1))) THEN $
        fully_contained_above = i
  ENDFOR

  print, 'slice is fully contained in grids after index:',  $
    fully_contained_above
  
  index_range = u_g.End_index-u_g.Start_index+1
  delta_distance = (u_g.Right_edge - $
                    u_g.Left_edge)/FLOAT(index_range)
;  print, size(u_g)
;  print,  'max:', min_left, max_right, scale_up
  this_level=0
  FOR i=fully_contained_above, num_of_grids-1 DO BEGIN
;      IF (grid_flag(i) gt 0 ) THEN BEGIN
; check this
      delta_dist   = delta_distance(*,i)
      data_points  = index_range(*,i)  
      slice_ori_l                = fltarr(3)
      slice_ori_l(const_sub)     = slice_ori(const_sub)
      slice_ori_l(other_subs) = [slice_coord(0),slice_coord(1)]
      slice_ori_r                = slice_ori_l 
      slice_ori_r(other_subs) = [slice_coord(2),slice_coord(3)]

; yep check again !!
      i_s = ROUND((slice_ori_l(const_sub) - 0.5*delta_dist(const_sub) -  $
                   u_g(i).Left_edge(const_sub))/delta_dist(const_sub) >0 $
        < (data_points(const_sub)-1))
      i_l = CEIL( (slice_ori_l(other_subs) - $
                   u_g(i).Left_edge(other_subs)) $
                    /delta_dist(other_subs) ) > 0
      i_r   = (FLOOR((slice_ori_r(other_subs) - $
                  u_g(i).Left_edge(other_subs))/ $
                    delta_dist(other_subs) - 1) < $
        (data_points(other_subs)-1) ) > 0
;      print, 'slize_ori:',  slice_ori_l, slice_ori_r
;      print, 'index', i,  i_s, i_r, i_l
;      temp_size  = data_points(other_subs)
;          data_slice = (extract_slice(temp,temp_size(0),temp_size(1), $
;                index(0),index(1),index(2),0.,0.,0.))
      Left_point = FIX(FLOAT(slice_size) *scale_up * $
                         ((u_g(i).Left_edge(other_subs)+ $
                           (i_l)*delta_dist(other_subs)-min_left) > 0.))
      Right_point= Left_point + $
        ROUND(FLOAT(slice_size) *scale_up * $
              float(i_r-i_l+1)*delta_dist(other_subs)) < (slice_size-1)
      sub_slice_size = Right_point-Left_point
;      print, 'left and right',  Left_point  , Right_point
      if ((min(i_r-i_l) ge 2) and (min(sub_slice_size) ge 3)) THEN BEGIN
          read_hdf, data_dir+u_g(i).baryon_file,temp,byname=dataset_name
          s = size(temp)
          if s(1) gt 3 THEN BEGIN
              CASE const_sub OF 
                  0: temp = REFORM(temp(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                  1: temp = REFORM(temp(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                  2: temp = REFORM(temp(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                  ELSE:
              ENDCASE
;              print, 'sub_slice_size', sub_slice_size
              sub_slice= $
                CONGRID(temp,sub_slice_size(0)+1,sub_slice_size(1)+1)
              
              slice_data(Left_point(0):(Right_point(0)),      $
                         Left_point(1):(Right_point(1))) = sub_slice

              min_data = MIN([min_data,MIN(sub_slice)])

              next_level = u_g((i+1)<num_of_grids-1).level 
              if (keyword_set(smooth) and $
                  ((this_level lt next_level) or (i eq (num_of_grids-1)))) $
                THEN BEGIN 
                  this_level = next_level
                  smooth_radius = ROUND(FLOAT(smooth) * $
                                        min(sub_slice_size/(i_r-i_l))) < $
                    ROUND(min([slice_size(0), slice_size(1)])/2)
                  if smooth_radius gt 2 then BEGIN
                      print, 'smoothing radius:', smooth_radius,$
                        ' next level',this_level
                      slice_data = smooth(slice_data,smooth_radius ,$
                                          /edge_truncate)
                  ENDIF
              ENDIF 
          ENDIF
      ENDIF
  ENDFOR

SECMIN = MIN(slice_data(where(slice_data gt min_data)))
slice_data = ALOG10(abs(slice_data) > secmin)*slice_data/(abs(slice_data) > 1.e-30)

END
