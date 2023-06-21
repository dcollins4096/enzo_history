FUNCTION grid_in_cube, grid_info, cube_center, cube_length, cube_size
; determine whether the specified grid  is at least partially in the
; cube  ;
; return for each grid 1 if it is and 0 if it is not in cube
; min left edge and scaling
  min_left  = $
    [cube_center(0),cube_center(1),cube_center(2)]-cube_length/2. > 0.
  max_right = $
    [cube_center(0),cube_center(1),cube_center(2)]+cube_length/2. < 1.

  scale_up    = min(1./(max_right-min_left))
  b = 0.3*cube_length/FLOAT(cube_size)
  print, 'b:',b
  index_range = grid_info.End_index-grid_info.Start_index+1
  delta_distance = (grid_info.Right_edge - $
                    grid_info.Left_edge)/(index_range)

  suff_res = delta_distance(0,*)
  suff_res = where(delta_distance(0,*) ge b)
  print, suff_res
  print, total(suff_res gt 0), 'grids have sufficient resolution' 
  if total(suff_res gt 0) lt 1 THEN BEGIN
      grid_in_cube = 0
      RETURN, grid_in_cube
  ENDIF
  new_i    = grid_info(suff_res)
;  print, new_i

  Left_edge  = new_i.Left_edge
  Right_edge = new_i.Right_edge

;  print, LEFT_EDGE
  help = ROUND(Left_Edge) < 0
  help = ((Right_edge(0,*)  gt min_left(0) )      and $
          (Right_edge(1,*)  gt min_left(1) )      and $
          (Right_edge(2,*)  gt min_left(2) )      and $
          (Left_edge(0,*)   lt max_right(0))      and $
          (Left_edge(1,*)   lt max_right(1))      and $
          (Left_edge(2,*)   lt max_right(2))          )

  grid_in = new_i(where(help gt 0))
  print, 'there are ', total(help),' grids at least partially in cube'
  return, grid_in
END
