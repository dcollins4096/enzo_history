pro finest_grid_containing_point, grid_info, point, grid_num, grid_level
; assumes grids are ordered from low to high resolution !
  dim_grid_info = size(grid_info)
  num_of_grids  = dim_grid_info(1)
  
  grid_num   = -1
  grid_level = -1

  ;; modification by mqk (02/20/2004)
  print,point
  for i=num_of_grids-1,0,-1 do begin
      found = (grid_info[i].left_edge  le point) and $
                  (grid_info[i].right_edge ge point)
      
      if min(found) then break
  endfor

  print,found

;  REPEAT BEGIN      
;      i = i -1 
;      print,i
;      FOUND = MIN((grid_info(i).left_edge  lt point) and $
;                  (grid_info(i).right_edge gt point))
;  ENDREP UNTIL FOUND or i lt 0

  IF min(FOUND) THEN BEGIN
      grid_num   = grid_info(i).num
      grid_level = grid_info(i).level
  ENDIF
  
  RETURN
END
  

  
      
