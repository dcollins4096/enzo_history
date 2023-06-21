pro order_grids, grid_info
; order grids from low to high resolution:
; and insert levels numbers
  IF (N_ELEMENTS(grid_info) lt 1) THEN BEGIN
      print, 'order_grids: no grid information in grid_info:',$
        grid_info
      RETURN
  ENDIF

  index_range = grid_info.End_index-grid_info.Start_index+1
  delta_distance = DOUBLE(grid_info.Right_edge - $
                    grid_info.Left_edge)/DOUBLE(index_range)
  order = (1.D/delta_distance(0,*))
; this is hardwired for now:
  refined_by = 2
  levels = ROUND((ALOG(DOUBLE(order))/ALOG(DOUBLE(refined_by))))
; normalize:
  levels = levels - Min(levels)
  grid_info = grid_info(sort(levels))
  grid_info(*).level = levels(sort(levels))
  num_of_levels = Max(levels)
  print, 'There are ', N_ELEMENTS(levels), ' grids with ', num_of_levels, $
    ' levels of refinement!'

;  print, 'histogram:'
;  print,'  0         10        20        30        40        50        60        70'
;  print,'   12345678901234567890123456789012345678901234567890123456789012345678901234567890'
;  for i=1, num_of_levels DO BEGIN
;;      print, 'level:',i
;      b = WHERE(levels eq i, count)
;      help_s = STRCOMPRESS(STRING(i))
;      for j = 1, count DO help_s = help_s+'#'
;      print, help_s+STRCOMPRESS(string(count))
;  ENDFOR

END
;.compile order_grids
