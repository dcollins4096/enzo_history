pro read_box_from_window, vec
; use device cursor to define a box
; vec is a four component vector that specifies left and right edge
; of a square !
  vec = dblarr(4)
  cursor, x1, y1, /norm, /down
  if (!err ne 4) then cursor, x2, y2, /norm, /down
  if (!err ne 4) then begin
      print, 'box', x1,y1,x2,y2
      vec(0) = min([x1,x2])
      vec(1) = min([y1,y2])
      vec(2) = max([x1,x2])
      vec(3) = max([y1,y2])
      side_length = max([vec(1)-vec(0), vec(3)-vec(2)])
      vec(2) = vec(0) + side_length
      vec(3) = vec(1) + side_length
  endif
END
