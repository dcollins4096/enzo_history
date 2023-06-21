FUNCTION grid_in_slice, grid_info, slice_ori, slice_coord
; determine whether the specified grid  is at least partially in the
; slice  ; defined by the 3 component vector slice_ori 
; return for each grid 1 if it is and 0 if it is not in slice
  axes = FIX(slice_ori gt 1.e-30)
  IF (Fix(TOTAL(axes)) gt 1) then BEGIN
      print, 'grid_in_slice: slice_ori:', slice_ori,' not allowed !'
      RETURN, -1
  ENDIF

;  print,'axes',axes
;  print,'slice_ori',slice_ori
;  print,'slice_coord',slice_coord

  const_sub =  where(axes eq max(axes))
  const_sub = const_sub(0)
  other_subs = where(axes eq 0)

;  print, 'const_sub', const_sub
;  print, 'other_subs', other_subs
  
  Left_edge  = grid_info.Left_edge(const_sub,*)
  o_Left_edge  = grid_info.Left_edge(other_subs,*)
  Right_edge = grid_info.Right_edge(const_sub,*)
  o_Right_edge  = grid_info.Right_edge(other_subs,*)

  cond1=(Left_edge  le slice_ori[const_sub])
  cond2=(Right_edge ge slice_ori[const_sub])
  cond3=reform(slice_coord[0] lt o_Right_edge[0,*])
  cond4=reform(slice_coord[1] lt o_Right_edge[1,*])
  cond5=reform(slice_coord[2] gt o_Left_edge[0,*])
  cond6=reform(slice_coord[3] gt o_Left_edge[1,*])

  grid_in_slice=(cond1 and cond2 and cond3 and cond4 and cond5 and cond6)


  print, 'there are '+string(fix(total(grid_in_slice)),format='(i0)')+' grids at least partially in slice'
  return, grid_in_slice
END
