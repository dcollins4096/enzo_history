PRO update_slice_text, widget_id, slice_size, title, pofi, SHOW_POFI=SHOW_POFI, $
                       SHOW_GRID=grid_info
  text_base = strarr(8)
  text_base(0) = title
  text_base(1) = ' slice size:'
  text_base(2:5) = ' '+STRCOMPRESS(string(slice_size))
  text_base(6)   = 'side length:'
  text_base(7)   = ' '+STRCOMPRESS(string(slice_size(2)-slice_size(0)))

  IF KEYWORD_SET(SHOW_POFI) THEN BEGIN
      IF ARG_PRESENT(grid_info) THEN $
        finest_grid_containing_point, grid_info, pofi(0:2), grid_num, grid_level
      
      text_base(1) = text_base(1) + '    point:'
      text_base(2) = text_base(2) + ''+ string(pofi(0))
      text_base(3) = text_base(3) + ''+ string(pofi(1))
      text_base(4) = text_base(4) + ''+ string(pofi(2))
      text_base(5) = text_base(5) + ''+ string(pofi(3))
      text_base(6) = text_base(6) + ''+ '    G'+ $
        STRCOMPRESS(string(grid_num),/REMOVE_ALL) + 'L' + $
        STRCOMPRESS(string(grid_level),/REMOVE_ALL)
  ENDIF

WIDGET_CONTROL, widget_id, SET_VALUE = text_base

RETURN
END
