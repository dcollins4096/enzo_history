pro update_shading
COMMON VOLUME_DATA, cube_data, shade_data, shades
COMMON ISO, verteces, polys, is_value, cube_dim, threeD_smoothing,shade_state
COMMON fields, list_str,mass_per_particle, divide_by_density, $
  var_index, center, slice_value, xy_sl_size, slice_size, slice_ori, $
   grid_info, data_dir, sds_num
COMMON widgets,  draw_area, draw_cm, TOP_BASE, base2, base3, slice_val_w, $
  slice_text, ISO_BASE, wISO, wSLICER , IS_VAL, peak_ids, peak_opt, snap_base

  construct_cube, shade_data, shade_state.sds, MIN=min_shade,$
    MAX=max_shade
  shade_state.max = max_shade
  shade_state.min = min_shade

  help,min_shade,max_shade

; set sliders for min and max
  widget_control, shade_state.min_s,$
    set_value=[DOUBLE(FLOOR(min_shade)), min_shade, $
               0.99*shade_state.max]
  widget_control, shade_state.max_s,$
    set_value=[1.01*min_shade,shade_state.max , $
               DOUBLE(CEIL(shade_state.max))]

  compute_em, draw_cm, xy_sl_size, shade_em=shade_state.this_state 
  display_iso
RETURN
END

