pro display_iso
COMMON VOLUME_DATA, cube_data, shade_data, shades
common ISO, verteces, polys, is_value, cube_dim, threeD_smoothing,shade_state
   IF  N_ELEMENTS(polys) gt 1 THEN $
     IF ((N_ELEMENTS(shades) gt 1) and (shade_state.this_state eq 1)) THEN $
     TV, POLYSHADE(verteces,polys,SHADES=shades,/T3D) $
   ELSE BEGIN 
       set_shading, LIGHT=shade_state.light_direction
       TV, POLYSHADE(verteces,polys,/T3D) 
   ENDELSE
RETURN
END


