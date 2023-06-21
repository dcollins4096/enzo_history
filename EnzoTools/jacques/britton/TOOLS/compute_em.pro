pro compute_em, draw_cm, xy_sl_size, SHADE_EM=SHADE_EM
; compute polygons for isosurface and determin shading values
COMMON VOLUME_DATA, cube_data, shade_data, shades
COMMON ISO, verteces, polys, is_value, cube_dim, threeD_smoothing,shade_state
COMMON options, verbose

s = size(cube_data)
if s(0) ne 3 THEN BEGIN
    Print, 'compute_em: Array is not 3 D !'
    RETURN
ENDIF
scale3, xrange = [0,s(1)], yrange = [0,s(2)], zrange = [0,s(3)]
shade_volume, cube_data, is_value, verteces, polys,/LOW, /VERBOSE

IF KEYWORD_SET(SHADE_EM) THEN BEGIN
    vert_size=size(verteces)
;    print,'vert_size:',vert_size
    shades = DBLARR(vert_size(4))

    i=0l
    jend = LONG(ROUND(N_ELEMENTS(shades)/3)-1)
    for j=0l, jend DO BEGIN
        shades(i:(i+2)) = shade_data(ROUND(verteces(0,j)),$
                                   ROUND(verteces(1,j)),$
                                   ROUND(verteces(2,j)))
        i=i+3 
    ENDFOR    
    min_shade = min(shade_data,max=max_shade)

    IF verbose THEN print, $
      'min and max values of temp on surface:', min(shades), max(shades)
    shades = BYTSCL(shades, TOP=!D.table_size,$
                    MIN=shade_state.min,MAX=shade_state.max)
    smoothr= MIN([60,FIX(N_elements(shades)/33.)])
    IF smoothr GT 1 THEN shades = SMOOTH(shades, smoothr, /EDGE_TRUNCATE)
    shade_max = shade_state.max & shade_min = shade_state.min
    draw_current_cm, draw_cm, xy_sl_size, shade_min, shade_max
END ELSE BEGIN
    shades = 0
    set_shading , LIGHT=shade_state.light_direction
    draw_current_cm, draw_cm, 0, 0, 0 ; erase cm window
ENDELSE

RETURN
END
