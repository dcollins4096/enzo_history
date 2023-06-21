pro show_image, NO_CM=nocm
; show current image from stack,
; update colormap, and display image information
@TOOLS/common_blocks.inc

widget_control, get_value = win_index, draw_area
; make it the current graphics window
wset, win_index

image = image_stack(0)
; set correct size of drawing area
if (image.xy_size ne xy_sl_size) then BEGIN
    widget_control, draw_area, $
      draw_xsize=image.xy_size, draw_ysize=image.xy_size
    xy_sl_size = image.xy_size
ENDIF

widget_control, get_value = win_index, draw_area
wset, win_index
min_data =  image.min_data
max_data =  image.max_data

;min_data=-2.0
;max_data=2.0

; check on the min and max values:
IF N_ELEMENTS(mi_ma_base) GT 0 THEN $
 IF widget_info(mi_ma_base_id,/VALID_ID) THEN BEGIN
   alive =  widget_info(mi_ma_base_id, /REALIZED)
   IF alive THEN BEGIN 
      widget_control, mi_ma_base_id, get_uvalue=state
      IF state.fixed THEN BEGIN
         min_data = state.mi_ma(0)
         max_data = state.mi_ma(1)
      ENDIF
   ENDIF
ENDIF

image_temp = BYTSCL(image.image_d(0:xy_sl_size-1,0:xy_sl_size-1), $
                    MIN=min_data, MAX=max_data, TOP=!d.table_size-2)

TV, image_temp
IF NOT KEYWORD_SET(NOCM) THEN BEGIN
   draw_current_cm, draw_cm, xy_sl_size, min_data, max_data
   update_slice_text, slice_text, image.slice_size, image.title, [0.,0.,0.,0.]
ENDIF

END
;.compile show_image
