PRO draw_current_cm, draw_widget_id, x_size, min_v, max_v
@TOOLS/common_blocks.inc
image = image_stack(0)

LABEL_COLOR =  !D.table_size 
label_num = 4
; store curent window index
old_wind = !d.window
if x_size gt 0 then widget_control, xsize = x_size, draw_widget_id
widget_control, get_value = win_index, draw_widget_id
wset, win_index
erase
rotation = !p.t
!p.t = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
min_v = FLOAT(min_v) & max_v = FLOAT(max_v)
IF min_v ne max_v THEN BEGIN
;    TVLCT, r_orig, g_orig, b_orig, /GET
;    r0  = r_orig
;    g0  = b_orig
;    b0  = g_orig
    ncolors = !d.table_size
    im_arr = INTARR(!d.x_size,!d.y_size)
    height = 25
;    width = 255
;    x_off = FiX((!d.x_size-width)/2)
    width = 500
    x_off = fix((!d.x_size-width)*0.75)
    y_off =  fix((!d.y_size-height)/1.5)
;    print, x_off, y_off
    im_arr[x_off:(x_off +width-1), y_off:(y_off+height-1)] = $
      bytscl(congrid(indgen(ncolors),width) # replicate(1, height), top = ncolors-1)
    
    loadct,4
    TVSCL, im_arr
;draw bounding box
    box = [[x_off,y_off], [x_off,y_off+height], [x_off+width,y_off+height], $
           [x_off+width,y_off], [x_off,y_off]]
    plots, box, /device, thick=1.5

; put label:
xyouts, 20, 20, image.title, /DEVICE, CHARSIZE=1.7, CHARTHICK=1.5
 print, 'z_current:', z_current
IF (z_current GE 0.) THEN xyouts, 20, 40, 'z='+string(z_current,FORMAT = '(1F7.4)'), /DEVICE, CHARSIZE=1.5, CHARTHICK=1.5
;write labels
    if (n_elements(min_v) eq 0 ) then min_v = 0.
    if (n_elements(max_v) eq 0 ) then max_v = 1.
    min_s = STRING(min_v, FORMAT = '(1F5.2)')
    help  = indgen(label_num)
    labels= STRING(min_v + (max_v-min_v)/FLOAT(label_num)*FLOAT(help+1),  FORMAT = '(1F5.2)')
    xyouts, x_off,y_off-20, min_s, ALIGNMENT = 0.5, /DEVICE,  CHARSIZE=2, CHARTHICK=2
    FOR i=0,label_num-1 DO BEGIN
        xyouts, x_off+(i+1)*width/N_elements(help),y_off-20, labels(i), ALIGNMENT = 0.5, /DEVICE,  CHARSIZE=2, CHARTHICK=2
; put tick marks
        plots, x_off+(i+1)*width/N_elements(help),y_off, /DEVICE
        plots, x_off+(i+1)*width/N_elements(help),y_off+ROUND(height/5), /CONTINUE, /DEVICE
    ENDFOR
    color = ncolors  -1
    cps = [0, ncolors-1]
ENDIF 
wset, old_wind
!p.t = rotation
END
; .compile draw_current_cm
