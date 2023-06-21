

; $Id: AMRcube_ORIG.pro,v 1.1 2006/12/14 06:18:25 bwoshea Exp $
;
; Copyright (c) 1994, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.
; (Of course, if you don't work for RSI, remove these lines or
;  modify to suit.)
;+
; NAME:
;	ROUTINE_NAME
;
; PURPOSE:
;	Tell what your routine does here.  I like to start with the words:
;	"This function (or procedure) ..."
;	Try to use the active, present tense.
;
; CATEGORY:
;	Put a category (or categories) here.  For example:
;	Widgets.
;
; CALLING SEQUENCE:
;	Write the calling sequence here. Include only positional parameters
;	(i.e., NO KEYWORDS). For procedures, use the form:
;
;	ROUTINE_NAME, Parameter1, Parameter2, Foobar
;
;	Note that the routine name is ALL CAPS and arguments have Initial
;	Caps.  For functions, use the form:
; 
;	Result = FUNCTION_NAME(Parameter1, Parameter2, Foobar)
;
;	Always use the "Result = " part to begin. This makes it super-obvious
;	to the user that this routine is a function!
;
; INPUTS:
;	Parm1:	Describe the positional input parameters here. Note again
;		that positional parameters are shown with Initial Caps.
;
; OPTIONAL INPUTS:
;	Parm2:	Describe optional inputs here. If you don't have any, just
;		delete this section.
;	
; KEYWORD PARAMETERS:
;	KEY1:	Document keyword parameters like this. Note that the keyword
;		is shown in ALL CAPS!
;
;	KEY2:	Yet another keyword. Try to use the active, present tense
;		when describing your keywords.  For example, if this keyword
;		is just a set or unset flag, say something like:
;		"Set this keyword to use foobar subfloatation. The default
;		 is foobar superfloatation."
;
; OUTPUTS:
;	Describe any outputs here.  For example, "This function returns the
;	foobar superflimpt version of the input array."  This is where you
;	should also document the return value for functions.
;
; OPTIONAL OUTPUTS:
;	Describe optional outputs here.  If the routine doesn't have any, 
;	just delete this section.
;
; COMMON BLOCKS:
;	BLOCK1:	Describe any common blocks here. If there are no COMMON
;		blocks, just delete this entry.
;
; SIDE EFFECTS:
;	Describe "side effects" here.  There aren't any?  Well, just delete
;	this entry.
;
; RESTRICTIONS:
;	Describe any "restrictions" here.  Delete this section if there are
;	no important restrictions.
;
; PROCEDURE:
;	You can describe the foobar superfloatation method being used here.
;	You might not need this section for your routine.
;
; EXAMPLE:
;	Please provide a simple example here. An example from the PICKFILE
;	documentation is shown below.
;
;	Create a PICKFILE widget that lets users select only files with 
;	the extensions 'pro' and 'dat'.  Use the 'Select File to Read' title 
;	and store the name of the selected file in the variable F.  Enter:
;
;		F = PICKFILE(/READ, FILTER = ['pro', 'dat'])
;
; MODIFICATION HISTORY:
; 	Written by:	Your name here, Date.
;	July, 1994	Any additional mods get described here.  Remember to
;			change the stuff above if you add a new keyword or
;			something!
;-



; DO NOT REMOVE THIS COMMENT: END HEADER
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN MAIN13


FUNCTION Rvalue, Eid
     WIDGET_CONTROL,Eid,GET_VALUE=value
     RETURN, value
END

PRO Sfile

;; we need to know the names of the fields if we want read from them
COMMON fields, list_str, var_index, cube_value, xy_sl_size, cube_size, cube_ori, $
   grid_info, data_dir, sds_num, draw_area, draw_cm, $
  MAIN13, base2, base3, cube_text

;; Get values
  WIDGET_CONTROL, CUBE_VAL, GET_VALUE=cube_value
  WIDGET_CONTROL, X_SIZE,    GET_VALUE=xy_sl_size
  WIDGET_CONTROL, left_x,     GET_VALUE=cube_size(0)
  WIDGET_CONTROL, left_y,     GET_VALUE=cube_size(1)
  WIDGET_CONTROL, right_x,     GET_VALUE=cube_size(2)
  WIDGET_CONTROL, right_y,     GET_VALUE=cube_size(3)

END

PRO update_cube_text, widget_id, cube_size
  text_base = strarr(7)
  text_base(0) = 'Window size:'
  text_base(1:4) = string(cube_size)
  text_base(5)   = 'side length:'
  text_base(6)   = string(cube_size(2)-cube_size(0))
WIDGET_CONTROL, widget_id, SET_VALUE = text_base
END

PRO draw_current_cm, draw_widget_id, x_size, min_v, max_v
label_num = 4
; store curent window index
old_wind = !d.window
widget_control, xsize = x_size, draw_widget_id
widget_control, get_value = win_index, draw_widget_id
wset, win_index
erase
TVLCT, r_orig, g_orig, b_orig, /GET
r0  = r_orig
g0  = b_orig
b0  = g_orig
ncolors = !d.table_size
im_arr = INTARR(!d.x_size,!d.y_size)
height = 25
x_off = FiX((!d.x_size-ncolors)/2)
y_off =  FiX((!d.y_size-height)/1.5)
print, x_off, y_off
im_arr(x_off:(x_off +ncolors-1), y_off:(y_off+height-1)) = $
   BYTSCL(INDGEN(ncolors) # REPLICATE(1, height), top = ncolors-1)
TVSCL, im_arr
;draw bounding box
box = [[x_off,y_off], [x_off,y_off+height], [x_off+ncolors,y_off+height], $
       [x_off+ncolors,y_off], [x_off,y_off]]
plots, box, /device
;write labels
if (n_elements(min_v) eq 0 ) then min_v = 0.
if (n_elements(max_v) eq 0 ) then max_v = 1.
min_s = STRING(min_v, FORMAT = '(1F5.2)')
help  = indgen(label_num)
labels= STRING(min_v + (max_v-min_v)/FLOAT(label_num)*FLOAT(help+1),  FORMAT = '(1F5.2)')
xyouts, x_off,y_off-15, min_s, ALIGNMENT = 0.5, /DEVICE, COLOR = [255,255,255]
FOR i=0,label_num-1 DO $
  xyouts, x_off+(i+1)*ncolors/N_elements(help),y_off-15, labels(i), ALIGNMENT = 0.5, /DEVICE, COLOR = [255,255,255]

color = ncolors  -1
cps = [0, ncolors-1]

wset, old_wind
END

PRO dump_ps, draw_id, draw_cm_id, un_sensitive_id, suggest_title


;filename = pickfile(get_path = path, TITLE = 'select name of PS file', FILTER = '*.ps')
filename = './test.ps'

widget_control, un_sensitive_id, sensitive =0 

base = widget_base(title='Dump PostScript file', /COLUMN)   

junk =  WIDGET_LABEL(base, VALUE='Plot Title:',  $
    FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1')
title_text = WIDGET_TEXT(base, /EDITABLE, UVALUE='PLOT_TITLE', $
                         VALUE = suggest_title, YSIZE=2)

i_date = 0
i_eps  = 1
opt_but = [  'Date', 'Encapsulated PS' ]
options = CW_BGROUP(base, /COLUMN,  /NO_REL, /NONEXCLUSIVE, /RETURN_NAME, $
    opt_but, UVALUE='OPTIONS',    LABEL_TOP='Options:', $
    FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1',$
    SET_VALUE = [1,0],  FRAME = 4)


PRINT_BUTTON = WIDGET_BUTTON( BASE, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='PRINT_BUTTON', $
      VALUE='PRINT')

CANCEL_BUTTON = WIDGET_BUTTON( BASE, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='CANCEL_BUTTON', $
      VALUE='CANCEL')


WIDGET_CONTROL, draw_id, GET_VALUE=win
WSET, win   
WINDOW, /PIXMAP, /FREE, xs = !d.x_size, ys=!d.y_size  ;Save window
wset, !d.window
DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, win]  ;Save it
backing = tvrd()
backing_num = !d.window
WIDGET_CONTROL, draw_cm_id, GET_VALUE=win2
WSET, win2
WINDOW, /PIXMAP, /FREE, xs = !d.x_size, ys=!d.y_size  ;Save CM window
wset, !d.window
DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, win2]  ;Save it
cm_backing = tvrd()
cm_backing_num = !d.window
wset, win
 
back_s    =  size(backing)   
cm_back_s =  size(cm_backing)   

comb_im   = bytarr(back_s(1), back_s(2)+cm_back_s(2)+80)
comb_im(0:(back_s(1)-1), cm_back_s(2):(cm_back_s(2)+back_s(2)-1)) = backing
comb_im(0:(back_s(1)-1), 0:(cm_back_s(2)-1)) = cm_backing
comb_im_s = size(comb_im)
print, size(comb_im)

WIDGET_CONTROL, base, /REALIZE

WHILE 1 DO BEGIN                ;Internal event loop   
    ev = WIDGET_EVENT([base])
    WIDGET_CONTROL,Ev.id,GET_UVALUE=Evu

    CASE Evu OF   
        'CANCEL_BUTTON':  BEGIN   
            print, 'event for cancel button'
            goto, all_done
        ENDCASE    
        
        'PRINT_BUTTON': BEGIN   
; read options
            WIDGET_CONTROL, options, get_value = buts
            WIDGET_CONTROL, title_text, get_value = title
            print, buts
            print, title
;             
            tvlct, r,g,b, /get
;            window, /free
;            tv, comb_im
            set_plot, 'PS'
            DEVICE, file = filename, /COLOR, BITS=8
            IF (buts(i_eps) gt 0) THEN BEGIN
                print, 'writing encapsulated PS file:', filename
                DEVICE, /ENCAPSULATED
            ENDIF $
            ELSE print, 'write standard PS file:', filename
            
            
; switch black and white      
;            white = where(comb_im EQ 255B)
;            black = where(comb_im EQ 0B)
;            print, size(white)
;            comb_im(white(*)) = 0B
;            comb_im(black(*)) = 255B
            tvlct, r,g,b
            TV, comb_im

            XYOUTS, 50,FIX(!d.y_size*0.95), title(0), /DEVICE, CHARSIZE = 1.3, COLOR=100
            XYOUTS, 50,FIX(!d.y_size*0.92), title(1), /DEVICE, CHARSIZE = 0.7, COLOR=100
            if (buts(i_date) gt 0) THEN BEGIN
                print, 'd.x_size:', !d.y_size, systime()
                XYOUTS, 0,0, systime(), $
                  /DEVICE, CHARSIZE = 0.8, COLOR=100, ALIGNMENT=0.
            END
            DEVICE, /CLOSE
            set_plot, 'x'

all_done:
            WIDGET_CONTROL, base, /DESTROY   
            widget_control, un_sensitive_id, sensitive = 1 
            return
        ENDCASE    
        ELSE: print, 'unspecified event', Evu
    ENDCASE   
    
ENDWHILE                        ;Event loop


END



PRO MAIN13_Event, Event
COMMON fields, list_str, var_index, cube_value, xy_sl_size,  cube_size, cube_ori, $
  grid_info, data_dir, sds_num, draw_area, draw_cm, $
  MAIN13, base2, base3, cube_text

WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev

CASE Ev OF 
    
    'BGROUP10': BEGIN
        CASE Event.Value OF
            0: begin 
                Print,'Button x  Pressed'
                i_s = 0
                cube_ori = [cube_value, 0, 0]
            end
            1: begin
                Print,'Button y  Pressed'
                i_s = 1
                cube_ori = [0, cube_value,  0]
            end
            2: begin
                Print,'Button z  Pressed'
                i_s = 2
                cube_ori = [0, 0, cube_value]
            end
            ELSE: Message,'Unknown button pressed'
        ENDCASE
    END
    'CUBE_VAL': BEGIN
        cube_value = Rvalue(Event.id)
        Print, 'Event for cube_val ', cube_value
        ind = where(cube_ori gt 1.e-5)
        cube_ori(ind) = cube_value
    END
    'X_SIZE': BEGIN
        xy_sl_size = RVALUE(Event.id)
        xy_sl_size = max([250, xy_sl_size])
        Print, 'Event for x_size:', xy_sl_size
        widget_control, draw_area, $
          draw_xsize=xy_sl_size, draw_ysize=xy_sl_size
        draw_current_cm, draw_cm, xy_sl_size, min_data, max_data

    END

    'VAR_LIST': BEGIN
        Print, 'Event for variable list'
        Print, 'index:', Event.index
        var_index = Event.index
        list_variables = [1,16,7,8,9,10,11,12,13,14,15,2,3,4,5,6]
        sds_num = list_variables(Event.index)
    END
    'HELP_BUTTON': BEGIN
        Print, 'Event for HELP'
    END
    'DRAW_BUTTON': BEGIN
        Print, 'Event for DRAW'
        widget_control, /hourglass
;      print, grid_info
        construct_cube_data, grid_info, cube_ori, cube_size, data_dir, $
          sds_num, [xy_sl_size,xy_sl_size], data
        sec_min_data = min(data(where(data ne min(data))))
        data = alog10(data > sec_min_data)
        TVSCL, data
        min_data = min(data)
        max_data = max(data)
        print, '2nd min data:', min_data, max_data
        draw_current_cm, draw_cm, xy_sl_size, min_data, max_data
        update_cube_text, cube_text, cube_size
        
    END
    'SHOW_GRIDS': BEGIN
       draw_grids, grid_info, cube_ori, cube_size, xy_sl_size
    END
    'CM_BUTTON': BEGIN
        Print, 'Event for Color Map'
        xloadct
    END
    'ZOOM_IN_BUTTON': BEGIN
        Print, 'Event for ZOOM'
        widget_control, base2, sensitive = 0
        widget_control, base3, sensitive = 0
        Q = CW_DEFSQUARE(draw_area)
        if (N_ELEMENTS(Q) ge 9) THEN BEGIN
            wind = FLOAT(Q(*,0:3))/FLOAT(xy_sl_size)
            left_x  = min(wind(0,*))
            left_y  = min(wind(1,*))
            right_x = max(wind(0,*))
            right_y = max(wind(1,*))
            print, 'wind:',wind
            windx = cube_size(2)-cube_size(0)
            windy = cube_size(3)-cube_size(1)
            print, 'windy:', windy
            cube_size(2) = cube_size(0) + windx*right_x
            cube_size(3) = cube_size(1) + windy*right_y
            cube_size(0) = cube_size(0) + windx*left_x
            cube_size(1) = cube_size(1) + windy*left_y
        ENDIF
        print, 'cube_size', cube_size
        widget_control, base2, sensitive = 1
        widget_control, base3, sensitive = 1
        update_cube_text, cube_text, cube_size

    END
    'ZOOM_OUT_BUTTON': BEGIN
        cube_size = [0.,0.,1.,1.]
        widget_control, /hourglass
;      print, grid_info
        construct_cube_data, grid_info, cube_ori, cube_size, data_dir, $
          sds_num, [xy_sl_size,xy_sl_size], cube_data
        TV, bytscl(alog10(cube_data>.1))
        update_cube_text, cube_text, cube_size
    END

    'INFO_BUTTON': BEGIN
        Print, 'Event for INFO'
        Print, 'sds_num:',sds_num
        Print, 'xy_sl_size:', xy_sl_size
        Print, 'size(grid_info)', size(grid_info)
        Print, 'cube_ori:', cube_ori
        Print, 'cube_size:', cube_size
        Print, 'data_dir:', data_dir
    END
    
    'PS_BUTTON': BEGIN
        PRINT, 'Event for dump PS file'
        dump_ps, draw_area, draw_cm, main13, $
          [list_str(var_index),'sidelength ='+ STRCOMPRESS(string(cube_size(2)-cube_size(0)))]
    END

    'DONE_BUTTON': BEGIN
        Print, 'Event for DONE'
        widget_control, /DESTROY, Event.top
    END
    ELSE: 
ENDCASE
END


; DO NOT REMOVE THIS COMMENT: END MAIN13
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.



PRO AMRcube, grid_in_info, data_in_dir ; GROUP=Group
COMMON fields, list_str, var_index,  cube_value, xy_sl_size,  cube_size, cube_ori, $
   grid_info, data_dir, sds_num, draw_area, draw_cm, $
  MAIN13, base2, base3, cube_text

; initial values
  cube_value = 0.282104
  cube_ori   = [cube_value,0,0]
  cube_size  = [0.,0.,1.,1.]
  xy_sl_size = 450

  grid_info = grid_in_info
  data_dir  = data_in_dir
  sds_num   = 1

;  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0
  group = 0
  junk   = { CW_PDMENU_S, flags:0, name:'' }

  MAIN13 = WIDGET_BASE(GROUP_LEADER=Group, $
      ROW=1, $
      MAP=1, FRAME = 4, $
      TITLE='LCA AMR Cuber', $
      UVALUE='MAIN13')

  BASE2 = WIDGET_BASE(MAIN13, $
      COLUMN=1, $
      MAP=1, FRAME = 4, $
      UVALUE='BASE2')

  BASE3 = WIDGET_BASE(MAIN13, $
      COLUMN=1, $
      MAP=1, FRAME = 4, $
      UVALUE='BASE2')

  BASE4 = WIDGET_BASE(MAIN13, $
      COLUMN=1, $
      MAP=1, FRAME = 4, $
      UVALUE='BASE2')

  Btns123 = [ $
    'x ', $
    'y ', $
    'z ' ]
  i_ind=  FIX(where(cube_ori gt 0))
  BGROUP10 = CW_BGROUP( BASE2, Btns123, $
      SET_VALUE = i_ind(0), $
      ROW=1, $
      /EXCLUSIVE, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      LABEL_TOP='Cube Specifics', $
      UVALUE='BGROUP10')

  CUBE_VAL = CW_FSLIDER( BASE2, $
      MAXIMUM=1.00000, $
      MINIMUM=0.00000, $
      TITLE='cube coord', $
      UVALUE='CUBE_VAL', $
      VALUE=cube_value, /EDIT)

  X_SIZE = CW_FIELD( BASE2,VALUE=xy_sl_size, $
      COLUMN=1, $
      INTEGER=1, $
      RETURN_EVENTS=1, $
      TITLE='x size (pixels)', $
      UVALUE='X_SIZE')

  CUBE_TEXT = WIDGET_TEXT(BASE2, $
     YSIZE = 7, UVALUE='CUBE_TEXT')
     

  list_str = ['Baryon Density', 'Temperature', 'e-    Density', 'HI    Density', 'HII   Density', 'HeI   Density', 'HeII  Density', 'HeIII Density', 'H-    Density', 'H2I   Density', 'H2II  Density', 'Total Energy', 'Gas   Energy', 'x-velocity', 'y-velocity', 'z-velocity']
  var_index= 0
  VAR_LIST = WIDGET_LIST( BASE2,VALUE=List_str, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='VAR_LIST', $
      YSIZE=10)
  widget_control, var_list, set_list_select = var_index

  HELP_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='HELP_BUTTON', $
      VALUE='HELP')

  DRAW_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='DRAW_BUTTON', $
      VALUE='DRAW')

  SHOW_GRIDS  =  WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='SHOW_GRIDS', $
      VALUE='SHOW GRIDS')

  CM_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='CM_BUTTON', $
      VALUE='Color Map')

  ZOOM_IN_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='ZOOM_IN_BUTTON', $
      VALUE='ZOOM_IN')

  ZOOM_OUT_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='ZOOM_OUT_BUTTON', $
      VALUE='ZOOM_OUT')

  INFO_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='INFO_BUTTON', $
      VALUE='INFO')

  PS_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='PS_BUTTON', $
      VALUE='DUMP PostScript')

  DONE_BUTTON = WIDGET_BUTTON( BASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='DONE_BUTTON', $
      VALUE='DONE')

  DRAW_AREA = WIDGET_DRAW(BASE4, $
      XSIZE = xy_sl_size, $
      YSIZE = xy_sl_size, $
      /BUTTON, /MOTION, $
      UVALUE = 'draw_area')

  DRAW_CM = WIDGET_DRAW(BASE4, $
      XSIZE = xy_sl_size, $
      YSIZE = 80, $
      /BUTTON, /MOTION, $
      UVALUE = 'draw_cm')

  WIDGET_CONTROL, MAIN13, /REALIZE
; now get the window number for the draw area
  draw_current_cm, draw_cm, xy_sl_size, 0., 1.
  update_cube_text, cube_text, cube_size
  
  widget_control, get_value = win_index, draw_area
; make it the current graphics window
  wset, win_index
  device, RETAIN=2

  XMANAGER, 'MAIN13', MAIN13


END
