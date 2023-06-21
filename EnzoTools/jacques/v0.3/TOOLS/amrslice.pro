;+
; NAME: AMRSlicer is the heart of Jaques who is enzo`s best friend.
; PURPOSE: Visualization of AMR data
; MODIFICATION HISTORY:
; Tom Abel :     June 1998 - 
;-
;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/17/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + sds_num now zero based
;; + DUMP_4D2 disabled
;; + INSIGHT disabled
;; + added RELOAD and BASE_NAME keyword
;; + added a pointer array to all the loaded grids. saves time when a
;;   new slice is requested, since only new grids have to be read from
;;   disk.
;;
FUNCTION Rvalue, Eid
     WIDGET_CONTROL,Eid,GET_VALUE=value
     RETURN, value
END

PRO TOP_BASE_Event, Event
@TOOLS/common_blocks.inc
COMMON VOLUME_DATA, cube_data, shade_data, shades 
; set the appropriate graphics window
  widget_control, get_value = win_index, draw_area
  wset, win_index

WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev

;IF ev.id eq draw THEN BEGIN
;    ev = WIDGET_EVENT([base, draw])
;    IF N_elements(ev.press) gt 0 THEN BEGIN 
;        if ev.press eq 2 then goto, zoom_in
;    END ELSE BEGIN

        CASE Ev OF 

            'draw_area':BEGIN
                if (substracting_vel lt 1) THEN BEGIN
                    if event.press ne 0 then begin
                        x = event.x/DOUBLE(xy_sl_size)*(slice_size(2)-slice_size(0))
                        y = event.y/DOUBLE(xy_sl_size)*(slice_size(3)-slice_size(1))
                        pofi = DBLARR(4)
                        pofi(0:2) = center

;                        pofi[3]   = image_stack[0].image_d[event.x-1,event.y-1]/$
;                          DOUBLE(!d.table_size)*$
;                          (image_stack(0).max_data-image_stack(0).min_data)+$
;                          image_stack(0).min_data

                        pofi[3]=image_stack[0].image_d[event.x-1,event.y-1]

                        os = where(slice_ori lt 1.e-15)
                        pofi(os) = [x,y]+slice_size(0:1)
                        update_slice_text, slice_text, image_stack(0).slice_size,$
                          image_stack(0).title, pofi, /SHOW_POFI, SHOW_GRID=grid_info
                    END 
                END ELSE BEGIN
                    IF N_elements(event.press) gt 0 THEN BEGIN 
                        if (event.press eq 1) then begin
                            x = event.x
                            y = event.y
                            if verbose then print, 'x, y, event.press:',x, y, event.press
                            cons_vel = [vel_one(x,y), vel_two(x,y)]
                            if verbose then print, 'cons_vel :',cons_vel 
                        endif                            
; cancel  ? 
                        if (event.press eq 4) then substracting_vel = 0
                    ENDIF
                ENDELSE
           END
           'IS value': BEGIN
               WIDGET_CONTROL,Event.id,GET_VALUE=tmp
               is_value=tmp(1)
           END
           'BTNS': BEGIN
               EVV = Event.Value
               CASE EVV OF
                   'Compute': BEGIN
                       create_view, Xmax=(cube_dim-1),$
                         YMAX=(cube_dim-1),ZMAX=cube_dim-1,$
                         AX=0.,AY=0.,AZ=0., $
                         WINX=xy_sl_size-1,WINY=xy_sl_size-1, ZOOM=0.57

                       Print, 'Event for compute polys'
                       compute_em, draw_cm, xy_sl_size,$
                         SHADE_EM=shade_state.this_state
                       display_iso
                   END
                   "Weigh": BEGIN
                       Print, 'Event for WEIGH'
                       Print, 'Note: that this makes only sense'
                       Print, 'if the volume data is a over density !'
                       zones_in = N_ELEMENTS(WHERE(Cube_data ge IS_VALUE))
                       side_length = DOUBLE(slice_size(2)-slice_size(0))
                       s = size(Cube_data)
                       zones_tot= s(5)
	               cube_dim = s(1)
                       del_x    = side_length/cube_dim
                       dens = 10.^Cube_data
                       mass_tot = TOTAL(dens) * del_x^3
                       mass_in  = TOTAL(dens(WHERE(Cube_data GT IS_VALUE)))* $
                              del_x^3.                       
                       mass_out = mass_tot - mass_in
                       PRINT, 'zones inside:', zones_in
                       PRINT, 'mass outside:',mass_out
                       PRINT, 'zones, mass total:',  zones_tot, mass_tot
                       PRINT, $
		          'mass inside (mean density * simulation volume (com Mpc^3):', $
                         mass_in
                       Print, 'this is in SolarMass per (boxlength/comoving Kpc)^3'
		       Print, 'so multiply by 2.76 10^11 M_sun/Mpc^3 * '+ $
                       'Omega_matter * h^2 (sim sidelength)^3'
		       Print, 'for first structure formation (128kpc):', mass_in*.128^3.*2.76e11/4.
                   END
                   '(Re)Read': BEGIN
                       widget_control, /hourglass
                       construct_cube, cube_data, sds_num,  MIN=min_v, MAX=max_v
                       max_v =  0.99*max_v
                       is_value = 0.4*max_v
                       is_slider = [min_v, is_value, max_v]
                       widget_control, IS_val, set_value = is_slider
                      
                   END
                   'IDL SLICER': BEGIN
;                       h_cube_data = PTR_NEW(cube_data, /NO_COPY)
;                       slicer3, h_cube_data, GROUP=TOP_BASE
                      slicer
                   END
                   ELSE: print, 'unrecognized button event'
               END
           ENDCASE
           'ORI': BEGIN
               display_iso
           END
           'VELBTNS': BEGIN
               EVV = Event.Value
               CASE EVV OF
                   "Read": BEGIN
; make sure the original velocity field will be plotted first
                       cons_vel = [0., 0.]
                       widget_control, /hourglass
                       read_velocities
                   END
                   "Hide": BEGIN
                       show_image
                   END
                   "Plot": BEGIN
                      IF N_ELEMENTS(vel_one) gt 1 THEN BEGIN
                         widget_control, /hourglass
                         show_image
                         XP=FIX((INDGEN(ROUND(xy_sl_size/vel.skip))+1) $
                                *vel.skip-1.)
                         U = DBLARR(N_ELEMENTS(XP),N_ELEMENTS(XP))
                         V = U
                         FOR i=0,N_ELEMENTS(XP)-1 DO BEGIN
                            U(*,i) =vel_one(XP,XP(i))-cons_vel(0)
                            V(*,i) =vel_two(XP,XP(i))-cons_vel(1)
                         ENDFOR
                         !p.noclip=0
                         print, max(vel_one-cons_vel(0)),max(vel_two-cons_vel(1))
                         !p.position=[0,0,xy_sl_size,xy_sl_size]
                         !p.ticklen = 0.
                         MYVELOVECT, U,V, XP, XP,$
                          LENGTH=vel.amplify, /NOERASE, /DEVICE, COLOR=vel.color
                      END ELSE BEGIN
                         print, 'no velocity data!'
                         print, 'must read in velocities first!'
                      END
                   END 
                   "Substract": Begin
                       if (N_ELEMENTS(vel_one) gt 1) then begin
                           if verbose then begin 
                               dum = STRARR(7)
                               dum(0) =  "Please klick at a point in "
                               dum(1) =  "the drawing window (left button)"
                               dum(2) =  "The velocity at that point will be "
                               dum(3) =  "substracted the next time you plot  "
                               dum(4) =  "the velocity vectors.  "
                               dum(5) =  "Cancel: right mouse button"
                               dum(6) =  "Press Reset to get original vel. field"
                               xdisplayfile, ' ', TEXT=dum, GROUP=TOP_BASE, $
                                 TITLE='INFORMATION', HEIGHT=7, WIDTH = 40
                           endif
                           substracting_vel = 1
                       END ELSE BEGIN
                           print, 'no velocity data!'
                           print, 'must read in velocities first!'
                       ENDELSE
                               
                       END
                   "Reset": Begin
                       cons_vel = [0., 0.]
                   END
                   ELSE: print, 'unrecognized button event'
               END 
           END
           'Vel_skip':BEGIN
               WIDGET_CONTROL,Event.id,GET_VALUE=skip_h
               vel.skip=skip_h
           END

           'Vel_amp':BEGIN
               WIDGET_CONTROL,Event.id,GET_VALUE=amplify_h
               vel.amplify =amplify_h(1)
           END
           "VEL_COL_SEL":BEGIN
               vel.color = event.value
           END
           
    'BGROUP10': BEGIN
        CASE Event.Value OF
            0: begin 
                If Verbose Then Print,'Button x  Pressed'
                i_s = 0
                slice_value = center(i_s)
                slice_ori = [slice_value, 0, 0]
                
            end
            1: begin
                If Verbose Then Print,'Button y  Pressed'
                i_s = 1
                slice_value = center(i_s)
                slice_ori = [0, slice_value,  0]
            end
            2: begin
                If Verbose Then Print,'Button z  Pressed'
                i_s = 2
                slice_value = center(i_s)
                slice_ori = [0, 0, slice_value]
            end
            ELSE: Message,'Unknown button pressed'
        ENDCASE
        windx = 0.5*(slice_size(2)-slice_size(0))
        os = where([0,1,2] ne Event.value)
        REPEAT BEGIN
            slice_size(0)   = center(os(0))-windx
            slice_size(1)   = center(os(1))-windx
            slice_size(2)   = center(os(0))+windx
            slice_size(3)   = center(os(1))+windx
            windx = 0.9*windx
        END UNTIL (min(slice_size ge 0.) and max(slice_size le 1.0))
        slider = [(slice_value-windx)>0., slice_value, (slice_value+windx)<1.]
        widget_control, slice_val_w, set_value = slider
    END
    'SLICE_VAL': BEGIN
        slice_value = (Rvalue(Event.id))(1)
        If Verbose Then Print, 'Event for slice_val ', slice_value
        ind = where(slice_ori gt 1.e-15)
        slice_ori(ind) = slice_value
    END
    'X_SIZE': BEGIN
        xy_sl_size = RVALUE(Event.id)
        xy_sl_size = max([250, xy_sl_size])
; make sure nobody tries to plot velocities again
        vel_one = 0.
        vel_two = 0.
        If Verbose Then Print, 'Event for x_size:', xy_sl_size
        widget_control, draw_area, $
          draw_xsize=xy_sl_size, draw_ysize=xy_sl_size
        widget_control, base2, scr_ysize=xy_sl_size+100
        draw_current_cm, draw_cm, xy_sl_size, min_data, max_data

    END

    'SMOOTH_VALUE': BEGIN
        WIDGET_CONTROL,Event.id, get_value=tmp
        smoothing = tmp(1)
    END

    'BGROUP_INTER': BEGIN
        CASE Event.Value OF
            0: begin 
                If Verbose Then Print,'Button none  Pressed'
                interpolate_i = 0
            end
            1: begin
                If Verbose Then Print,'Button linear  Pressed'
                interpolate_i = 1
            end
            2: begin
                If Verbose Then Print,'Button quintic  Pressed'
                interpolate_i = 2
            end
            ELSE: Message,'Unknown button pressed'
        ENDCASE
        IF image_stack(0).title NE '' THEN BEGIN 
           widget_control, /hourglass
           render_image, result, SECMIN=image_stack(0).min_data
           image_stack(0).image_d(0:xy_sl_size-1,0:xy_sl_size-1) = result
           show_image, /NO_CM
        ENDIF 
    END

    'AN_CM': BEGIN
        ANNOTATE, DRAWABLE=draw_cm
    END

    'AN_DRAW': BEGIN
        ANNOTATE, DRAWABLE=draw_area
    END

    'EDITOR': BEGIN
        myeditor
    END
    
    'Calculator': BEGIN
        mycalc
    END
    
    'XPALETTE': BEGIN
        xpalette, GROUP=TOP_BASE
    END

    'VAR_LIST': BEGIN
        If Verbose Then Print, 'Event for variable list'
        If Verbose Then Print, 'index:', Event.index
        var_index = Event.index
;        list_variables = [1,16,7,8,9,10,11,12,13,14,15,2,3,4,5,6]
;        sds_num = list_variables(Event.index)
;        sds_num = Event.index+1
        sds_num=event.index
   END
    'HELP': BEGIN
        If Verbose Then Print, 'Event for HELP'
        print, 'Sorry no HELP available yet ...'
        print, 'check out, the README, WHATS_NEW, and TODO files'
        print, 'in the distribution for more information.'
    END
    'About': BEGIN
        dum = STRARR(5)
        dum(0) =  "Jaques is enzo's best friend"
        dum(1) =  "Author: Tom Abel        1998"
        dum(2) =  "Jaques is distributed in the hope it is helpful...."
        dum(3) =  "For scientific app's I ask for an Acknowledgement"
        dum(4) =  "send bug reports to tabel@cfa.harvard.edu"
        dummy  = dialog_message(dum, /INFORMATION)
    END
    'DRAW_BUTTON': BEGIN
        if verbose then  Print, 'Event for DRAW'
        widget_control, /hourglass
        construct_image
        show_image
    END
    'SHOW_CHART': BEGIN
        draw_grids, 2, $
          draw_cm, grid_info, slice_ori, slice_size, xy_sl_size
    END
    'SHOW_BOXES': BEGIN
        draw_grids, 1, $
          draw_cm, grid_info, slice_ori, slice_size, xy_sl_size
    END
    'SHOW_PEAKS': BEGIN
        show_peaks
    END
    'CM_BUTTON': BEGIN
        If Verbose Then Print, 'Event for Color Map'
        xloadct, GROUP=TOP_BASE
    END
  'VIEW_HIERARCHY': BEGIN
        If Verbose Then Print, 'Event for VIEW_HIERARCHY'
        if N_ELEMENTS(grid_info(0)) gt 0 THEN $
            myeditor, grid_info(0).hier_file $
          ELSE message, 'No grid info'
        END

  'VIEW_PARAMETER': BEGIN
        If Verbose Then Print, 'Event for VIEW_PARAMETER'

        if N_ELEMENTS(grid_info(0)) gt 0 THEN BEGIN
            par_file = strmid(grid_info(0).hier_file, 0, $
                              STRPOS(grid_info(0).hier_file,'.hierarchy'))
            myeditor, par_file
        END ELSE message, 'No grid info'

        END

    'ZOOM_IN_BUTTON': BEGIN
        If Verbose Then Print, 'Event for ZOOM'
        windx = DOUBLE(slice_size(2)-slice_size(0))
        if (EVENT.value eq 'ZOOM IN.pick') THEN BEGIN
 zoom_in:        
            widget_control, base2, sensitive = 0
            widget_control, base3, sensitive = 0
            Q = CW_DEFSQUARE(draw_area, xy_sl_size, slice_size(2)-slice_size(0))
            if (N_ELEMENTS(Q) ge 9) THEN BEGIN
                wind = DOUBLE(Q(*,0:3))/DOUBLE(xy_sl_size)
                left_x  = min(wind(0,*))
                left_y  = min(wind(1,*))
                right_x = max(wind(0,*))
                right_y = max(wind(1,*))
                if verbose then print, 'wind:',wind
                windx = slice_size(2)-slice_size(0)
                windy = slice_size(3)-slice_size(1)
                if verbose then print, 'windy:', windy
                slice_size(2) = DOUBLE(slice_size(0) + windx*right_x)
                slice_size(3) = DOUBLE(slice_size(1) + windy*right_y)
                slice_size(0) = DOUBLE(slice_size(0) + windx*left_x)
                slice_size(1) = DOUBLE(slice_size(1) + windy*left_y)
; redefine center:
                os = where(slice_ori lt 1.e-15)
                center(os) = DOUBLE(slice_size(2:3)+slice_size(0:1))/2.D
            ENDIF
        END ELSE BEGIN 
            CASE Event.value OF 
                'ZOOM IN.by 2': factor = 1.D/2.D
                'ZOOM IN.by 4': factor = 1.D/4.D
                'ZOOM IN.by 10': factor = 1./10.D
                'ZOOM IN.by 100': factor = 1./100.D
                'ZOOM IN.by 1000': factor = 1./1000.D
                ELSE:
            ENDCASE
            os = where(slice_ori lt 1.e-15)
            windx = factor*0.5D*DOUBLE(slice_size(2)-slice_size(0))
            REPEAT BEGIN
                slice_size(0)   = center(os(0))-windx
                slice_size(1)   = center(os(1))-windx
                slice_size(2)   = center(os(0))+windx
                slice_size(3)   = center(os(1))+windx
                windx = 0.9*windx
            END UNTIL (min(slice_size ge 0.) and max(slice_size le 1.0))
        ENDELSE
        if verbose then print, 'slice_size', slice_size, 'side_length:', slice_size(2)-slice_size(0)
        slider = [(slice_value-windx)>0., slice_value, (slice_value+windx)<1.]
        widget_control, slice_val_w, set_value = slider

        widget_control, base2, sensitive = 1
        widget_control, base3, sensitive = 1

    END
    'ZOOM_OUT_BUTTON': BEGIN
        widget_control, /hourglass
        side_length = DOUBLE(slice_size(2)-slice_size(0))
        mid_point   = DOUBLE(slice_size(2:3)+slice_size(0:1))/2.D
        old_sl = DOUBLE(slice_size)
        factor =  1.D
        IF (NOT (EVENT.value EQ 'ZOOM OUT.top')) THEN BEGIN
           CASE EVENT.value OF 
              'ZOOM OUT.by 2': BEGIN
                 factor =  2.D
              END
              'ZOOM OUT.by 5': BEGIN
                 factor =  5.D
              END
              'ZOOM OUT.by 10': BEGIN
                 factor = 10.D
              END
              ELSE: Message,'Unknown button pressed'
           ENDCASE
           slice_size(0) = MAX([0.D,mid_point(0)-factor/2.D*side_length])
           slice_size(1) = MAX([0.D,mid_point(1)-factor/2.D*side_length])
           slice_size(2:3) = (slice_size(0:1)+factor*side_length) < 1.D
        END ELSE slice_size = [0.D,0.D,1.D,1.D] 

        windx = slice_size(2)-slice_size(0)
        slider = DOUBLE([(slice_value-windx)>0., slice_value, (slice_value+windx)<1.])
        widget_control, slice_val_w, set_value = slider
;        construct_image
;        show_image
    END

    'INFO_BUTTON': BEGIN
        Print, 'Event for INFO'
        Print, 'sds_num:',sds_num
        Print, 'xy_sl_size:', xy_sl_size
        Print, 'size(grid_info)', size(grid_info)
        Print, 'slice_ori:', slice_ori
        Print, 'slice_size:', slice_size
        Print, 'data_dir:', data_dir
    END
    
    'PS_BUTTON': BEGIN
        IF VERBOSE THEN PRINT, 'Event for dump PS file'
        dump_ps, draw_area, draw_cm, TOP_BASE, $
          [list_str(var_index),'sidelength ='+ $
           STRCOMPRESS(string(slice_size(2)-slice_size(0)))]
    END

    'GIF_BUTTON': BEGIN
        IF VERBOSE THEN PRINT, 'Event for dump GIF image'
        dump_ps, draw_area, draw_cm, TOP_BASE, $
          [list_str(var_index),'sidelength ='+ $
           STRCOMPRESS(string(slice_size(2)-slice_size(0)))], /GIF
    END

    'TIFF_BUTTON': BEGIN
        IF VERBOSE THEN PRINT, 'Event for dump TIFF image'
        dump_tiff, draw_area, draw_cm, TOP_BASE, $
          [list_str(var_index),'sidelength ='+ $
           STRCOMPRESS(string(slice_size(2)-slice_size(0)))]
    END

    'FILE_BUTTON': BEGIN
        If Verbose Then Print, 'Event for FILE'
        widget_control, base2, sensitive = 0
        widget_control, base3, sensitive = 0
        file_name = pickfile(FILTER = '*.hierarchy', $
                             Path=data_dir+'../',   $
                             TITLE  = 'select hierarchy files', $
                             /MUST_EXIST, GET_PATH = data_dir, /NOCONFIRM)
        READ_GRID_INFO, data_dir, file_name, Grid_in_info, $
          list_str, mass_per_particle, divide_by_density, REDSHIFT = z_current
        grid_info = grid_in_info
        widget_control, base2, sensitive = 1
        widget_control, base3, sensitive = 1
    END

    'PEAK_BUTTONS': BEGIN
        If Verbose Then Print, 'Event for PEAK_BUTTONS'
        this_one = WHERE(peak_ids eq FIX(Event.id))
        this_one = this_one(0)
        current_peak = this_one
        FOR i=0,max_peak_nr-1 DO  $
          widget_control, peak_ids(i), set_button=0
        widget_control, peak_ids(current_peak), set_button=1

        print, 'pressed:', this_one
        center = DOUBLE(peaks(0:2, this_one))
        this_sub = where(slice_ori gt 1.e-15)
        slice_ori(this_sub(0)) = center(this_sub(0))
        slice_value = center(this_sub(0))
        windx = DOUBLE(slice_size(2)-slice_size(0))
        os = where(slice_ori lt 1.e-15)
        REPEAT BEGIN
            slice_size(0)   = center(os(0))-windx/2.D
            slice_size(1)   = center(os(1))-windx/2.D
            slice_size(2)   = center(os(0))+windx/2.D
            slice_size(3)   = center(os(1))+windx/2.D
            windx = 0.9*windx
        END UNTIL (min(slice_size ge 0.) and max(slice_size le 1.0))
        slider = [(slice_value-windx)>0., slice_value, (slice_value+windx)<1.]
        widget_control, slice_val_w, set_value = slider

    END


    'READ_PEAKS': BEGIN
        If Verbose Then Print, 'Event for READ_PEAKS'
        widget_control, base2, sensitive = 0
        widget_control, base3, sensitive = 0
        help = ['PICK PEAK FILE', 'Select existing peak file or', $
                'specify a new file name to create your own list of peaks']
         file_name = mypickfile(GROUP=TOP_BASE, FILTER = '*.peaks', $
                             Path=data_dir,   $
                             TITLE  = 'Read peak file', $
                             GET_PATH = peak_in_dir, HELP_TEXT = help)

         widget_control,peak_opt, map = 0
         IF file_name ne '' THEN BEGIN 
             exists = (findfile(file_name))
             print, file_name
             print, 'exists:',exists
             if  (exists(0) eq '') then myeditor, file_name, /NEW $
             ELSE BEGIN
                 current_peak_file = file_name
                 read_peaks, file_name, peaks, N_peaks, max_peak_nr
                 if N_peaks gt 0 then begin 
                     junk_t = STRCOMPRESS(STRING(INDGEN(N_peaks)), /REMOVE_ALL)
                     FOR i =N_peaks-1, max_peak_nr-1 DO $
                       widget_control, peak_ids(i), map = 0, $
                       set_value = '' , set_button=0
                     FOR i =0, N_peaks-1 DO $
                       widget_control, peak_ids(i), map = 1, set_value = junk_t(i) 
                 endif
             ENDELSE
             widget_control,peak_opt, map = 1
         ENDIF           
         widget_control, base2, sensitive = 1
         widget_control, base3, sensitive = 1
     END

    'REREAD_PEAKS':BEGIN
        read_peaks, current_peak_file, peaks, N_peaks, max_peak_nr
    END

    'EDIT_PEAKS': BEGIN
        If Verbose Then Print, 'Event for Edit Peaks'
        IF current_peak_file eq '' THEN BEGIN 
            PRINT, 'INFO: No current peak file'
            PRINT, 'Use "read/create peaks" in File menu !'
        END ELSE BEGIN
            myeditor, current_peak_file
            read_peaks, current_peak_file, peaks, N_peaks, max_peak_nr
            junk_t = STRCOMPRESS(STRING(INDGEN(N_peaks)), /REMOVE_ALL)
            FOR i =0, max_peak_nr-1 DO $
              widget_control, peak_ids(i), map = 0, set_value = '' 
            FOR i =0, N_peaks-1 DO $
              widget_control, peak_ids(i), map = 1, set_value = junk_t(i) 
        ENDELSE
    END


    'FIND_BUTTON': BEGIN
;        widget_control, base2, sensitive = 0
;        widget_control, base3, sensitive = 0
        widget_control, /hourglass

        find_max_value, EVENT.value, grid_info, data_dir,sds_num, xyz_vector
	center = DOUBLE(xyz_vector)
        this_sub = where(slice_ori gt 1.e-15)
        slice_ori(this_sub(0)) = center(this_sub(0))
        slice_value = center(this_sub(0))
        windx = slice_size(2)-slice_size(0)
        os = where(slice_ori lt 1.e-15)
        REPEAT BEGIN
            slice_size(0)   = center(os(0))-windx/2.
            slice_size(1)   = center(os(1))-windx/2.
            slice_size(2)   = center(os(0))+windx/2.
            slice_size(3)   = center(os(1))+windx/2.
            windx = 0.9*windx
        END UNTIL (min(slice_size ge 0.D) and max(slice_size le 1.0D))
        slider = [(slice_value-windx)>0., slice_value, (slice_value+windx)<1.]
        widget_control, slice_val_w, set_value = slider

        if verbose then print, 'max. value at:',  xyz_vector
        WIDGET_CONTROL, base2, sensitive = 1
        WIDGET_CONTROL, base3, sensitive = 1
     END

     'SLICER_BUTTON': BEGIN
        !p.t3d = 1.
         t3d, /reset
         !p.t(*)=0.
        !p.t3d = 0.
         WIDGET_CONTROL, base2,    MAP = 1
         WIDGET_CONTROL, ISO_BASE, MAP = 0
         WIDGET_CONTROL, base3, sensitive = 1
         WIDGET_CONTROL, BASE2, /REALIZE
         WIDGET_CONTROL, wISO,sensitive=1
         WIDGET_CONTROL, wSLICER,sensitive=0
         show_image
     END
         
     'CUBE_DIM': BEGIN
         WIDGET_CONTROL, Event.id, GET_VALUE = cube_dim
     END

     'ISO_SMOOTH_VALUE':BEGIN
         WIDGET_CONTROL, Event.id, GET_VALUE = tmp
         threeD_smoothing=tmp(1)
     END

     'LIGHT_DIRECTION':BEGIN
         aha =  event.value##[1,0,0]
         shade_state.light_direction(0)=-aha(2)
         aha =  event.value##[0,1,0]
         shade_state.light_direction(1)=-aha(2)
         aha =  event.value##[0,0,1]
         shade_state.light_direction(2)=-aha(2)
         display_iso
     END

     'SHADING_MODEL': BEGIN
         EVV = Event.Value
         CASE EVV OF
             "Light Source": BEGIN
                 widget_control, shade_state.base(1), map=0
                 widget_control, shade_state.base(0), map=1,/update
                 shade_state.this_state=0
             END
             "Data set": BEGIN
                 widget_control, shade_state.base(0), map=0
                 widget_control, shade_state.base(1), map=1,/update
                 shade_state.this_state=1
             END 
             ELSE: print, 'unrecognized button event'
         END    
     END 

     'ISO_VAR_LIST': BEGIN
        If Verbose Then Print, 'Event for variable list'
        If Verbose Then Print, 'index:', Event.index
;        shade_state.sds = Event.index+1
        shade_state.sds = Event.index
     END

     'SHADE_MINIMUM': BEGIN
         widget_control, Event.id, get_value=help_val
         shade_state.min=help_val(1)
         widget_control, /hourglass
         compute_em, draw_cm, xy_sl_size,$
           SHADE_EM=shade_state.this_state
         display_iso
     END

     'SHADE_MAXIMUM': BEGIN
         widget_control, Event.id, get_value=help_val
         shade_state.max=help_val(1)
         widget_control, /hourglass
         compute_em, draw_cm, xy_sl_size,$
           SHADE_EM=shade_state.this_state
         display_iso
     END

     'DATA_SHADING_READ_BUTTON': BEGIN
; read data for shading
         update_shading
     END 

     'ISO_BUTTON': BEGIN
;        widget_control, base2, sensitive = 0
         widget_control, vel.baseid, MAP = 0
         widget_control, base3, sensitive = 0
         widget_control, base2,    MAP = 0
         widget_control, ISO_BASE, MAP = 1
         WIDGET_CONTROL, ISO_BASE, /REALIZE
         widget_control, wISO, sensitive=0
         widget_control, wSLICER, sensitive=1
         widget_control, vel.menuid, sensitive=1
         widget_control, /hourglass
; initially have only light shading
         shade_state.this_state = 0
         widget_control, base2, sensitive = 1
     END

     'NO_VEL':BEGIN
         widget_control, vel.baseid, MAP = 0
         widget_control, vel.menuid, set_value='Vel. vectors', $
           SET_UVALUE='VEL_BUTTON'
     END
     'VEL_BUTTON': BEGIN
         widget_control, wSLICER,sensitive=1
         widget_control, wISO,   sensitive=1
         widget_control, ISO_BASE, MAP = 0
         widget_control, vel.baseid, MAP = 1
         widget_control, base3, sensitive = 1
         WIDGET_CONTROL,vel.baseid , /REALIZE
         widget_control, vel.menuid, set_value='hide vector', $
           SET_UVALUE='NO_VEL'
         widget_control, /hourglass
     END

    '4D2_BUTTON': BEGIN
        message,/info,'DUMP_4D2 currently disabled!'

;         If Verbose Then Print, 'Event for dump 4d2'
;         widget_control, base2, sensitive = 0
;         widget_control, base3, sensitive = 0
;         cube_center = slice_ori
;         cube_center(where(slice_ori  lt 1.e-30)) = $
;           [slice_size(0)+0.5*(slice_size(2)-slice_size(0)), $
;            slice_size(1)+0.5*(slice_size(3)-slice_size(1))]
;         cube_length = slice_size(2)-slice_size(0)
;         dump_4d2_data, cube_center, cube_length
;         if verbose then print, 'done with dump 4d2 data'
;         widget_control, base2, sensitive = 1
;         widget_control, base3, sensitive = 1
     END

    'JEANS_BUTTON': BEGIN
        If Verbose Then Print, 'Event for JEANS'
        widget_control, base2, sensitive = 0
        widget_control, base3, sensitive = 0
        
;        construct_slice_data, grid_info, slice_ori, slice_size, data_dir, $
;          sds_num, [xy_sl_size,xy_sl_size], data, SMOOTH = smoothing
        construct_interpolated_slice_data, grid_info, slice_ori, slice_size, data_dir, $
             sds_num, [xy_sl_size,xy_sl_size], data, SMOOTH = smoothing
        data_min =  min(data)
        sec_min  = max([1.e-10,min(data(where(data gt data_min)))])
        data = data > sec_min

        temp_num = 1+where(list_str eq 'Temperature')
         
        construct_slice_data, grid_info, slice_ori, slice_size, data_dir, $
          temp_num(0), [xy_sl_size,xy_sl_size], temp_data, SMOOTH = smoothing
        data_t_min =  min(temp_data)
;        temp_data  = 

        M_jeans = 5.5 + ALOG10( 1/(data^0.5 < 1.e15) *  $
          (temp_data > 1.e-4)^1.5)
        sec_min  =  min(M_jeans(where(m_jeans gt min(M_jeans))))
        TVSCL, M_jeans > sec_min
        min_data = min(M_jeans)
        max_data = max(M_jeans)
        draw_current_cm, draw_cm, xy_sl_size, min_data, max_data

        widget_control, base2, sensitive = 1
        widget_control, base3, sensitive = 1

    END

    'INSIGHT': BEGIN
        message,/info,'INSIGHT currently not available!'
        widget_control, snap_base,sensitive=1
;        insight, /INDEXED_COLOR
    END

    'SNAP_BUTTON': BEGIN
        WIDGET_CONTROL, draw_area, GET_VALUE=win
        WSET, win   
        WINDOW, /PIXMAP, /FREE, xs = !d.x_size, ys=!d.y_size ;Save window
        wset, !d.window
        DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, win] ;Save it
        backing = tvrd()
        backing_num = !d.window
        WIDGET_CONTROL, draw_cm, GET_VALUE=win2
        WSET, win2
        WINDOW, /PIXMAP, /FREE, xs = !d.x_size, ys=!d.y_size ;Save CM window
        wset, !d.window
        DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, win2] ;Save it
        cm_backing = tvrd()
        cm_backing_num = !d.window
        wset, win
        tvlct, r,g,b, /GET
       
        IF (N_ELEMENTS(snap_counter) eq 0) THEN snap_counter=0 $
        ELSE snap_counter = snap_counter + 1

        if event.value eq 'SNAP.separate' THEN BEGIN
            data_name = STRCOMPRESS([image_stack(0).title + '_'+ $
                                     STRCOMPRESS(STRING(snap_counter)), $
                                     image_stack(0).title+'_CM'+ '_'+ $
                                     STRCOMPRESS(STRING(snap_counter))],$
                                    /REMOVE_ALL)
            INSPUT, backing, cm_backing, /IMAGE, REPLACE=2, $
              NAME=[image_stack(0).title,image_stack(0).title+'CM'], $
              RED=r, GREEN=g, BLUE=b, GROUP=TOP_BASE 
        ENDIF ELSE BEGIN
            back_s    =  size(backing)   
            cm_back_s =  size(cm_backing)   
            comb_im   = bytarr(back_s(1), back_s(2)+cm_back_s(2))
            comb_im(0:(back_s(1)-1),cm_back_s(2):(cm_back_s(2)+back_s(2)-1))$
              = backing
            comb_im(0:(back_s(1)-1), 0:(cm_back_s(2)-1)) = cm_backing
            data_name = STRCOMPRESS(image_stack(0).title+ '_'+ $
              STRCOMPRESS(STRING(snap_counter), /REMOVE_ALL))
            INSPUT, comb_im, /IMAGE, REPLACE=2, $
              NAME=data_name, $
              RED=r, GREEN=g, BLUE=b, GROUP=TOP_BASE 
        ENDELSE
    END

    'SET_BUTTON': BEGIN
       mi_ma =  [0.,1.]
       IF  image_stack(0).title NE '' THEN mi_ma = [image_stack(0).min_data,$
                                                    image_stack(0).max_data]
       widget_control, /hourglass
       set_mi_ma, TOP_BASE, mi_ma_base_id, MI_MA=mi_ma
    END

    'REFRESH_BUTTON': BEGIN
       widget_control, /hourglass
       show_image
    END

    'KEEP_BUTTON': BEGIN
        image_stack(0).image_d(0:xy_sl_size-1,0:xy_sl_size-1) = tvrd()
        FOR i=N_images-1,1,-1 DO image_stack(i) = image_stack(i-1)
    END

    'IMAGE_STACK': BEGIN
        im_index   = FIX(EVENT.value)
        if image_stack(im_index).title ne '' THEN BEGIN
            image_stack(0) = image_stack(im_index)
            show_image
        END ELSE print, 'Empty image !'
        widget_control, event.id, set_value = 0
    END

    'DONE_BUTTON': BEGIN
        If Verbose Then Print, 'Event for DONE'
        widget_control, /DESTROY, Event.top
    END
    ELSE:
ENDCASE
;ENDELSE
;ENDIF

END


pro AMRslice, grid_in_info, data_in_dir,$
              reload=reload, base_name=base_name
@TOOLS/common_blocks.inc

  device, retain=2,  BYPASS_TRANSLATION=0 



  if n_elements(grid_in_info) eq 0 then reload=1
  if n_elements(data_in_dir) eq 0 then reload=1
  if n_elements(list_str) eq 0 then reload=1

  if keyword_set(reload) then begin

;      file_name = pickfile(FILTER = '*.hierarchy', $
;                              Path=data_in_dir,   $
;                           TITLE  = 'select hierarchy files', $
;                           /MUST_EXIST, GET_PATH = data_in_dir, /NOCONFIRM)

      if n_elements(base_name) eq 0 then begin
          data_in_dir='/data/enzo/miniqso/P128_G128/PPMdumps_2grid/z=15.50/grids/'
          base_name='miniqso_multi_128_128_2100'
      endif

      READ_GRID_INFO, $
        data_in_dir+base_name, $
        Grid_in_info, $
        list_str, $
        mass_per_particle, $
        divide_by_density, $
        REDSHIFT = z_current
  ENDIF

  

;; mqk: create a pointer array to hold all grids
  ngrids=n_elements(grid_in_info)
  allgrids=ptrarr(ngrids)


  verbose = 1
;; initial values
  slice_value = 0.50
  center = [0.5, slice_value, 0.5]
  slice_ori   = [0,slice_value,0]

  slice_size  = [0.D,0.D,1.D,1.D]
  xy_sl_size = 700
  smoothing  = 0.0
  interpolate_i = 0

  grid_info = grid_in_info
  data_dir  = data_in_dir
  sds_num   = 1
  var_index= 1

  cons_vel = [0., 0.]
  substracting_vel = 0
; prepare Imagestack
  N_images = 4
  Image_structure = {image, $
                     title:          '', $
                     hier_file:      '', $
                     image_d:        DBLARR(1000, 1000), $
                     xy_size:         1,  $
                     min_data:       0.,  $ 
                     max_data:       0.,  $ 
                     slice_size:     DBLARR(4)   $ 
                    }
  
; peak stuff
  max_peak_nr = 4
  N_peaks = 0
  peaks = [0.5,0.5,0.5]
  current_peak_file = ''
  current_peak      = 0
; iso stuff
  cube_dim = 64
  threeD_smoothing = 0.0
  shade_state= {shade_state, main_base:LONG(0), this_state:0, base:LONG([0,0]), $
                sds:0, light_direction:[0.,0.,1.], $
                min:0.,max:1.,min_s:LONG(0),max_s:LONG(0), xa:0.,za:0.}
  
  
  help_image = {image, '', grid_info(0).hier_file, $
                DBLARR(1000, 1000), xy_sl_size, 0., 1., slice_size}
  image_stack = REPLICATE(help_image, N_images)
  
;  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0
  group = 0
  junk   = { CW_PDMENU_S, flags:0, name:'' }
  
  TOP_BASE = WIDGET_BASE(GROUP_LEADER=Group, $
                         ROW=1, $
                         MAP=1, FRAME = 4, $
                         TITLE="JACQUES      -      Enzo's Best friend" , $
                         UVALUE='TOP_BASE', MBAR=wMenuBar, TLB_FRAME_ATTR=1)
  
  parent = WIDGET_BASE(TOP_BASE, $
                       MAP=1, FRAME=4, $
                       UVALUE='BASE2')
  
  BASE2 = WIDGET_BASE(parent, $
                      /COLUMN, MAP=1, $ ; SCR_XSIZE=220,Y_SCROLL_SIZE=XY_SL_SIZE+100, $
                      UVALUE='BASE2')
  
  ISO_BASE = WIDGET_BASE(parent, /COLUMN, MAP=0, UVALUE = 'VALUE_BASE')
  
  
  BASE5 = WIDGET_BASE(TOP_BASE, $
                      COLUMN=1, $
                      MAP=1,  $
                      UVALUE='BASE5')
  
  
  BASE3 = WIDGET_BASE(BASE5, $
                      ROW=1, FRAME = 4, $
                      MAP=1, $
                      UVALUE='BASE3')
  
  BASE4 = WIDGET_BASE(BASE5, $
                      COLUMN=1, $
                      MAP=1, FRAME = 4, $
                      UVALUE='BASE4')
  
  Btns123 = [ $
             'x ', $
             'y ', $
             'z ' ]
  i_ind=  FIX(where(slice_ori gt 0))
  BGROUP10 = CW_BGROUP( BASE2, Btns123, $
                        SET_VALUE = i_ind(0), $
                        ROW=1, $
                        /EXCLUSIVE, $
                        LABEL_TOP='Slice Specifics', $
                        UVALUE='BGROUP10')

; peak stuff
  junk_t = STRARR(max_peak_nr)
  junk_t(*) = ''
  hold_base = widget_base(base2, ROW=3, Frame = 2)
  junk = widget_label(hold_base, value='Peaks')
  n = n_elements(junk_t)
  peak_ids = lonarr(n)
  junk_base = lonarr(n)
  for i = 0, n-1 do begin
     junk_base(i) = WIDGET_BASE(hold_base,/ROW,/NONEXCLUSIVE,map=0 )
     peak_ids(i) = WIDGET_BUTTON(junk_base(i), $
                                 value=junk_t(i), UVALUE='PEAK_BUTTONS',$
                                 /DYNAMIC_RESIZE )
  endfor
  junk_t = STRCOMPRESS(STRING(INDGEN(max_peak_nr)), /REMOVE_ALL)
  FOR i =0, N_peaks-1 DO $
   widget_control, peak_ids(i), map = 1, set_value = junk_t(i), set_button=0
  
  peak_opt    = Widget_base(hold_base, /ROW, MAP = 0)
  peak_update = WIDGET_BUTTON(peak_opt, $
                              VALUE='Reread', UVALUE='REREAD_PEAKS') 
  peak_edit = WIDGET_BUTTON(peak_opt, VALUE='Edit', UVALUE='EDIT_PEAKS') 
  
  junk = CW_MYFSLIDER( BASE2, $
                       MAXIMUM=1.00000, $
                       MINIMUM=0.00001, $
                       TITLE='slice coord', $
                       UVALUE='SLICE_VAL', $
                       FORMAT='(F13.9)',    $
                       VALUE=slice_value, /EDIT)
  
  slice_val_w = junk
  
  Btns123 = [ $
             'none', $
             'lin.', $
             '5th' ]

  BGROUP10 = CW_BGROUP( BASE2, Btns123, $
                        SET_VALUE = interpolate_i, $
                        ROW=1, $
                        /EXCLUSIVE, /NO_RELEASE,$
                        LABEL_TOP='interpolation', $
                        UVALUE='BGROUP_INTER')
  
  
  SLICE_TEXT = WIDGET_TEXT(BASE2, $
                           YSIZE = 8, UVALUE='SLICE_TEXT')
  
  
  VAR_LIST = WIDGET_DROPLIST( BASE2, VALUE = List_str, $
                              UVALUE='VAR_LIST'   )
  widget_control,var_list, set_droplist_select = var_index
  
  VEL_BASE = WIDGET_BASE(BASE2, /COLUMN, MAP=0, UVALUE = 'VEL_BASE',FRAME=1)
  
  
  DRAW_BUTTON = WIDGET_BUTTON( BASE3, $
                               UVALUE='DRAW_BUTTON', $
                               VALUE='DRAW')
  
  desc = [{CW_PDMENU_S, 1, name: 'ZOOM IN'},  $
          {CW_PDMENU_S, 0, name: 'by 2'},  $
          {CW_PDMENU_S, 0, name: 'by 4'},  $
          {CW_PDMENU_S, 0, name: 'by 10'},  $
          {CW_PDMENU_S, 0, name: 'by 100'},  $
          {CW_PDMENU_S, 0, name: 'by 1000'},  $
          {CW_PDMENU_S, 0, name: 'pick'} ]
  
  ZOOM_IN_BUTTON = CW_PDMENU( BASE3, desc, /RETURN_FULL_NAME, $
                              UVALUE='ZOOM_IN_BUTTON')
  
  desc = [{CW_PDMENU_S, 1, name: 'ZOOM OUT'},  $
          {CW_PDMENU_S, 0, name: 'by 2'},  $
          {CW_PDMENU_S, 0, name: 'by 5'},  $
          {CW_PDMENU_S, 0, name: 'by 10'},  $
          {CW_PDMENU_S, 0, name: 'top'} ]
  
  ZOOM_OUT_BUTTON = CW_PDMENU( BASE3, desc, /RETURN_FULL_NAME, $
                               UVALUE='ZOOM_OUT_BUTTON')
  
  desc = [{CW_PDMENU_S, 1, name: 'FIND MAX'},  $
          {CW_PDMENU_S, 0, name: 'from the finest levels'},  $
          {CW_PDMENU_S, 0, name: 'from all'} ]
  FIND_BUTTON = CW_PDMENU( BASE3, desc, /RETURN_FULL_NAME, $
                           UVALUE='FIND_BUTTON')
  
  desc = [{CW_PDMENU_S, 1, name: 'SNAP'},  $
          {CW_PDMENU_S, 0, name: 'separate'},  $
          {CW_PDMENU_S, 0, name: 'combined'} ]
  SNAP_BASE =  CW_PDMENU( BASE3, desc, /RETURN_FULL_NAME,$
                          UVALUE='SNAP_BUTTON')
  widget_control, snap_base,sensitive=0
  
  SET_BUTTON = WIDGET_BUTTON( BASE3, $
                              UVALUE='SET_BUTTON', $
                              VALUE='MI MA')
  
  REFRESH_BUTTON = WIDGET_BUTTON( BASE3, $
                                  UVALUE='REFRESH_BUTTON', $
                                  VALUE='REFRESH')
  
  KEEP_BUTTON = WIDGET_BUTTON( BASE3, $
                               UVALUE='KEEP_BUTTON', $
                               VALUE='KEEP')
  
  junk_t = STRCOMPRESS(STRING(INDGEN(N_images)), /REMOVE_ALL)
  junk = CW_BGROUP(base3, column = N_images, /EXCLUSIVE,  $
                   junk_t, UVALUE = 'IMAGE_STACK',  $
                   BUTTON_UVALUE = b_uvalues, SET_VALUE = 0)
  
  
  DRAW_AREA = WIDGET_DRAW(BASE4, $
                          XSIZE = xy_sl_size, $
                          YSIZE = xy_sl_size,   $
                          /BUTTON, /MOTION, RETAIN=2, $
                          UVALUE = 'draw_area')
  
  DRAW_CM = WIDGET_DRAW(BASE4, $
                        XSIZE = xy_sl_size, $
                        YSIZE = 80,  $
                        /BUTTON, /MOTION, RETAIN=2, /TRACKING_EVENTS,$
                        UVALUE = 'draw_cm')
  
  wBase = WIDGET_BASE(TOP_BASE, /COLUMN)
  
  wFileMenu = WIDGET_BUTTON(wMenuBar, VALUE='File', /MENU)
  wLoadItem = WIDGET_BUTTON(wFileMenu, VALUE='Open', UVALUE='FILE_BUTTON')
  wLoadItem = WIDGET_BUTTON(wFileMenu, VALUE='Read/Create peaks', UVALUE='READ_PEAKS')   
  wLoadItem = WIDGET_BUTTON(wFileMenu, VALUE='Save as PS', UVALUE='PS_BUTTON')
  wLoadItem = WIDGET_BUTTON(wFileMenu, VALUE='Save as Gif', UVALUE='GIF_BUTTON')
  wLoadItem = WIDGET_BUTTON(wFileMenu, VALUE='Save as TIFF', UVALUE='TIFF_BUTTON')
  wLoadItem = WIDGET_BUTTON(wFileMenu, $
                            VALUE='Dump hdf cube (4d2)', UVALUE='4D2_BUTTON')
  wExitItem = WIDGET_BUTTON(wFileMenu, VALUE='Exit', UVALUE='DONE_BUTTON')
  
  wMenu = WIDGET_BUTTON(wMenuBar, VALUE='View', /MENU)
  wISO  = WIDGET_BUTTON(wMenu, VALUE='Isosurface', UVALUE='ISO_BUTTON')
  wSLICER  = WIDGET_BUTTON(wMenu, VALUE='Slicer', UVALUE='SLICER_BUTTON')
  widget_control, wSLICER,sensitive=0
  wVEL  = WIDGET_BUTTON(wMenu, VALUE='Vel. vectors', UVALUE='VEL_BUTTON')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Grid chart', UVALUE='SHOW_CHART')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Grid boxes', UVALUE='SHOW_BOXES')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Peaks', UVALUE='SHOW_PEAKS')   
  wItem = WIDGET_BUTTON(wMenu, VALUE='Color Map', UVALUE='CM_BUTTON')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Data info', UVALUE='INFO_BUTTON')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Hierarchy file', UVALUE='VIEW_HIERARCHY')
  wItem = WIDGET_BUTTON(wMenu, $
                        VALUE='Data parameter file',UVALUE='VIEW_PARAMETER')
  
  wMenu = WIDGET_BUTTON(wMenuBar, VALUE='Tools', /MENU)
  wItem = WIDGET_BUTTON(wMenu, VALUE='Annotate CM', UVALUE='AN_CM')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Annotate Graph', UVALUE='AN_DRAW')   
  wItem = WIDGET_BUTTON(wMenu, VALUE='Examine Text file', UVALUE='EDITOR')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Calculator', UVALUE='Calculator')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Xpalette', UVALUE='XPALETTE')
  wItem = WIDGET_BUTTON(wMenu, VALUE='Insight', UVALUE='INSIGHT')
  
  
  wHelpMenu = WIDGET_BUTTON(wMenuBar, VALUE='Help', /HELP, /MENU)
  wInfoItem = WIDGET_BUTTON(wHelpMenu, VALUE='Help', UVALUE='HELP')
  wInfoItem = WIDGET_BUTTON(wHelpMenu, VALUE='About...', UVALUE='About')
  
; ------------------------------ ISO SURFACE STUFF ---------------------------
  if n_elements(min_v) eq 0 and n_elements(is_value) ne 0 then $
    undefine,is_value

  IS_VAL = CW_myFSLIDER( ISO_BASE, $
                         MAXIMUM=max_v, $
                         MINIMUM=min_v,$
                         UVALUE = 'IS value', $
                         TITLE='IS value', $
                         VALUE=is_value, /EDIT)
  BTNS = ["Compute", 'Weigh','(Re)Read','IDL SLICER']
  BGROUP = CW_BGROUP(ISO_BASE, BTNS,  $
                     COLUMN=2, $
                     UVALUE='BTNS', $
                     /RETURN_NAME)
  
; variable shading models
  isojunk = widget_base(ISO_BASE,/COLUMN,FRAME=3)
  junk = CW_BGROUP(isojunk, /ROW, /EXCLUSIVE, /NO_REL, LABEL_TOP='Shading',$
                   ['Light Source', 'Data set'], $
                   UVALUE='SHADING_MODEL', SET_VALUE=0, $
                   /RETURN_NAME)
  shade_state.main_base = WIDGET_BASE(isojunk)
  for i=0,1 DO shade_state.base(i) = WIDGET_BASE(shade_state.main_base,$
                                                 /COLUMN,/FRAME)
  fbase = WIDGET_BASE(shade_state.base(0))
  fbasetwo = WIDGET_BASE(shade_state.base(1),/column)
  
  temp_index =  (where(STRPOS(List_str,'emper') ne -1))
  IF N_ELEMENTS(temp_index) gt 0 THEN temp_index = temp_index(0) $
  ELSE temp_index = 0
  
; light model base
  help = indgen(5)+1
  while max(help) lt !D.table_size-3 DO help=(1.1*help)
  help = Round(help/1.1)
  junk = CW_ARCBALL(fbase, LABEL='light direction', $
                    size=150, $
                    UVALUE = 'LIGHT_DIRECTION', RETAIN = 2, /UPDATE)
  
  junk = WIDGET_DROPLIST(fbasetwo, VALUE = List_str, $
                         UVALUE='ISO_VAR_LIST')
  widget_control, junk, set_droplist_select=temp_index
;  shade_state.sds = temp_index+1
  shade_state.sds = temp_index
  shade_state.min_s = CW_myFSLIDER(fbasetwo, $
                                   MAXIMUM=5., $
                                   MINIMUM=0.1,$
                                   UVALUE = 'SHADE_MINIMUM', $
                                   TITLE='shade min', $
                                   VALUE=1.)
  nextbase = WIDGET_BASE(fbasetwo,/COLUMN)
  shade_state.max_s = CW_myFSLIDER(nextbase,    $
                                   MAXIMUM=5., $
                                   MINIMUM=0.1,$
                                   UVALUE = 'SHADE_MAXIMUM', $
                                   TITLE='shade max', $
                                   VALUE=1.)
  nextbase = WIDGET_BASE(nextbase,/COLUMN)
  BTNS = 'Read'
  BGROUP = WIDGET_BUTTON(nextbase, VALUE=BTNS,  $
                         UVALUE='DATA_SHADING_READ_BUTTON')
  
  ori = CW_ORIENT(ISO_BASE, FRAME = 3, TITLE='ORIENTATION', $
                  UVALUE = 'ORI')
  
  ISO_SMOOTH_VALUE_W = CW_myFSLIDER(ISO_BASE, FORMAT='(F3.1)', $
                                    MAXIMUM=10., $
                                    MINIMUM=0.,$
                                    UVALUE = 'ISO_SMOOTH_VALUE', $
                                    TITLE='smooth (0=off)', $
                                    VALUE=threeD_smoothing, /EDIT)
  
  junk = CW_FIELD(ISO_BASE,VALUE=cube_dim, $
                  COLUMN=1, $
                  INTEGER=1, $
                  RETURN_EVENTS=1, $
                  TITLE='Cube size', $
                  UVALUE='CUBE_DIM')
  
  
  
; -----------------------------------------------------------------------
  
; --------------------- VELOCITY STUFF ----------------------------------
  def_skip = FIX(xy_sl_size/100) > 5
  def_amp  = 2.
  vel_structur = {velocity, skip:0, amplify:0., $
                  color : 0, baseid:0, menuid:0}
  vel ={velocity, def_skip, def_amp, $
        !d.table_size-1, VEL_BASE, wVEL}
  junk = WIDGET_LABEL(VEL_BASE,VALUE='velocity vectors')
  junk = WIDGET_SLIDER( VEL_BASE, $
                        /DRAG,     $
                        MAXIMUM=25, $
                        MINIMUM=5,$
                        UVALUE = 'Vel_skip', $
                        TITLE='Skip factor', $
                        VALUE=vel.skip)
  junk = CW_myFSLIDER( VEL_BASE, $
                       MAXIMUM=10., $
                       MINIMUM=0.1,$
                       UVALUE = 'Vel_amp', $
                       TITLE='amplify ', $
                       VALUE=def_amp)
  junk = CW_MYCLR_INDEX(VEL_BASE, start_color = 0, NCOLORS=!D.table_size, $
                        xsize=100,UVALUE="VEL_COL_SEL", FRAME = 1, LABEL="Color")
  BTNS = ['Read','Substract','Plot', 'Reset', 'Hide']
  BGROUP = CW_BGROUP(VEL_BASE, BTNS,  $
                     COLUMN=3, $
                     UVALUE='VELBTNS', $
                     /RETURN_NAME)
; -----------------------------------------------------------------------
  
  widget_control, shade_state.base(1), MAP=0
  
  WIDGET_CONTROL, TOP_BASE, /REALIZE

;; set colormap
  device, bypass_translation=0
  loadct,4

; now get the window number for the draw area
  draw_current_cm, draw_cm, xy_sl_size, 0., 1.
  update_slice_text, slice_text, slice_size, ''
  widget_control, get_value = win_index, draw_area
; make it the current graphics window
  wset, win_index
;  DEVICE, BYPASS=0

  XMANAGER, 'TOP_BASE', TOP_BASE, EVENT_HANDLER='TOP_BASE_Event'
  

;; free allgrids pointers
  ptr_free,allgrids
  
END
; .compile TOOLS/AMRslice.pro
