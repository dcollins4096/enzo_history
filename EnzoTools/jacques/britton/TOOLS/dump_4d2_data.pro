;pro dump_4d2_data, grid_info, cube_center, cube_length, data_dir, cube_size

pro  dump_the_cube, grid_info, list_str, cube_center, cube_length, data_dir, selected, cube_size, file_name

; check first whether file path and name are valid
get_lun, unit
openw, unit, file_name, /DELETE, ERROR = myerr
;printf,unit, '?'
print, 'myerr:', myerr
if (myerr NE 0) then BEGIN
    print, 'there is a problem with your filename :', file_name
    PRINT, !ERR_string
    RETURN
ENDIF
close, unit
free_lun, unit
  s_selected = size(selected)
  entries = s_selected(2)
  take_log = selected(2,*)
  take_dex = selected(1,*)
  dump_it  = selected(0,*)
  these    = take_log or take_dex or dump_it
  print, these
  help     = indgen(entries)
  d_sds    = help(where(these gt 0))
  print, 'd_sds:', d_sds
  s_d_sds  = size(d_sds)
; reorder variables such that vector quantities are first
  i=-1
  found_one = 0 
  order    = d_sds
  
  FOR i=0,s_d_sds(3)-1 DO BEGIN
      check = STRPOS(list_str(d_sds(i)),'elocity') 
      if ((check gt -1) and (found_one eq 0)) THEN BEGIN
          found_one = 1
          order = shift(d_sds,-i)
          print, 'reorderd data set to have vector fields first!'
      ENDIF
  ENDFOR

;  take_log= [0,0, 0, 1, 1] 
  
  dim_grid_info = size(grid_info)
  num_of_grids  = dim_grid_info(1)
  
; min left edge and scaling
  min_left  = $
    [cube_center(0),cube_center(1),cube_center(2)]-cube_length/2. > 0.
  max_right = $
    [cube_center(0),cube_center(1),cube_center(2)]+cube_length/2. < 1.
  
  scale_up    = min(1./(max_right-min_left))
  
; determine grids that have sufficient res and are at least partially
; in the cube 
  u_g = grid_in_cube(grid_info, cube_center, cube_length, cube_size)
  s_u_g = size(u_g)
  if (s_u_g(0) eq 0) THEN RETURN
  
  index_range = u_g.End_index-u_g.Start_index+1
  delta_distance = (u_g.Right_edge - $
                    u_g.Left_edge)/DOUBLE(index_range)
  dim_grid_info = size(u_g)
  num_of_grids  = dim_grid_info(1)
  print, size(u_g)
;  print,  'max:', min_left, max_right, scale_up
; set up data array
  cube_data = dblarr(cube_size, cube_size, cube_size)
  
  s_order = size(order)
  FOR sds_count=0, s_order(1)-1 DO BEGIN ; iterate over number of SDS
      sds_num = order(sds_count)
      cube_data(*) = 0.
      FOR i=0,num_of_grids-1 DO BEGIN
; check this
          delta_dist   = delta_distance(*,i)
          data_points  = index_range(*,i)  
          i_l = CEIL((min_left - u_g(i).Left_edge)/delta_dist >0 $
                     < (data_points-1))
          i_r   = FLOOR(((max_right -u_g(i).Left_edge)/ $
                         delta_dist - 1) < (data_points-1) ) > 0
;        print, 'index', i, i_r, i_l
          Left_point = FIX(DOUBLE(cube_size) *scale_up * $
                           ((u_g(i).Left_edge+ $
                             (i_l)*delta_dist)-min_left) > 0.)
          Right_point= Left_point + $
            FIX(DOUBLE(cube_size) *scale_up * $
                double(i_r-i_l+1)*delta_dist) < (cube_size-1)
          sub_cube_size = Right_point-Left_point
          print, 'left and right',  Left_point  , Right_point
          if ((min(i_r-i_l) ge 3) and (min(sub_cube_size) ge 3)) THEN BEGIN
              read_hdf, data_dir+u_g(i).baryon_file, sds_num+1, temp, dummy_l
              label = list_str(sds_num)
              temp = REFORM(temp(i_l(0):i_r(0),i_l(1):i_r(1),i_l(2):i_r(2)))
              
              print, 'sub_cube_size', sub_cube_size
              sub_cube= CONGRID(temp, $
                sub_cube_size(0)+1,sub_cube_size(1)+1,sub_cube_size(2)+1)
              cube_data(Left_point(0):Right_point(0),      $
                        Left_point(1):Right_point(1),      $
                        Left_point(2):Right_point(2)       ) = sub_cube
          ENDIF
      ENDFOR
      
      if take_log(order(sds_count)) then BEGIN
          print, 'take log10 of ' + label +' and write it to ' +file_name
          label = 'logof'+ label
          sec_min_data = min(cube_data(where(cube_data ne min(cube_data))))
          cube_data = alog10(cube_data > sec_min_data)
      ENDIF

      if take_dex(order(sds_count)) then BEGIN
          print, 'take dex of ' + label +' and write it to ' +file_name
          label = 'dex of '+ label
          cube_data = 10^(cube_data)  
      ENDIF
      
; append to file
print,' adding ',label, ' to ',file_name
  if sds_count eq 0 then hdfu = HDF_OPEN(file_name, /ALL)
  HDF_DFSD_SETINFO, DIMS = $ 
    [cube_size,cube_size,cube_size], /CLEAR, /DOUBLE, label = label
  HDF_DFSD_ADDDATA, file_name, DOUBLE(cube_data)
      
      
ENDFOR
HDF_CLOSE, hdfu
RETURN
END


pro dump_4d2_data, cube_center, cube_length
@TOOLS/common_blocks.inc

boldfont = '-*-fixed-bold-*-normal--*-*-*-*-*-*-*-*'
    num_dim = 3
; rough check of input data
    IF (N_ELEMENTS(grid_info) lt 1) THEN BEGIN
        print, 'construct_3D_data: no grid information in grid_info:',$
          grid_info
        STOP
    ENDIF

; Initial values
; filename for output file
file_name = 'for4d2'
file_path = './'
cube_size = 128
    
fbase = widget_base(title='Dump Data Cube', /ROW)  
fbase1 = widget_base(fbase, /COLUMN, /FRAME)  
templabel   = widget_label(fbase1, $
                           value = 'Select data sets to dump', $
                           FONT = boldfont)
fbase2 = widget_base(fbase1, /COLUMN, /FRAME)  
fbase3 = widget_base(fbase, /COLUMN, /FRAME)  

s = size(list_str(where(list_str ne ''))) 
entries = s(3)
options = 3
junk = LONARR(entries)
basejunk = LONARR(entries)
FOR i=0,entries-1 DO BEGIN
    b_uvalues = indgen(entries)+i*options
    basejunk(i) = widget_base(fbase2, /ROW) 
    templabel   = widget_label(basejunk(i), $
                   value = STRING(list_str(i), FORMAT = '(A20)'), $
                   FONT = boldfont)
    junk(i) = CW_BGROUP(basejunk(i), column = entries, /NONEXCLUSIVE,  $
    ['dump','dex','log10'], UVALUE = STRING(i), $
      BUTTON_UVALUE = b_uvalues)

ENDFOR

dimwid  = CW_FIELD(FBASE3,  VALUE=cube_size, TITLE ='Cube dimensions', $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
                       UVALUE = 'GRIDDIM', /ALL_EVENTS, /INTEGER)

fbase4 = widget_base(fbase3, /COLUMN, FRAME = 2)  
filewid = CW_FIELD(FBASE4,  VALUE=file_name, TITLE ='File name', $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
                       UVALUE = 'FILENAME', /STRING)
browseit = WIDGET_BUTTON( FBASE4, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='BROWSE', $
      VALUE='BROWSE')


dumpit = WIDGET_BUTTON( FBASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='DUMP', $
      VALUE='DUMP')

didit = WIDGET_BUTTON( FBASE3, $
      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
      UVALUE='didit', $
      VALUE='DONE')

WIDGET_CONTROL, fbase, /REALIZE
;WSHOW, win

WHILE 1 DO BEGIN                ;Internal event loop   
    ev = WIDGET_EVENT(fbase)
    WIDGET_CONTROL,Ev.id,GET_UVALUE=Event
    CASE event OF 
        'didit'  : BEGIN
            WIDGET_CONTROL, fbase, /DESTROY   
            RETURN
        END
        'GRIDDIM': BEGIN 
            WIDGET_CONTROL, Ev.id, GET_VALUE = cube_size
            print, 'cube_size:', cube_size
            IF cube_size gt 512 THEN BEGIN 
                print, cube_size, ' is huge !'
            ENDIF
            IF cube_size lt 32 THEN BEGIN 
                print, cube_size, ' is very small !'
            ENDIF
        END
        'FILENAME' : BEGIN
            WIDGET_CONTROL, Ev.id, GET_VALUE = file_name
        END            
        'BROWSE' : BEGIN
            new_name = PICKFILE(FILE=file_name, GET_PATH = new_path, $
                                 GROUP = FBASE, /WRITE, PATH = file_path, $
                                 FILTER = '*.hdf')
            if strlen(new_name) gt 0 THEN BEGIN
                WIDGET_CONTROL, filewid, SET_VALUE=new_name
                file_path = new_path
                file_name = new_name
            ENDIF
        END
        'DUMP' : BEGIN
;
; Do it
;
  selected = BYTARR(options, entries)
  for i=0, entries-1 DO BEGIN
      WIDGET_CONTROL,  GET_VALUE = theval , junk(i)
      print, i, theval
      selected(*,i) = theval                
  ENDFOR
;  print, 'selected:', selected

  IF TOTAL(selected) gt 0 THEN BEGIN
      widget_control, Fbase, sensitive = 0
      widget_control, /hourglass
      dump_the_cube, grid_info, list_str, cube_center, $
        cube_length, data_dir, selected, cube_size, file_name 
      widget_control, Fbase, sensitive = 1
  ENDIF ELSE BEGIN
      print, 'You did not select any fields to dump !!'
     END 
  END
        ELSE : 
    ENDCASE 
    
ENDWHILE

  END
