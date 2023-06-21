pro construct_interpolated_slice_data, $
              slice_data, SECMIN=SECMIN
; PLOT DERIVED QUANTITIES such as angular momentum, Entropy and Pressure ...
@common_blocks.inc
; using grid_info read in data (found in data_dir directory) of all
; necessary grids for a slice ; specified by the 3 component vector
; slice_ori. slice_size has 4 components specifying the left and right
; edge of the slice
; note that slice_ori has to have two elements set to zero since we
; only do orthogonal slices.
; 
interpolate =  interpolate_i
; should I log it ?
logit = 1
IF (strpos(list_str(var_index), '-velo') GE 0) THEN logit = 0 

angular = 0 & & pressure = angular & & entropy = angular 
CASE 1 OF
   (strpos(list_str(var_index), 'Angul') GE 0):     angular =  1
   (strpos(list_str(var_index), 'Entropy') GE 0):   entropy =  1
   (strpos(list_str(var_index), 'Pressure') GE 0):  pressure = 1
   ELSE:
ENDCASE
  IF angular THEN BEGIN
     IF verbose THEN print, 'con_angular_momentum...: Using center:', center
  ENDIF                         ; if angular
;ok find variables:
  FOR i=0,N_ELEMENTS(list_str)-1 DO BEGIN
     l =  list_str(i)
     CASE 1 OF
        (strpos(l, 'x-vel') GE 0): xvel_sds =  i+1
        (strpos(l, 'y-vel') GE 0): yvel_sds =  i+1
        (strpos(l, 'z-vel') GE 0): zvel_sds =  i+1
        (strpos(l, 'Density') GE 0): D_sds =  i+1
        (strpos(l, 'Temperature') GE 0): T_sds =  i+1
        ELSE:
     ENDCASE
  ENDFOR
  
  num_dim = 3
; rough check of input data
  slice_data = DBLARR(xy_sl_size, xy_sl_size)
;  slice_data(*,*) = 0.
  axes = (slice_ori gt 1.e-30)
  const_sub =  where(axes eq max(axes))
  const_sub = const_sub(0)
  other_subs = where(axes eq 0)
  
  IF (Fix(TOTAL(axes)) gt 1) then BEGIN
      print, 'construct_angular_momentum_slice_data: slice_ori:', $
        slice_ori,' not allowed !'
      STOP
  ENDIF
  IF (N_ELEMENTS(grid_info) lt 1) THEN BEGIN
      print, 'construct_angular_momentum_slice_data: no grid information in grid_info:',$
        grid_info
      STOP
  ENDIF
  dim_grid_info = size(grid_info)
  num_of_grids  = dim_grid_info(1)

; min left edge and scaling
  min_left  = DOUBLE([slice_size(0),slice_size(1)])
  max_right = DOUBLE([slice_size(2),slice_size(3)]) 
  print, min_left
  scale_up    = min(1.D/(max_right-min_left))

  index_range    = grid_info.End_index-grid_info.Start_index+1
  delta_distance = DOUBLE(grid_info.Right_edge - $
                    grid_info.Left_edge)/DOUBLE(index_range)

  IF angular THEN BEGIN 
; find velocity of center:
     ml =  -1
     this_grid =  1
     FOR i =0,num_OF_grids-1 DO BEGIN
        if  (((grid_info(i).left_edge(0) - center(0)) GT 0.) AND $
             ((grid_info(i).left_edge(1) - center(1)) GT 0.) AND $
             ((grid_info(i).left_edge(2) - center(2)) GT 0.) AND $
             ((-grid_info(i).right_edge(0)+center(0)) LT 0.) AND $
             ((-grid_info(i).right_edge(1)+center(1)) LT 0.) AND $
             ((-grid_info(i).right_edge(2)+center(2)) LT 0.)) THEN $
         IF (grid_info(i).level Gt ml) THEN BEGIN 
           this_grid = i 
           ml =  grid_info(i).level 
        ENDIF
     ENDFOR
     print, 'this_grid:', this_grid, grid_info(0:5).level
; read that grid
     read_hdf, data_dir+grid_info(this_grid).baryon_file, xvel_sds, tempx
     read_hdf, data_dir+grid_info(this_grid).baryon_file, yvel_sds, tempy
     read_hdf, data_dir+grid_info(this_grid).baryon_file, zvel_sds, tempz
     ix =  ROUND((center(0)-grid_info(this_grid).left_edge(0))*delta_distance(this_grid))
     iy =  ROUND((center(1)-grid_info(this_grid).left_edge(1))*delta_distance(this_grid))
     iz =  ROUND((center(2)-grid_info(this_grid).left_edge(2))*delta_distance(this_grid))
     cenXvel =  tempx(ix,iy,iz)
     cenYvel =  tempy(ix,iy,iz)
     cenZvel =  tempz(ix,iy,iz)
     IF verbose THEN print, 'central velocity:', cenXvel, cenYvel, cenZvel
  ENDIF                         ; if angular 

; use only grids that will show up on the image:
  dxts = DOUBLE([xy_sl_size,xy_sl_size])/(max_right-min_left)
; determine grids that are at least partially in the slice
  grid_flag      = grid_in_slice(grid_info, DOUBLE(slice_ori), DOUBLE(slice_size))
  index_range    = grid_info.End_index-grid_info.Start_index+1
  delta_distance = DOUBLE(grid_info.Right_edge - $
                    grid_info.Left_edge)/DOUBLE(index_range)
  del_pix =   DOUBLE(max([xy_sl_size,xy_sl_size]))*delta_distance/(max_right(0)-min_left(0))*DOUBLE(index_range)
; initialize minimum of data 
  min_data = 1.e30
  u_g = grid_info(where(( $
        del_pix(other_subs(0),*) GE .5) AND (grid_flag gt 0)  )) 

; print,  del_pix
  dim_grid_info = size(u_g)
  num_of_grids  = dim_grid_info(1)

  fully_contained_above =  0
  FOR i =0, num_of_grids-1 DO BEGIN
     if ((u_g(i).Left_edge(other_subs(0)) le (min_left(0)))$
        and  (u_g(i).Left_edge(other_subs(1)) le (min_left(1))) $
         and (u_g(i).Right_edge(other_subs(0)) ge (max_right(0))) $
         and (u_g(i).Right_edge(other_subs(1)) ge (max_right(1)))) THEN $
      fully_contained_above = i
  ENDFOR
  
  print, 'slice is fully contained in grids after index:',  $
   fully_contained_above
  
  index_range = u_g.End_index-u_g.Start_index+1
  delta_distance = (u_g.Right_edge - $
                    u_g.Left_edge)/DOUBLE(index_range)
;  print, size(u_g)
;  print,  'max:', min_left, max_right, scale_up
  x_stuetz =  DBLARR(1000000)
  y_stuetz =  DBLARR(1000000)
  z_stuetz =  DBLARR(1000000)
  dx_stuetz = DBLARR(1000000)
     
  n_stuetz = 0L
  this_level=0
  le_o    = u_g(0:num_of_grids-1).left_edge(other_subs(0))
  le_t    = u_g(0:num_of_grids-1).left_edge(other_subs(1))
  ri_o    = u_g(0:num_of_grids-1).right_edge(other_subs(0))
  ri_t    = u_g(0:num_of_grids-1).right_edge(other_subs(1))
  
  slice_ori_l                = dblarr(3)
  slice_ori_l(const_sub)     = slice_ori(const_sub)
  slice_ori_l(other_subs) = min_left
  slice_ori_r                = slice_ori_l 
  slice_ori_r(other_subs) = max_right
; 
  Num_pix_x = xy_sl_size
  Num_pix_y = xy_sl_size
  
  
  FOR i=fully_contained_above, num_of_grids-1 DO BEGIN
; check this
     delta_dist   = delta_distance(*,i)
     data_points  = index_range(*,i)  
; yep check again !!
     i_s = ROUND((slice_ori_l(const_sub) - 0.5*delta_dist(const_sub) -  $
                  u_g(i).Left_edge(const_sub))/delta_dist(const_sub) >0 $
                 < (data_points(const_sub)-1))
     i_l =   (FLOOR((slice_ori_l(other_subs) - $
                  u_g(i).Left_edge(other_subs)) $
                 /delta_dist(other_subs)) -1) > 0
     i_r   = ((CEIL((slice_ori_r(other_subs) - $
                     u_g(i).Left_edge(other_subs))/ $
                    delta_dist(other_subs)) ) < $
              (data_points(other_subs)-1))  > 0
     if ((min(i_r-i_l) ge 1)  ) THEN BEGIN
        CASE 1 OF
           angular: BEGIN
              read_hdf, data_dir+u_g(i).baryon_file, xvel_sds, tempx
              read_hdf, data_dir+u_g(i).baryon_file, yvel_sds, tempy
              read_hdf, data_dir+u_g(i).baryon_file, zvel_sds, tempz
              s = size(tempx)
              if (min([s(3),s(2),s(1)]) ge 2) THEN BEGIN
                 CASE const_sub OF 
                    0: tempx = REFORM(tempx(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                    1: tempx = REFORM(tempx(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                    2: tempx = REFORM(tempx(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                    ELSE:
                 ENDCASE
                 CASE const_sub OF 
                    0: tempy = REFORM(tempy(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                    1: tempy = REFORM(tempy(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                    2: tempy = REFORM(tempy(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                    ELSE:
                 ENDCASE
                 CASE const_sub OF 
                    0: tempz = REFORM(tempz(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                    1: tempz = REFORM(tempz(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                    2: tempz = REFORM(tempz(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                    ELSE:
                 ENDCASE
              END
           END
           (entropy OR pressure): BEGIN 
              read_hdf, data_dir+u_g(i).baryon_file, D_sds, tempD
              read_hdf, data_dir+u_g(i).baryon_file, T_sds, tempT
              s = size(tempD)
              if (min([s(3),s(2),s(1)]) ge 2) THEN BEGIN
                 CASE const_sub OF 
                    0: tempD = REFORM(tempD(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                    1: tempD = REFORM(tempD(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                    2: tempD = REFORM(tempD(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                    ELSE:
                 ENDCASE
                 CASE const_sub OF 
                    0: tempT = REFORM(tempT(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                    1: tempT = REFORM(tempT(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                    2: tempT = REFORM(tempT(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                    ELSE:
                 ENDCASE
              ENDIF
           END                  ; entropy or pressure
           ELSE: BEGIN          ; 'normal' quantity
              read_hdf, data_dir+u_g(i).baryon_file, sds_num, temp
              s = size(temp)
              if (min([s(3),s(2),s(1)]) ge 2) THEN BEGIN
                 CASE const_sub OF 
                    0: temp = REFORM(temp(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                    1: temp = REFORM(temp(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                    2: temp = REFORM(temp(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                    ELSE: 
                 ENDCASE
              END  ELSE BEGIN print, 'warning: very small grid '
              ENDELSE
              IF  (divide_by_density(var_index) gt 0) THEN BEGIN
                 read_hdf, data_dir+u_g(i).baryon_file, 1, tempD
                 s = size(tempD)
                 if (min([s(3),s(2),s(1)]) ge 2) THEN BEGIN
                    CASE const_sub OF 
                       0: tempD = REFORM(tempD(i_s, i_l(0):i_r(0),i_l(1):i_r(1)))
                       1: tempD = REFORM(tempD(i_l(0):i_r(0),i_s,i_l(1):i_r(1)))
                       2: tempD = REFORM(tempD(i_l(0):i_r(0), i_l(1):i_r(1),i_s))
                       ELSE:
                    ENDCASE
                 ENDIF
              ENDIF
              
           END           
        ENDCASE
        
        Leftedgeh = DOUBLE([u_g(i).Left_edge(other_subs(0)),$
                     u_g(i).Left_edge(other_subs(1))])
;          print, 'size(temp):',size(temp), i_l, i_r
        IF (min([s(2),s(1)]) ge 2) THEN BEGIN
           FOR is=i_l(0),i_r(0) DO BEGIN 
              FOR js=i_l(1),i_r(1) DO BEGIN 
;print, 'running through grid'
                 x_stuetz_h = Leftedgeh(0) + $
                  DOUBLE(is)*delta_dist(other_subs(0))+ $
                  0.5D*delta_dist(other_subs(0))
                 y_stuetz_h = Leftedgeh(1) + $
                  DOUBLE(js)*delta_dist(other_subs(1))+ $
                  0.5D*delta_dist(other_subs(1))
; make sure there won't be a finer grid close to this location: 
; but don't worry for last grid:
                 con = 0 & & con_o = .1 & &  con_t=0.1
                 if (i lt num_of_grids-1) then begin
                    con_o = (le_o((i+1):(num_of_grids-1)) lt x_stuetz_h) and $
                     (le_t((i+1):(num_of_grids-1)) lt y_stuetz_h) 
                    con_t = (ri_o((i+1):(num_of_grids-1)) gt x_stuetz_h) and $
                     (ri_t((i+1):(num_of_grids-1)) gt y_stuetz_h)
;                          stop
                 endif
                 CASE 1 OF 
                    angular: BEGIN
                       radius =  DBLARR(3)
                       radius(const_sub) =  slice_ori(const_sub)
                       radius(other_subs) =  [x_stuetz_h, y_stuetz_h]
                       radius =  radius - center
                       jv =  CROSSP(radius, [tempx(is-i_l(0),js-i_l(1))-cenXvel,$
                                             tempy(is-i_l(0),js-i_l(1))-cenYvel,$
                                             tempz(is-i_l(0),js-i_l(1))-cenZvel])
                       jbar =  sqrt(total(jv*jv))
                       rbar =  sqrt(total(radius*radius))
                       jpr  =  jbar/rbar
                    END
                    pressure: $
                     jpr =  tempD(is-i_l(0),js-i_l(1))*tempT(is-i_l(0),js-i_l(1))
                    entropy:$
                     jpr =  tempD(is-i_l(0),js-i_l(1))^(-2./3.)*tempT(is-i_l(0),js-i_l(1))
                    ELSE: BEGIN 
                       jpr = temp(is-i_l(0),js-i_l(1))
                       IF  (divide_by_density(var_index) gt 0) THEN BEGIN
                          jpr = jpr/(tempD(is-i_l(0),js-i_l(1) >  1.e-20))/mass_per_particle(var_index)
                       ENDIF
                    END
                 ENDCASE
;                 if (((interpolate EQ 0) AND (i EQ 0)) OR $
                 if ( $
                      ((max(con_o and con_t) lt 1) AND (abs(jpr) GT 1.e-30))) THEN begin
                    x_stuetz(n_stuetz) = x_stuetz_h
                    y_stuetz(n_stuetz) = y_stuetz_h
                    z_stuetz(n_stuetz) = jpr
                    dx_stuetz(n_stuetz) =  delta_dist(0)
                    n_stuetz = n_stuetz+1
;                    IF (n_stuetz GT 16000) THEN GOTO, enough
                 endif
              ENDFOR            ; loop over single grid 
           ENDFOR               ; loop over single grid 
        END  ELSE print, 'warning grid too small !', i, i_r, i_l                   ; gird large enough?
     ENDIF
ENDFOR ; over grids ...

enough:
print, 'n_stuetz:', n_stuetz
x_stuetz = x_stuetz(0:n_stuetz-1)
y_stuetz = y_stuetz(0:n_stuetz-1)
z_stuetz = z_stuetz(0:n_stuetz-1)

SECMIN = MIN(z_stuetz(where(z_stuetz gt min(z_stuetz))))
missed =  1
IF logit THEN BEGIN 
   z_stuetz =  ALOG10(abs(z_stuetz) > abs(secmin))*z_stuetz/(abs(z_stuetz) > 1.e-30)
   print, 'construct_interpolated_slice_data: Took log of data value'
   missed = 0
ENDIF

secmin =  min(z_stuetz)

render_image, result, SECMIN=secmin, MISSING_VALUE=missed
slice_data = result

RETURN
END
;.compile construct_interpolated_slice_data


PRO render_image, result, SECMIN=secmin, MISSING_VALUE=missed
@common_blocks.inc
   interpolate =  interpolate_i
; interpolate all values
; get triangles:

   IF (NOT KEYWORD_SET(SECMIN)) THEN secmin = min(z_stuetz)
   
   IF verbose THEN  print, 'max before griding:',max(z_stuetz), min(z_stuetz)

;  print, 'b:',b
   result = FLTARR(xy_sl_size,xy_sl_size)
   IF N_ELEMENTS(missed) GT 0 THEN BEGIN
      IF missed THEN BEGIN 
         result(*) = 0. 
         IF verbose THEN print, 'render_image: set default value to 0.'
      END ELSE BEGIN
         result(*) = secmin
         IF verbose THEN print, 'render_image: set default value to ', secmin
      ENDELSE
   ENDIF

;   gs =  DOUBLE([(slice_size(2)-slice_size(0))/DOUBLE(xy_sl_size), $
;          (slice_size(3)-slice_size(1))/DOUBLE(xy_sl_size)])

   limits =  FLOAT([0., 0., 1., 1.])
   IF (interpolate NE 0) THEN BEGIN
      gs =  FLOAT([1./FLOAT(xy_sl_size), 1./FLOAT(xy_sl_size)])
      tx_stuetz =     FLOAT((x_stuetz-slice_size(0))/(slice_size(2)-slice_size(0))) 
      ty_stuetz =     FLOAT((y_stuetz-slice_size(1))/(slice_size(3)-slice_size(1))) 

      ind =  where((tx_stuetz GE 0. ) AND (tx_stuetz LE 1. ) AND (ty_stuetz GE 0. ) AND (ty_stuetz LE 1. ))
      tx_stuetz =    tx_stuetz(ind)
      ty_stuetz =    ty_stuetz(ind)
      tz_stuetz  =    z_stuetz(ind)
      TRIANGULATE, tx_stuetz, ty_stuetz, triangles, b
   ENDIF

   if (interpolate eq 2) then result = $
    TRIGRID(tx_stuetz,ty_stuetz, FLOAT(tz_stuetz), triangles, INPUT=result, $
            gs, limits, /quintic, min_value=secmin) <  max(z_stuetz)
   if (interpolate eq 1) then result = $
    TRIGRID(tx_stuetz,ty_stuetz, FLOAT(tz_stuetz), triangles, INPUT=result, $
            gs, limits, min_value=secmin)  <  max(z_stuetz)
   if (interpolate eq 0) then BEGIN 
; no interpolation:  
      dxts =  DOUBLE(xy_sl_size)/(DOUBLE(slice_size(2))-DOUBLE(slice_size(0)))
      dyts =  DOUBLE(xy_sl_size)/(DOUBLE(slice_size(3))-DOUBLE(slice_size(1)))
      FOR i=0L,n_stuetz-1L DO BEGIN
         xmi =  ROUND((x_stuetz(i)-0.5D*dx_stuetz(i)-DOUBLE(slice_size(0)))*dxts) >  0  < (xy_sl_size-1)
         xma =  ROUND((x_stuetz(i)+0.5D*dx_stuetz(i)-DOUBLE(slice_size(0)))*dxts) < (xy_sl_size-1) > 0
         ymi =  ROUND((y_stuetz(i)-0.5D*dx_stuetz(i)-DOUBLE(slice_size(1)))*dyts) > 0 < (xy_sl_size-1)
         yma =  ROUND((y_stuetz(i)+0.5D*dx_stuetz(i)-DOUBLE(slice_size(1)))*dyts) < (xy_sl_size-1) > 0
         IF ((ymi GE yma) or (xmi GE xma)) THEN result(xmi,ymi) = z_stuetz(i) $
         ELSE BEGIN
            hh =  DBLARR((xma-xmi+1),(yma-ymi+1))
            hh =  float(z_stuetz(i))
            result(xmi:xma,ymi:yma) = hh(*,*)
         ENDELSE
      ENDFOR
   ENDIF
   
   IF verbose THEN print, 'max after griding:',max(result), min(result)
   
   RETURN
END
;.compile construct_interpolated_slice_data
