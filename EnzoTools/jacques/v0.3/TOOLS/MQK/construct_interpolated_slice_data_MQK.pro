;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/17/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + read data by dataset_name, not sds_num

pro construct_interpolated_slice_data, $
              slice_data, SECMIN=SECMIN
; PLOT DERIVED QUANTITIES such as angular momentum, Entropy and Pressure ...
@TOOLS/common_blocks.inc
; using grid_info read in data (found in data_dir directory) of all
; necessary grids for a slice ; specified by the 3 component vector
; slice_ori. slice_size has 4 components specifying the left and right
; edge of the slice
; note that slice_ori has to have two elements set to zero since we
; only do orthogonal slices.
; 

;; set dataset_name
  dataset_name=list_str[var_index]
  print,'dataset_name='+dataset_name

;; determine some flags
  angular=0
  entropy=0
  pressure=0
  case 1 of
      (strpos(list_str[var_index], 'Angul') ge 0): angular=1
      (strpos(list_str[var_index], 'Entropy') ge 0): entropy=1
      (strpos(list_str[var_index], 'Pressure') ge 0): pressure=1
      else:
  endcase
  
  interpolate =  interpolate_i
  
; should I log it ?
  logit = 1
  if(strpos(dataset_name,'-velo') ge 0) then logit=0 
  
  if angular then $
    if verbose then print, 'con_angular_momentum...: Using center:', center

  num_dim = 3
; rough check of input data
  slice_data = dblarr(xy_sl_size, xy_sl_size)
  axes = (slice_ori gt 1.e-30)
  const_sub =  where(axes eq max(axes))
  const_sub = const_sub[0]
  other_subs = where(axes eq 0)
  
  if(fix(total(axes)) gt 1) then $
    message,'slice_ori = '+string(slice_ori,format='(3(i2,2x))')+'not allowed!'

  if(n_elements(grid_info) lt 1) then $
      message,'No grid information in grid_info:'

  dim_grid_info=size(grid_info)
  num_of_grids=dim_grid_info[1]


;; diagnostics
  print,'dim_grid_info: ',dim_grid_info
  print,'num_of_grids: ',num_of_grids
  print,' '
  print,'dim: ',grid_info[0].dim
  print,'start_index: ',grid_info[0].start_index
  print,'end_index: ',grid_info[0].end_index
  print,'left_edge: ',grid_info[0].left_edge
  print,'right_edge: ',grid_info[0].right_edge


; min left edge and scaling
  min_left  = double([slice_size[0],slice_size[1]])
  max_right = double([slice_size[2],slice_size[3]]) 
;  print, 'min_left: ',min_left
;  print, 'max_right: ',max_right
  scale_up    = min(1.D/(max_right-min_left))
  
  index_range    = grid_info.End_index-grid_info.Start_index+1
  delta_distance = double(grid_info.Right_edge - $
                          grid_info.Left_edge)/double(index_range)
  
; find velocity of center:
  if angular then begin 
      ml =  -1
      this_grid =  1
      for i =0,num_OF_grids-1 do begin
          if  (((grid_info[i].left_edge[0] - center[0]) gt 0.) and $
               ((grid_info[i].left_edge[1] - center[1]) gt 0.) and $
               ((grid_info[i].left_edge[2] - center[2]) gt 0.) and $
               ((-grid_info[i].right_edge[0]+center[0]) lt 0.) and $
               ((-grid_info[i].right_edge[1]+center[1]) lt 0.) and $
               ((-grid_info[i].right_edge[2]+center[2]) lt 0.)) then begin 
              if (grid_info[i].level gt ml) then begin 
                  this_grid = i 
                  ml =  grid_info[i].level 
              endif
          endif
      endfor
      print, 'this_grid:', this_grid, grid_info[0:5].level

; read that grid
      read_grid,grid_info[this_grid].baryon_file,tempx,$
        byname='x-velocity'
      read_grid,grid_info[this_grid].baryon_file,tempy,$
        byname='y-velocity'
      read_grid,grid_info[this_grid].baryon_file,tempz,$
        byname='z-velocity'

      ix =  round((center[0]-grid_info[this_grid].left_edge[0])*delta_distance[this_grid])
      iy =  round((center[1]-grid_info[this_grid].left_edge[1])*delta_distance[this_grid])
      iz =  round((center[2]-grid_info[this_grid].left_edge[2])*delta_distance[this_grid])

      cenXvel =  tempx[ix,iy,iz]
      cenYvel =  tempy[ix,iy,iz]
      cenZvel =  tempz[ix,iy,iz]

      if verbose then print, 'central velocity:', cenXvel, cenYvel, cenZvel
  endif                         ; if angular 
  
; use only grids that will show up on the image:
  dxts = double([xy_sl_size,xy_sl_size])/(max_right-min_left)
; determine grids that are at least partially in the slice
  grid_flag      = grid_in_slice(grid_info, double(slice_ori), double(slice_size))
  index_range    = grid_info.End_index-grid_info.Start_index+1
  delta_distance = double(grid_info.Right_edge - $
                          grid_info.Left_edge)/double(index_range)
  del_pix =   double(max([xy_sl_size,xy_sl_size]))*delta_distance/(max_right[0]-min_left[0])*double(index_range)
; initialize minimum of data 
  min_data = 1.e30
  u_g = grid_info[where((del_pix[other_subs[0],*] ge .5) and $
                        (grid_flag gt 0))]
  
; print,  del_pix
  dim_grid_info = size(u_g)
  num_of_grids  = dim_grid_info[1]
  
  fully_contained_above =  0
  for i =0, num_of_grids-1 do begin
      if ((u_g[i].Left_edge[other_subs[0]] le (min_left[0]))$
          and  (u_g[i].Left_edge[other_subs[1]] le (min_left[1])) $
          and (u_g[i].Right_edge[other_subs[0]] ge (max_right[0])) $
          and (u_g[i].Right_edge[other_subs[1]] ge (max_right[1]))) Then $
        fully_contained_above = i
  endfor
  
  print, 'slice is fully contained in grids after index:',  $
    fully_contained_above
  
;; diagnostics
;  print,'dim: ',u_g[0].dim
;  print,'start_index: ',u_g[0].start_index
;  print,'end_index: ',u_g[0].end_index
;  print,'left_edge: ',u_g[0].left_edge
;  print,'right_edge: ',u_g[0].right_edge


  index_range = u_g.End_index-u_g.Start_index+1
  delta_distance = (u_g.Right_edge - $
                    u_g.Left_edge)/double(index_range)

;  print, size(u_g)
;  print,  'max:', min_left, max_right, scale_up
  x_stuetz =  dblarr(1000000l)
  y_stuetz =  dblarr(1000000l)
  z_stuetz =  dblarr(1000000l)
  dx_stuetz = dblarr(1000000l)
  
  n_stuetz = 0L
  this_level=0
  le_o    = u_g[0:num_of_grids-1].left_edge[other_subs[0]]
  le_t    = u_g[0:num_of_grids-1].left_edge[other_subs[1]]
  ri_o    = u_g[0:num_of_grids-1].right_edge[other_subs[0]]
  ri_t    = u_g[0:num_of_grids-1].right_edge[other_subs[1]]
  
  slice_ori_l                = dblarr(3)
  slice_ori_l[const_sub]     = slice_ori[const_sub]
  slice_ori_l[other_subs]    = min_left
  slice_ori_r                = slice_ori_l 
  slice_ori_r[other_subs]    = max_right
; 
  Num_pix_x = xy_sl_size
  Num_pix_y = xy_sl_size
  
  i=0l
  for i=fully_contained_above, num_of_grids-1 do begin

; check this
      delta_dist   = delta_distance[*,i]
      data_points  = index_range[*,i]  

;; diagnostics
;      print,'slice_ori_l: ',slice_ori_l
;      print,'slice_ori_r: ',slice_ori_r
;      print,'delta_dist: ',delta_dist
;      print,'data_points: ',data_points

; yep check again !!
      i_s = round((slice_ori_l[const_sub] - 0.5*delta_dist[const_sub] -  $
                   u_g[i].Left_edge[const_sub])/delta_dist[const_sub] >0 $
                  < (data_points[const_sub]-1))
      i_l =   (floor((slice_ori_l[other_subs] - $
                      u_g[i].Left_edge[other_subs]) $
                     /delta_dist[other_subs]) -1) > 0
      i_r   = ((ceil((slice_ori_r[other_subs] - $
                      u_g[i].Left_edge[other_subs])/ $
                     delta_dist[other_subs]) ) < $
               (data_points[other_subs]-1))  > 0

;      print,'i_s=',i_s
;      print,'i_l=',i_l
;      print,'i_r=',i_r

      if ((min(i_r-i_l) ge 1)  ) then begin
          case 1 of
              angular: begin
                  read_grid, u_g[i].baryon_file,tempx,$
                    byname='x-velocity'
                  read_grid, u_g[i].baryon_file,tempy,$
                    byname='y-velocity'
                  read_grid, u_g[i].baryon_file,tempz,$
                    byname='z-velocity'

                  s = size(tempx)
                  if (min([s[3],s[2],s[1]]) ge 2) then begin
                      case const_sub of 
                          0: tempx = reform(tempx[i_s, i_l[0]:i_r[0],i_l[1]:i_r[1]])
                          1: tempx = reform(tempx[i_l[0]:i_r[0],i_s,i_l[1]:i_r[1]])
                          2: tempx = reform(tempx[i_l[0]:i_r[0], i_l[1]:i_r[1],i_s])
                          else:
                      endcase
                      case const_sub of 
                          0: tempy = reform(tempy[i_s, i_l[0]:i_r[0],i_l[1]:i_r[1]])
                          1: tempy = reform(tempy[i_l[0]:i_r[0],i_s,i_l[1]:i_r[1]])
                          2: tempy = reform(tempy[i_l[0]:i_r[0], i_l[1]:i_r[1],i_s])
                          else:
                      endcase
                      case const_sub of 
                          0: tempz = reform(tempz[i_s, i_l[0]:i_r[0],i_l[1]:i_r[1]])
                          1: tempz = reform(tempz[i_l[0]:i_r[0],i_s,i_l[1]:i_r[1]])
                          2: tempz = reform(tempz[i_l[0]:i_r[0], i_l[1]:i_r[1],i_s])
                          else:
                      endcase
                  end
              end
              (entropy or pressure): begin 
                  read_grid, u_g[i].baryon_file,tempD,$
                    byname='Density'
                  read_grid, u_g[i].baryon_file,tempT,$
                    byname='Temperature'

                  s = size(tempD)
                  if (min([s[3],s[2],s[1]]) ge 2) then begin
                      case const_sub of 
                          0: tempD = reform(tempD[i_s, i_l[0]:i_r[0],i_l[1]:i_r[1]])
                          1: tempD = reform(tempD[i_l[0]:i_r[0],i_s,i_l[1]:i_r[1]])
                          2: tempD = reform(tempD[i_l[0]:i_r[0], i_l[1]:i_r[1],i_s])
                          else:
                      endcase
                      case const_sub of 
                          0: tempT = reform(tempT[i_s, i_l[0]:i_r[0],i_l[1]:i_r[1]])
                          1: tempT = reform(tempT[i_l[0]:i_r[0],i_s,i_l[1]:i_r[1]])
                          2: tempT = reform(tempT[i_l[0]:i_r[0], i_l[1]:i_r[1],i_s])
                          else:
                      endcase
                  endif
              end               ; entropy or pressure
              else: begin       ; 'normal' quantity
;                  print,'reading '+dataset_name
                  read_grid, u_g[i].baryon_file,temp,$
                    byname=dataset_name

;                  help,temp

                  s = size(temp)
                  if (min([s[3],s[2],s[1]]) ge 2) then begin
                      case const_sub of 
                          0: temp = reform(temp[i_s, i_l[0]:i_r[0],i_l[1]:i_r[1]])
                          1: temp = reform(temp[i_l[0]:i_r[0],i_s,i_l[1]:i_r[1]])
                          2: temp = reform(temp[i_l[0]:i_r[0], i_l[1]:i_r[1],i_s])
                          else: 
                      endcase
                  end else begin 
                      print, 'warning: very small grid '
                  endelse
                  if(divide_by_density[var_index] gt 0) then begin
                      read_grid, u_g[i].baryon_file,tempD,$
                        byname='Density'

                      s = size(tempD)
                      if (min([s[3],s[2],s[1]]) ge 2) then begin
                          case const_sub of 
                              0: tempD = reform(tempD[i_s, i_l[0]:i_r[0],i_l[1]:i_r[1]])
                              1: tempD = reform(tempD[i_l[0]:i_r[0],i_s,i_l[1]:i_r[1]])
                              2: tempD = reform(tempD[i_l[0]:i_r[0], i_l[1]:i_r[1],i_s])
                              else:
                          endcase
                      endif
                  endif
                  
              end           
          endcase
          
          Leftedgeh = double([u_g[i].Left_edge(other_subs[0]),$
                              u_g[i].Left_edge(other_subs[1])])
          if (min([s[2],s[1]]) ge 2) then begin
              for is=i_l[0],i_r[0] do begin 
                  for js=i_l[1],i_r[1] do begin 
;print, 'running through grid'
                      x_stuetz_h = Leftedgeh[0] + $
                        double(is)*delta_dist[other_subs[0]]+ $
                        0.5d*delta_dist[other_subs[0]]
                      y_stuetz_h = Leftedgeh[1] + $
                        double(js)*delta_dist[other_subs[1]]+ $
                        0.5d*delta_dist[other_subs[1]]
; make sure there won't be a finer grid close to this location: 
; but don't worry for last grid:

; mqk: this caused problems - remove for now
                      con_o=0
                      con_t=0

;                      con = 0 & & con_o = .1 & &  con_t=0.1
;                      if (i lt num_of_grids-1) then begin
;                          con_o = (le_o[(i+1):(num_of_grids-1)] lt x_stuetz_h) and $
;                            (le_t[(i+1):(num_of_grids-1)] lt y_stuetz_h) 
;                          con_t = (ri_o[(i+1):(num_of_grids-1)] gt x_stuetz_h) and $
;                            (ri_t[(i+1):(num_of_grids-1)] gt y_stuetz_h)
;;                          stop
;                      endif
                      case 1 of 
                          angular: begin
                              radius =  dblarr(3)
                              radius[const_sub] =  slice_ori[const_sub]
                              radius[other_subs] =  [x_stuetz_h, y_stuetz_h]
                              radius =  radius - center
                              jv =  crossp(radius, [tempx[is-i_l[0],js-i_l[1]]-cenXvel,$
                                                    tempy[is-i_l[0],js-i_l[1]]-cenYvel,$
                                                    tempz[is-i_l[0],js-i_l[1]]-cenZvel])
                              jbar =  sqrt(total(jv*jv))
                              rbar =  sqrt(total(radius*radius))
                              jpr  =  jbar/rbar
                          end
                          pressure: $
                            jpr =  tempD[is-i_l[0],js-i_l[1]]*tempT[is-i_l[0],js-i_l[1]]
                          entropy:$
                            jpr =  tempD[is-i_l[0],js-i_l[1]]^(-2./3.)*tempT[is-i_l[0],js-i_l[1]]
                          else: begin 
                              jpr = temp[is-i_l[0],js-i_l[1]]
                              if(divide_by_density[var_index] gt 0) then begin
                                  jpr = jpr/(tempD[is-i_l[0],js-i_l[1] >  1.e-20])/mass_per_particle[var_index]
                              endif
                          end
                      endcase
;                 if (((interpolate EQ 0) AND (i EQ 0)) OR $
                      if ( $
                           ((max(con_o and con_t) lt 1) AND (abs(jpr) GT 1.e-30))) then begin
                          x_stuetz[n_stuetz] = x_stuetz_h
                          y_stuetz[n_stuetz] = y_stuetz_h
                          z_stuetz[n_stuetz] = jpr
                          dx_stuetz[n_stuetz] =  delta_dist[0]
                          n_stuetz = n_stuetz+1
;                    IF (n_stuetz GT 16000) THEN GOTO, enough
                      endif
                  endfor        ; loop over single grid 
              endfor            ; loop over single grid 
          end else print, 'warning grid too small !', i, i_r, i_l ; grid large enough?
      endif
  endfor                        ; over grids ...
  
enough:
  print, 'n_stuetz:', n_stuetz
  x_stuetz = x_stuetz[0:n_stuetz-1]
  y_stuetz = y_stuetz[0:n_stuetz-1]
  z_stuetz = z_stuetz[0:n_stuetz-1]
  
  secmin = min(z_stuetz[where(z_stuetz gt min(z_stuetz))])
  missed =  1
  if logit then begin 
      z_stuetz =  alog10(abs(z_stuetz) > abs(secmin))*z_stuetz/(abs(z_stuetz) > 1.e-30)
      print, 'construct_interpolated_slice_data: Took log of data value'
      missed = 0
  endif
  
  secmin =  min(z_stuetz)

  
  render_image, result, secmin=secmin, missing_value=missed
  slice_data = result

  return
end
;.compile construct_interpolated_slice_data

