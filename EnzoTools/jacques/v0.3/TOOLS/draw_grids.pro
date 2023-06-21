pro draw_grids, style, draw_cm, grid_info, slice_ori, slice_coord, xy_sl_size

if style eq 2  then  erase

; determine grids in slice
grid_flag = grid_in_slice(grid_info, slice_ori, slice_coord)
axes = (slice_ori gt 1.e-30)
const_sub = where(axes gt 1e-10)
other_subs = where(axes eq 0)
slice_size = [ xy_sl_size,  xy_sl_size]

min_left  = [slice_coord(0),slice_coord(1)]
max_right = [slice_coord(2),slice_coord(3)]
scale_up    = min(1./(max_right-min_left))

; use only grids that will show up on the image:
b = 0.9*(max_right-min_left)/FLOAT(slice_size)
; determine grids that are at least partially in the slice
grid_flag = grid_in_slice(grid_info, slice_ori, slice_coord)
index_range = grid_info.End_index-grid_info.Start_index+1
delta_distance = (grid_info.Right_edge - $
                  grid_info.Left_edge)/(index_range-1)

u_g = grid_info(FIX(where((delta_distance(other_subs,*) ge max(b)) $
                          and (grid_flag gt 0)  ) ))
index_range = u_g.End_index-u_g.Start_index+1
delta_distance = (u_g.Right_edge - $
                  u_g.Left_edge)/(index_range)

dim_grid_info = size(u_g)
num_of_grids  = dim_grid_info(1)

for i=0, num_of_grids-1 DO BEGIN
      delta_dist   = delta_distance(*,i)
      data_points  = index_range(*,i)  
      slice_ori_l                = fltarr(3)
      slice_ori_l(const_sub)     = slice_ori(const_sub)
      slice_ori_l(other_subs)    = [slice_coord(0),slice_coord(1)]
      slice_ori_r                = slice_ori_l 
      slice_ori_r(other_subs)    = [slice_coord(2),slice_coord(3)]

      i_s = ROUND((slice_ori_l(const_sub) - $
                   u_g(i).Left_edge(const_sub))/delta_dist(const_sub) >0 $
        < data_points(const_sub)-1)
      i_l = CEIL( (slice_ori_l(other_subs) - $
                   u_g(i).Left_edge(other_subs)) $
                    /delta_dist(other_subs) ) > 0
      i_r   = (FLOOR((slice_ori_r(other_subs) - $
                  u_g(i).Left_edge(other_subs))/ $
                    delta_dist(other_subs) - 1) < $
        (data_points(other_subs)-1) ) > 0
;      print, 'slize_ori:',  slice_ori_l, slice_ori_r
;      print, 'index', i,  i_s, i_r, i_l
;      temp_size  = data_points(other_subs)
;          data_slice = (extract_slice(temp,temp_size(0),temp_size(1), $
;                index(0),index(1),index(2),0.,0.,0.))
      L_p = FIX(FLOAT(slice_size) *scale_up * $
                         ((u_g(i).Left_edge(other_subs)+ $
                           (i_l)*delta_dist(other_subs)-min_left) > 0.))
      R_p= L_p + $
        Round(FLOAT(slice_size) *scale_up * $
              float(i_r-i_l+1)*delta_dist(other_subs)) < (slice_size-1)
  
    box = [[L_p(0), L_p(1)], $
           [L_p(0), R_p(1)], $
           [R_p(0), R_p(1)], $
           [R_p(0), L_p(1)], $
           [L_p(0), L_p(1)]]

    CASE style OF 
        1 : BEGIN
            plots, box, /device, LINESTYLE = 1
            xyouts, R_p(0), R_p(1), STRCOMPRESS(STRING(u_g(i).level)), $
              CHARSIZE = 0.7, /DEVICE, ALIGNMENT = 1.0
        END
        ELSE: BEGIN 
            colors = indgen(16)*10+20
            polyfill, box+[-1.,-1.,-1.,1.,1.,1.,1.,-1.,-1.,-1.],$
              /device, col = colors(u_g(i).level)
            next_level = 0
            if i lt num_of_grids-1 then next_level = u_g(i+1).level $
            else next_level =200
            current_level = u_g(i).level 
            if next_level gt current_level then BEGIN
;    xyouts, R_p(0), R_p(1), STRCOMPRESS(STRING(u_g(i).num)), $
                xyouts, R_p(0), R_p(1), STRCOMPRESS(STRING(u_g(i).level)), $
                  CHARSIZE = 0.7, /DEVICE, ALIGNMENT = 1.0
            ENDIF
        END
    END
ENDFOR

; for the drawe chart erase colormap window 
IF style ne 1 THEN BEGIN
    old_wind = !d.window
    widget_control, get_value = win_index, draw_cm
    wset, win_index
    erase
    wset, old_wind
ENDIF
END

