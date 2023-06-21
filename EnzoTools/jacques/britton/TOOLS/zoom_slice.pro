pro zoom_slice, grid_info, slice_ori, data_dir, sds_num, slice_size

window, xsize = FIX(slice_size(0)),  ysize = FIX(slice_size(1))
!err = 0
wind = [0.,0.,1.,1.]
while (!ERR ne 4) DO begin  ; repeat until right button is pressed
    construct_slice_data,grid_info,slice_ori, wind, data_dir, $
      sds_num, slice_size, slice_data
    TV, bytscl(alog10(slice_data>.1))
     read_box_from_window, box
;normalize pixel positions
;     windx = wind(2)-wind(0)
     windy = max([wind(3)-wind(1),wind(2)-wind(0)])
     print, 'windy:', windy
  wind =  [wind(0) + windy*box(0,0)/slice_size(0), $
           wind(1) + windy*box(0,1)/slice_size(1), $
           wind(0) + windy*box(1,0)/slice_size(0), $
           wind(1) + windy*box(1,1)/slice_size(1)]
  if (!err eq 2) then wind=[0.,0.,1.,1.]
  print, 'wind', wind
ENDWHILE
END

