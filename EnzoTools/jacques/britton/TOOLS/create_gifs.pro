PRO create_gifs, filename, data_dir,sds_num, xy_sl_size, $
          mass_per_particle,divide_by_density, $
          slice_size, grid_info, data_in_dir, gif_dir, xyz_max, minimum, $
          maximum
; create gifs of slices centered around the highest density peak  
FOR i=0,0 DO BEGIN
    slice_ori    = [0.,0.,0.]
    slice_ori(i) = xyz_max(i)

    other_subs   = where(slice_ori eq 0.)
    print, other_subs
    slice_coord    = DBLARR(4)
    slice_coord(0:1)  = xyz_max(other_subs)-0.5*slice_size
    slice_coord(2:3)  =  slice_coord(0:1)+slice_size

    construct_slice_data, grid_info, slice_ori, slice_coord, data_dir, $
      sds_num, [xy_sl_size,xy_sl_size], data
    data_min =  min(data)
    if data_min ge 0. THEN BEGIN 
        sec_min_data = min(data(where(data ne data_min)))
        data = alog10(data > sec_min_data)
    END
    IF  (divide_by_density(sds_num) gt 0) THEN BEGIN
        print, 'also get density data to get ratio'
        print, 'particle mass:', mass_per_particle(var_index)
        construct_slice_data, grid_info, slice_ori, slice_coord, data_dir, $
          1, [xy_sl_size,xy_sl_size], data_den, SMOOTH = smoothing
        sec_min_data = min(data_den(where(data_den ne min(data_den))))
        data_den = alog10(data_den > sec_min_data)
        print, 'divide by particle mass:',mass_per_particle(var_index)
        data = data - data_den - alog10(mass_per_particle(var_index))
    ENDIF
    if i eq 0 then file_name=filename+'x'
    if i eq 1 then file_name=filename+'y'
    if i eq 2 then file_name=filename+'z'

    TV, BYTSCL(data, MIN=minimum,MAX=maximum)
    print, 'write gif to ', file_name
    WRITE_GIF, file_name, TVRD()

END

END

PRO animate
; hard wired animation ...
data_in_dir = '/local/tabel/First/f128L12/'
file_base   = '128_halo1_'
gif_dir     = '/local/tabel/First/GIFs/f128L12/'

num = [40, 41, 42, 43, 44, 45, 46, 47, 49, 51, 52, 53, 54]
tot = N_elements(num)

xy_sl_size = 500

window,0, xsize=xy_sl_size, ysize=xy_sl_size
image = BYTARR(xy_sl_size,xy_sl_size,tot)
xinteranimate, set=[xy_sl_size,xy_sl_size,tot],/SHOWLOAD
FOR j=0, tot-1 DO BEGIN 
    this = num(j)
; DENSITY 
    proj = 'x'
    gif_base = 'BD128.'
    filename = gif_dir+gif_base + $
      STRCOMPRESS(STRING(this, FORMAT='(I4.4)')) +'.gif.'+proj
    read_gif, filename, image1,R,G,B
;    TVLCT, R,G,B
;    TV, image1
    image(*,*,j) = image1
ENDFOR

FOR j=0, tot-1 DO xinteranimate, FRAME=j, image=image(*,*,j)
xinteranimate, /keep_pixmaps


END


;animate
;STOP



; MAIN
;@TOOLS/tools.com
;data_in_dir = '/local/tabel/First/first128/first128_halo1/'
data_in_dir = '/local/tabel/First/first128/'
;file_base   = 'RedshiftOutput'
file_base   = '128_halo1_'

gif_dir     = '/local/tabel/First/GIFs/r128L12/'

num = [51, 52, 53, 54, 55, 56, 57, 59, 60, 61]
tot = N_elements(num)

xy_sl_size = 500
;slice_size = 0.05
slice_size = 0.001

window,0, xsize=xy_sl_size, ysize=xy_sl_size

FOR j=0, tot-1 DO BEGIN 
    this = num(j)
;    data_this_dir = data_in_dir+ 'time'+ $
;      STRCOMPRESS(STRING(this, FORMAT='(I2.2)'))+'/'
    data_this_dir = data_in_dir
    file = data_this_dir +  $
      file_base + STRCOMPRESS(STRING(this, FORMAT='(I4.4)')) + '.hierarchy'
    READ_GRID_INFO, file, Grid_info, $
        list_str, mass_per_particle, divide_by_density   

    
   find_max_value, 'FIND MAX.from the finest levels',$
      grid_info, data_this_dir, 1, xyz_max

; DENSITY 
loadct, 15
    minimum = 2.
    maximum = 9.

    sds_num = 1
    gif_base = 'BD128.'
    filename = gif_dir+gif_base+STRCOMPRESS(STRING(this, FORMAT='(I4.4)')) +'.gif.'
    CREATE_GIFS, filename, data_this_dir, sds_num, xy_sl_size, $
      mass_per_particle,divide_by_density, $
      slice_size, grid_info, data_this_dir, gif_dir, xyz_max, minimum, $
      maximum

;; DM DENSITY 
;    minimum = -1.
;    maximum = 8.3

;    sds_num = 16+1
;    gif_base = 'DM128.'
;    filename = gif_dir+gif_base+STRCOMPRESS(STRING(this, FORMAT='(I4.4)')) +'.gif.'
;    CREATE_GIFS, filename, data_this_dir, sds_num, xy_sl_size, $
;      mass_per_particle,divide_by_density, $
;      slice_size, grid_info, data_this_dir, gif_dir, xyz_max, minimum, $
;      maximum

;; TEMPERATURE
;loadct, 3
;    minimum = 1.6
;    maximum = 3.

;    sds_num = 16
;    gif_base = 'T128.'
;    filename = gif_dir+gif_base+STRCOMPRESS(STRING(this, FORMAT='(I4.4)')) +'.gif.'
;    CREATE_GIFS, filename, data_this_dir, sds_num, xy_sl_size, $
;      mass_per_particle,divide_by_density, $
;      slice_size, grid_info, data_this_dir, gif_dir, xyz_max, minimum, $
;      maximum

ENDFOR

END
;.run TOOLS/create_gifs.pro
