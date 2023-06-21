pro construct_cube, data, sds, MIN=min_data, MAX=max_data
; constructs cube 
; this is the wrapper around construct_3D_data used to compute
; fractions and do other special data dependent actions 
; returns the data cube
@TOOLS/common_blocks.inc


cube_center = slice_ori
cube_center(where(slice_ori  lt 1.e-30)) = $
  [slice_size(0)+0.5*(slice_size(2)-slice_size(0)), $
   slice_size(1)+0.5*(slice_size(3)-slice_size(1))]
cube_length = slice_size(2)-slice_size(0)

construct_3D_data,   grid_info, cube_center, cube_length, $
  data_dir, sds, cube_dim, data, SMOOTH = threeD_smoothing
; fudge problems with density and temperature data
if (sds eq 1) or (sds eq 16) then data = data > 1.e-10
data_min =  min(data)
if data_min ge 0. THEN BEGIN 
    sec_min_data = min(data(where(data ne data_min)))
    data = alog10(data > sec_min_data)
END
IF  (divide_by_density(sds-1) gt 0) THEN BEGIN
    if verbose then BEGIN 
        print, 'also get density data to get ratio'
        print, 'particle mass:', mass_per_particle(var_index)
    END
    construct_3D_data,   grid_info, cube_center, cube_length, $
           data_dir, 1, cube_dim, data_den, SMOOTH = threeD_smoothing
    sec_min_data = min(data_den(where(data_den ne min(data_den))))
    data_den = alog10(data_den > sec_min_data)
    if verbose then $
      print, 'divide by particle mass:',mass_per_particle(sds-1)
    data = data - data_den - alog10(mass_per_particle(sds-1))
ENDIF

;if (STRPOS(list_str(var_index),'Dark') gt -1) then $
;  data = alog10(data>1.e-1)

; return min and max if requested 
min_data = min(data, max=max_data)

RETURN
END
; .compile construct_cube

