;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/17/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + call construct_3D_data with dataset name (list_str[sds]), not sds
;;   number.
;;

pro construct_cube, data, sds, MIN=min_data, MAX=max_data
; constructs cube 
; this is the wrapper around construct_3D_data used to compute
; fractions and do other special data dependent actions 
; returns the data cube
@TOOLS/common_blocks.inc


  cube_center = slice_ori
  cube_center[where(slice_ori  lt 1.e-30)] = $
    [slice_size[0]+0.5*(slice_size[2]-slice_size[0]), $
     slice_size[1]+0.5*(slice_size[3]-slice_size[1])]
  cube_length = slice_size[2]-slice_size[0]
  
;;  print,'cube_center: ',cube_center
;;  print,'cube_length: ',cube_length

  construct_3D_data, grid_info, cube_center, cube_length, $
    data_dir, list_str[sds], cube_dim, data, $
    smooth = threeD_smoothing

; fudge problems with density and temperature data
  if (sds eq 1) or (sds eq 16) then data = data > 1.e-10
  data_min =  min(data)
  if data_min ge 0. then begin 
      sec_min_data = min(data(where(data ne data_min)))
      data = alog10(data > sec_min_data)
  end
  if  (divide_by_density(sds-1) gt 0) then begin
      if verbose then begin 
          print, 'also get density data to get ratio'
          print, 'particle mass:', mass_per_particle(var_index)
      end
      construct_3D_data,   grid_info, cube_center, cube_length, $
        data_dir, list_str[1], cube_dim, data_den, smooth = threeD_smoothing
      sec_min_data = min(data_den(where(data_den ne min(data_den))))
      data_den = alog10(data_den > sec_min_data)
      if verbose then $
        print, 'divide by particle mass:',mass_per_particle(sds-1)
      data = data - data_den - alog10(mass_per_particle(sds-1))
  endif
  
;if (STRPOS(list_str(var_index),'Dark') gt -1) then $
;  data = alog10(data>1.e-1)
  
; return min and max if requested 
  min_data = min(data, max=max_data)
  
  return
end
; .compile construct_cube

