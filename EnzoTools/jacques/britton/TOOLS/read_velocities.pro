pro read_velocities
@TOOLS/common_blocks.inc
; find velocity fields         
   sds_in =  sds_num
  found_at = 0 
  FOR i=0,N_ELEMENTS(list_str)-1 DO BEGIN
      check = STRPOS(list_str(i),'elocity') 
      if ((check gt -1) and (found_at eq 0)) THEN found_at = i
  ENDFOR
; assume it is 3D and the following are y and z
  os = where(slice_ori lt 1.e-15)
  IF (N_ELEMENTS(os) eq 2) THEN BEGIN
     sds_num =  found_at + os(0) + 1
     var_index =  sds_num - 1
     construct_interpolated_slice_data, vel_one
     sds_num =  found_at + os(1) + 1
     var_index =  sds_num - 1
     construct_interpolated_slice_data, vel_two
  END         
  sds_num =  sds_in
  var_index =  sds_num - 1
  RETURN
END 
;.compile read_velocities.pro

