pro construct_image
; constructs slice and displays it
; returns the image structure
@TOOLS/common_blocks.inc

help,var_index
help,divide_by_density

if  (verbose and (divide_by_density(var_index) gt 0)) THEN BEGIN
   print, 'also get density data to get ratio'
   print, 'particle mass:', mass_per_particle(var_index)
ENDIF

construct_interpolated_slice_data, data, SECMIN=sec_min 
print, '1st min', sec_min

min_data = min(data, MAX=max_data)
if verbose then print, '2nd min data,max data:',min_data, max_data

image_d = DBLARR(1000, 1000)
image_d(0:xy_sl_size-1,0:xy_sl_size-1) = data
image_stack(0).slice_size = slice_size
;shift image_stack
for j=N_images-1,1,-1 DO image_stack(j)=image_stack(j-1)

image_stack(0) = {image, list_str(var_index), grid_info(0).hier_file, $
         image_d, xy_sl_size, min_data, max_data, DOUBLE(slice_size)}

END

;.compile construct_image
