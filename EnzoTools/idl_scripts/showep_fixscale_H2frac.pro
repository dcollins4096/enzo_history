pro showep_fixscale_H2frac, file_name, outfile_name

; user sets these
colortable = 4
fixed_upper = -2.0
fixed_lower = -5.0
fixed_divisions = 3



file_id = H5F_OPEN(file_name)
dataset_id=H5D_OPEN(file_id,'/H2I_Density')
HI_data = H5D_READ(dataset_id)
dataspace_id = H5D_GET_SPACE(dataset_id)

dimensions = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)

print, 'dimensions: ', dimensions

H5S_CLOSE, dataspace_id
H5D_CLOSE, dataset_id
H5F_CLOSE, file_id


file_id = H5F_OPEN(file_name)
dataset_id=H5D_OPEN(file_id,'/Density')
density_data = H5D_READ(dataset_id)
dataspace_id = H5D_GET_SPACE(dataset_id)

dimensions = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)

print, 'dimensions: ', dimensions

H5S_CLOSE, dataspace_id
H5D_CLOSE, dataset_id
H5F_CLOSE, file_id


loadct,colortable


logdata = alog10(HI_data/(0.76*density_data) + 0.000000000001)

;logdata = alog10(data * 8.349e-17 )

print, 'minimum of log data is:  ', min(logdata)
print, 'maximum of log data is:  ',max(logdata)

;bigdimx = dimensions[0]+150
;bigdimy = dimensions[1]+150

bigdimx = 512
bigdimy = 512

print, 'big dimensions', bigdimx, bigdimy

;window, xsize=bigdimx, ysize=bigdimy
;erase

set_plot,'z'
device,set_resolution=[bigdimx,bigdimy]
erase
scaled_data = bytscl(logdata, min=fixed_lower, max=fixed_upper, top=255)


tv,scaled_data,0,0

image=TVRD()

print, size(image)

image3D = bytarr(3,bigdimx,bigdimy) 
tvlct,r,g,b,/get

image3d[0,*,*] = r[image]
image3d[1,*,*] = g[image]
image3d[2,*,*] = b[image]

write_jpeg,outfile_name,image3d,true=1,quality=100


end
