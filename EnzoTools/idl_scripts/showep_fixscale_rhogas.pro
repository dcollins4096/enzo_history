pro showep_fixscale_rhogas, file_name, outfile_name,caption,caption2

; user sets these
colortable = 4
fixed_upper = 15.0
fixed_lower = 11.0
fixed_divisions = 4



file_id = H5F_OPEN(file_name)
dataset_id=H5D_OPEN(file_id,'/Density')
data = H5D_READ(dataset_id)
dataspace_id = H5D_GET_SPACE(dataset_id)

dimensions = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)

print, 'dimensions: ', dimensions

H5S_CLOSE, dataspace_id
H5D_CLOSE, dataset_id
H5F_CLOSE, file_id

loadct,colortable


logdata = alog10(data + 0.0001)

;logdata = alog10(data * 8.349e-17 )

print, 'minimum of log data is:  ', min(logdata)
print, 'maximum of log data is:  ',max(logdata)

;bigdimx = dimensions[0]+150
;bigdimy = dimensions[1]+150

bigdimx = 768
bigdimy = 768

print, 'big dimensions', bigdimx, bigdimy

;window, xsize=bigdimx, ysize=bigdimy
;erase

set_plot,'z'
device,set_resolution=[bigdimx,bigdimy]
erase
scaled_data = bytscl(logdata, min=fixed_lower, max=fixed_upper, top=255)

print, 'min, max of scaled data is: ',min(scaled_data),max(scaled_data)
print, 'fixed_lower, fixed_upper are: ',fixed_lower, fixed_upper
print, ''
tv,scaled_data,180,128

cbarrange = [fixed_lower,fixed_upper]

cbarposition =  [0.1, 0.1, 0.85, 0.15]  

;colorbar=replicate(1b,20) # bindgen(256)
colorbar,ncolors=256,position=cbarposition,RANGE = cbarrange,vertical=1,divisions=fixed_divisions

xyouts, 300,670,'Log Density',/device,charsize=1.5, color=255
xyouts, 190,95,caption,/device,charsize=1.5,color=255
xyouts, 450,95,caption2,/device,charsize=1.5,color=255


image=TVRD()

print, size(image)

image3D = bytarr(3,bigdimx,bigdimy) 
tvlct,r,g,b,/get

image3d[0,*,*] = r[image]
image3d[1,*,*] = g[image]
image3d[2,*,*] = b[image]

write_jpeg,outfile_name,image3d,true=1,quality=100


end
