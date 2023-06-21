pro makebigimage, outfile_name, slice

print, 'opening file'

file_id = H5F_OPEN('DD_Dark_Matter_Density_0016')

print, 'opening dataset'

dataset_id=H5d_open(file_id,'/Dark_Matter_Density')

print, 'reading dataset'

data = h5d_read(dataset_id)

print, 'done reading dataset, closing file'

h5d_close, dataset_id
h5f_close, file_id

print, 'logging data'

logdata = alog10(data + 0.1)

print, 'getting slice'

slice = logdata(1:1024,*,*)

print, 'making plot'

set_plot,'z'
device,set_resolution=[1024,1024]

data_dim = size(slice)
print, 'image info: ',data_dim

min1 = MIN(slice)
max1 = MAX(slice)

print, 'min, max = ', min1, max1

loadct, 0

TVSCL, slice

image = TVRD()
image3D = bytarr(3,1024,1024)
tvlct,r,g,b,/get

image3d[0,*,*] = r[image]
image3d[1,*,*] = g[image]
image3d[2,*,*] = b[image]

write_jpeg,outfile_name,image3d,true=1,quality=100

END
