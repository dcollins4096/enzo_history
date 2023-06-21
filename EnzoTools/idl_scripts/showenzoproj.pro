pro showenzoproj, file_name, outfile_name, number

file_id = H5F_OPEN(file_name)

if (number eq 0) then dataset_id=H5D_OPEN(file_id,'/Dark_Matter_Density')
if (number eq 1) then dataset_id=H5D_OPEN(file_id,'/Density')
if (number eq 2) then dataset_id=H5D_OPEN(file_id,'/Electron_Density')
if (number eq 3) then dataset_id=H5D_OPEN(file_id,'/H2II_Density')
if (number eq 4) then dataset_id=H5D_OPEN(file_id,'/H2I_Density')
if (number eq 5) then dataset_id=H5D_OPEN(file_id,'/HII_Density')
if (number eq 6) then dataset_id=H5D_OPEN(file_id,'/HI_Density')
if (number eq 7) then dataset_id=H5D_OPEN(file_id,'/HM_Density')
if (number eq 8) then dataset_id=H5D_OPEN(file_id,'/HeIII_Density')
if (number eq 9) then dataset_id=H5D_OPEN(file_id,'/HeII_Density')
if (number eq 10) then dataset_id=H5D_OPEN(file_id,'/HeI_Density')
if (number eq 11) then dataset_id=H5D_OPEN(file_id,'/Level')
if (number eq 12) then dataset_id=H5D_OPEN(file_id,'/Metal_Density')
if (number eq 13) then dataset_id=H5D_OPEN(file_id,'/SZ_Kinetic')
if (number eq 14) then dataset_id=H5D_OPEN(file_id,'/SZ_Y_Effect')
if (number eq 15) then dataset_id=H5D_OPEN(file_id,'/Star_Density')
if (number eq 16) then dataset_id=H5D_OPEN(file_id,'/Temp_X_Ray_Weighted')
if (number eq 17) then dataset_id=H5D_OPEN(file_id,'/X_Ray_Luminosity')
if (number eq 18) then dataset_id=H5D_OPEN(file_id,'/Z1_Metal_Density')
if (number eq 19) then dataset_id=H5D_OPEN(file_id,'/Z2_Metal_Density')

print, 'test ', file_name, number
print, 'dataset id ', dataset_id

data = H5D_READ(dataset_id)

nummem = H5G_GET_NMEMBERS(file_id,'/')

print, 'there are ',nummem,' members'

membername = H5G_GET_MEMBER_NAME(file_id,'/',number)
print, 'we are opening member number ',membername

dataspace_id = H5D_GET_SPACE(dataset_id)

dimensions = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)

print, 'dimensions: ', dimensions

H5S_CLOSE, dataspace_id
H5D_CLOSE, dataset_id
H5F_CLOSE, file_id


;WINDOW, XSIZE=dimensions[0], YSIZE=dimensions[1]

set_plot,'z'
device,set_resolution=[dimensions[0],dimensions[1]]

data_dim = size(data)
print, 'image info: ',data_dim
;window

if (number eq 0) then data2 = ALOG10(data+1.0)
if (number eq 1) then data2 = ALOG10(data+1.0)
if (number eq 2) then data2 = ALOG10(data+1.0)
if (number eq 3) then data2 = ALOG10(data+1.0)
if (number eq 4) then data2 = ALOG10(data+1.0)
if (number eq 5) then data2 = ALOG10(data+1.0)
if (number eq 6) then data2 = ALOG10(data+1.0)
if (number eq 7) then data2 = ALOG10(data+1.0)
if (number eq 8) then data2 = ALOG10(data+1.0)
if (number eq 9) then data2 = ALOG10(data+1.0)
if (number eq 10) then data2 = ALOG10(data+1.0)
if (number eq 11) then data2 = data
if (number eq 12) then data2 = ALOG10(data+1.0)
if (number eq 13) then data2 = data
if (number eq 14) then data2 = data
if (number eq 15) then data2 = ALOG10(data+1.0)
if (number eq 16) then data2 = ALOG10(data+1.0)
if (number eq 17) then data2 = ALOG10(data+1.0)
if (number eq 18) then data2 = ALOG10(data+1.0)
if (number eq 19) then data2 = ALOG10(data+1.0)

min1 = MIN(data2)
max1 = MAX(data2)

print,  'min, max =', min1,max1

TVSCL, data2

image=TVRD()
image3D = bytarr(3,dimensions[0],dimensions[1]) 
tvlct,r,g,b,/get

image3d[0,*,*] = r[image]
image3d[1,*,*] = g[image]
image3d[2,*,*] = b[image]

write_jpeg,outfile_name,image3d,true=1,quality=100


END 
