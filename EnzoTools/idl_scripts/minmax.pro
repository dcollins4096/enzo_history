pro minmax, infile

print, 'opening file'

file_id = H5F_OPEN(infile)

print, 'opening dataset'

dataset_id=H5d_open(file_id,'/Density')

print, 'reading dataset'

data = h5d_read(dataset_id)

print, 'done reading dataset, closing file'

h5d_close, dataset_id
h5f_close, file_id

min1 = MIN(data)
max1 = MAX(data)

print, 'min, max = ', min1, max1


END
