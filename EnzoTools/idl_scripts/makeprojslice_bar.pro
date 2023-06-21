pro makeprojslice_bar, infile_name, outfile_proj, outfile_slice, slicenum

print, 'input:           ',infile_name
print, 'projection file: ',outfile_proj
print, 'slice file:      ',outfile_slice
print, 'projection #:    ',slicenum

print, 'opening file ', infile_name

file_id = H5F_OPEN(infile_name)

print, 'opening dataset'

dataset_id=H5d_open(file_id,'/Density')

print, 'reading dataset'

data = h5d_read(dataset_id)

print, 'done reading dataset, closing file'

h5d_close, dataset_id
h5f_close, file_id

print, 'making projection'

projection = data(0,*,*)

for ii=1, 1023 do begin

;print, ii

projection = projection + data(ii,*,*)

endfor

print, 'making slice'

slice = data(slicenum,*,*)

slice = slice/1024.0
projection = projection/1024.0

print, 'logging data'

logprojection = alog10(projection + 0.00001)
logslice = alog10(slice + 0.00001)

print, 'writing files: ', outfile_proj, outfile_slice

; preemptive close
close, 1
close, 2

openw, 1, outfile_proj
openw, 2, outfile_slice

writeu, 1, logprojection
writeu, 2, logslice

close, 1
close, 2

END
