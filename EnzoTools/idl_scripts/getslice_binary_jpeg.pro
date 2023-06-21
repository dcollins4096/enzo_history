pro getslice, file_name, outfile_name, size

openr,1,file_name
data = fltarr(size,size);
readu,1,data
close,1

logdata = alog10(data+1.0)


;bigslice = rebin(smallslice,1,512,512)

;window, xsize=size,ysize=size

print, 'min,max = ',min(logdata),max(logdata)



set_plot,'z'
device,set_resolution=[size,size]

loadct,4

TVSCL, logdata
image=TVRD()
image3D = bytarr(3,size, size) 
tvlct,r,g,b,/get 


image3d[0,*,*] = r[image]
image3d[1,*,*] = g[image]
image3d[2,*,*] = b[image]

write_jpeg,outfile_name,image3d,true=1,quality=100


END 
