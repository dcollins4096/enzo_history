pro show_projection, file_name, sds_num, RANGE = range

  icolor_bar = 0
  IF (sds_num lt 0) THEN BEGIN
      sds_num = -sds_num
      icolor_bar = 1
  ENDIF
 
  DFSD_SETINFO, /RESTART
  DFSD_GETINFO, file_name, NSDS=NumSDS

  IF (sds_num gt numSDS) THEN BEGIN
      print, ' read_hdf: Cannot read sds#', sds_num
      print, ' there are only ', numsds,' scientific data sets in ', file_name
      print, ' returning zero valued fltarr'
      data = fltarr(3,3,3)
      RETURN
  ENDIF

  FOR i=1,sds_num DO BEGIN
      r=0.
      l=0.
      HDF_DFSD_GETINFO, file_name, label = l, dims=d, type=t,range=r
;      print,i, l
  ENDFOR
  print, ' read ',l, ' data from ', file_name
  hdf_DFSD_GETDATA, file_name, data

data_dim = size(data)
print, 'image info: ',data_dim

if (sds_num eq 1) then data2 = ALOG10(data)
if (sds_num eq 2) then data2 = ALOG10(data)
if (sds_num eq 3) then data2 = ALOG10(data+1.0e12)
if (sds_num eq 4) then data2 = data
if (sds_num eq 5) then data2 = data
if (sds_num eq 6) then data2 = data
if (sds_num eq 7) then data2 = data
if (sds_num eq 8) then data2 = data
if (sds_num eq 9) then data2 = ALOG10(data+1.0e9)

min1 = MIN(data2)
max1 = MAX(data2)

n1 = data_dim(1)
n2 = data_dim(2)

print, 'min, max =', MIN(data2), MAX(data2)

IF (icolor_bar eq 1) THEN BEGIN
  nadd = 40
  n3 = n2+nadd
  data3 = FLTARR(n1, n3)
  data3(*, *) = min1
  print, size(data3)
  data3(0:n1-1, nadd:n3-1) = data2

  FOR i=1, n1 DO BEGIN
    data3(i-1, 0:nadd-5) = min1 + (i-0.5)*(max1-min1)/n1
  ENDFOR

  data2 = data3

ENDIF

;WINDOW, XSIZE=n1, YSIZE=n2
TVSCL, data2

;show3, data2
;surface, data2
;contour, data2

;write_gif, 'image.gif', TVRD()


;if (sds_num eq 3) then $
;  TVSCL, ALOG10(REBIN(data>1.e12, data_dim(1)*2,data_dim(2)*2, /SAMPLE))+12. $
;  else $
;  TVSCL, ALOG10(REBIN(data, data_dim(1)*2,data_dim(2)*2, /SAMPLE))

END

