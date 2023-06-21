PRO  read_hdf, file_name, sds_num, data, l
; read hdf file with all data matching field_names of each grid 
; and return it in data_structure
;	See if there is anything there to read

  l = ' '
  dummy = FindFile(file_name , COUNT = exist)  
  if dummy(0) eq '' then BEGIN
      print, 'read_hdf: the file "', file_name, '" does not exist!'
      print, ' returning zero valued fltarr'
      data = fltarr(3,3,3)
      RETURN
  ENDIF
  sd_id = HDF_SD_START(file_name, /READ)
  HDF_SD_FILEINFO, sd_id, NumSDS, attributes

  IF (sds_num gt numSDS) THEN BEGIN
      print, ' read_hdf: Cannot read sds#', sds_num
      print, ' there are only ', numsds,' scientific data sets in ', file_name
      print, ' returning zero valued fltarr'
      data = fltarr(3,3,3)
      RETURN
  ENDIF

  n = 0
  FOR i=1,sds_num DO BEGIN
      r=0.
      l=0.
      sds_id = HDF_SD_SELECT(sd_id, n)
      WHILE (HDF_SD_ISCOORDVAR(sds_id)) DO BEGIN
        HDF_SD_ENDACCESS, sds_id
        n = n + 1
        sds_id = HDF_SD_SELECT(sd_id, n)
      ENDWHILE
      n = n + 1
      IF (i NE sds_num) THEN HDF_SD_ENDACCESS, sds_id
  ENDFOR
  HDF_SD_GETINFO, sds_id, label = l, dims=d, type=t
;      print,i, l
  print, ' read ',l, ' data from ', file_name
  HDF_SD_GETDATA, sds_id, data
  HDF_SD_ENDACCESS, sds_id
  HDF_SD_END, sd_id

  RETURN
END

