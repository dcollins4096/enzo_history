pro STR_TO_FLT_ARR, string, flt_arr
; extract floating point array from a string that are separated by spaces
flt_arr = FLOAT(STR_SEP(strcompress(string),' '))
END
pro STR_TO_INT_ARR, string, int_arr
; extract integers from a string that are separated by spaces
;remove leading and trailing spaces
int_arr = FIX(STR_SEP(strcompress(string), ' '))
END
pro READ_GRID_INFO, file_name, grid_info, REDSHIFT=REDSHIFT
; extract interesting information from hierarchy file 

; 3d data:
  data_dim  = 3
  max_num_of_grids = 1000
  grid_num_max = 0
  grid_level = 0
  info = {grid, $
          num:          0,   $
          hier_file:   '',   $
          level:        0,   $
          dim:          INTARR(data_dim),     $
          Start_index:  INTARR(data_dim),     $
          End_index:    INTARR(data_dim),     $
          Left_edge:    DBLARR(data_dim),     $
          Right_edge:   DBLARR(data_dim),     $
          num_baryon_fields: 0 ,              $
          num_particle: LONG(0) ,              $
          baryon_file:  '',              $
          particle_file:''               $
         }
  temp_grid_info  = Replicate(info, max_num_of_grids) 
  file_in_string = ''
  help_string    = ''
  Grid_string    = 'Grid ='
  print, 'read grid info from ', file_name
 h_num_baryon_fields = 0
  h_baryon_file = ''
  i = 0
  OPENR, unit, file_name, /get_lun
  WHILE (NOT eof(unit)) DO begin
      i = i+ 1
;      print, i, found_flag, help_string
      readf, unit, help_string
      str_parts  = STR_SEP(help_string, ' = ')
; remove leading and trailing zeros'
      key_word   = STRTRIM(str_parts(0),2)
      IF (key_word eq 'Grid') THEN BEGIN
; extract grid number
          value      = STRTRIM(str_parts(1),2)
          grid_num   = FIX(value)
          grid_num_max = max(grid_num, grid_num_max)
;          print, 'Grid number ', grid_num
          WHILE ((strlen(help_string) ge 1) and(NOT eof(unit))) DO begin
; digest AMR terminology
              i = i+ 1
              str_parts  = STR_SEP(help_string, ' = ')
; remove leading and trailing zeros'
              key_word   = STRTRIM(str_parts(0),2)
              value      = STRTRIM(str_parts(1),2)
              CASE key_word OF 
;                  'GridRank'      : Grid_rank  = FIX(value) 
                  'GridDimension' : str_to_int_arr, value, h_dim
                  'GridStartIndex': str_to_int_arr, value, h_Start_index
                  'GridEndIndex'  : str_to_int_arr, value, h_End_index
                  'GridLeftEdge'  : str_to_flt_arr, value, h_Left_edge
                  'GridRightEdge' : str_to_flt_arr, value, h_Right_edge
                  'NumberOfBaryonFields': h_num_baryon_fields = value
'BaryonFileName': h_baryon_file = value
                  'NumberOfParticles' : h_num_particle  = LONG(value)
                  'ParticleFileName'  : h_particle_file = value
                  ELSE: 
              ENDCASE
              readf, unit, help_string
          ENDWHILE
; store grid parameters 
          temp_grid_info(Grid_num-1) =                    $
            {grid, grid_num, file_name,  grid_level, h_dim,           $
             h_Start_index, h_End_index,                  $
             h_Left_edge, h_Right_edge,                   $
             h_num_baryon_fields, h_num_particle,         $
             h_baryon_file, h_particle_file}
      ENDIF      
  ENDWHILE
  CLOSE, unit
  FREE_LUN, unit
; return grid parameters
  Grid_info = temp_grid_info(0:grid_num_max-1)

  dot_is_at = StrPos(file_name, '.hierarchy')
  file_name = STRMID(file_name, 0, dot_is_at)

; if asked for get redshift
  if KEYWORD_SET(REDSHIFT) THEN BEGIN
      OPENR, unit, file_name, /get_lun
      WHILE (NOT eof(unit)) DO begin
          readf, unit, help_string
          str_parts  = STR_SEP(help_string, ' = ')
; remove leading and trailing zeros'
          key_word   = STRTRIM(str_parts(0),2)
          if (key_word eq 'CosmologyCurrentRedshift') then $
            REDSHIFT  = FLOAT(STRTRIM(str_parts(1),2))
      ENDWHILE
      CLOSE, unit
      FREE_LUN, unit
  ENDIF


  RETURN
END

PRO  hdf_info, file_name
; print out info on the sds in the specified hdf file
  l = ' '
  dummy = FindFile(file_name , COUNT = exist)  
  if dummy(0) eq '' then BEGIN
      print, 'read_hdf: the file "', file_name, '" does not exist!'
      RETURN
  ENDIF
  sd_id = HDF_SD_START(file_name, /READ)
  HDF_SD_FILEINFO, sd_id, NumSDS, attributes
  print, '# of global attr:', attributes
  n = 0
  FOR i=1,numsds-1 DO BEGIN
      r=0.
      l=0.
      sds_id = HDF_SD_SELECT(sd_id, n)
      WHILE (HDF_SD_ISCOORDVAR(sds_id)) DO BEGIN
          HDF_SD_ENDACCESS, sds_id
          n = n + 1
          i = i + 1
          sds_id = HDF_SD_SELECT(sd_id, n)
      ENDWHILE
      HDF_SD_GETINFO, sds_id, name = nam, label = l, dims=d, type=t
      print, i,' ', nam, STRCOMPRESS(STRING(d))+' '+STRING( t)
      n = n + 1
      HDF_SD_ENDACCESS, sds_id
  ENDFOR
  HDF_SD_ENDACCESS, sds_id
  HDF_SD_END, sd_id
  RETURN
END

PRO read_particle_data, file_name, sds_start_index, N_p, data, data_name
; start_index is the sds starting index where to search for particle data
; this is needed since enzo is not consistent with his output !
  l = ' '
  dummy = FindFile(file_name , COUNT = exist)  
  if dummy(0) eq '' then BEGIN
      print, 'read_hdf: the file "', file_name, '" does not exist!'
      RETURN
  ENDIF
  sd_id = HDF_SD_START(file_name, /READ)
  HDF_SD_FILEINFO, sd_id, NumSDS, attributes
;  print, '# of global attr:', attributes
  sds_id = HDF_SD_SELECT(sd_id, sds_start_index)
  HDF_SD_GETINFO, sds_id, name = data_name
  WHILE ((strpos(data_name, 'ticle') le 0) and $
         (sds_start_index lt numsds)) DO  BEGIN
      HDF_SD_ENDACCESS, sds_id
      sds_start_index=sds_start_index+1
      sds_id = HDF_SD_SELECT(sd_id, sds_start_index)
      HDF_SD_GETINFO, sds_id, name = data_name
  ENDWHILE
      
;  print, sds_start_index, data_name
  HDF_SD_GETDATA, sds_id, t_data
  HDF_SD_ENDACCESS, sds_id
  ds = size(t_data)
  data = FLTARR(8,ds(1)) 
  data(0,*) = t_data(*)
  for j=0,6 DO BEGIN
      sds_id = HDF_SD_SELECT(sd_id, sds_start_index+j)
      HDF_SD_GETINFO, sds_id, name = data_name
      print,  sds_start_index+j, data_name
      HDF_SD_GETDATA, sds_id, t_data
      HDF_SD_ENDACCESS, sds_id
      data(j,*) = t_data(*)
  ENDFOR      
  N_P = N_elements(data(0,*))
  HDF_SD_END, sd_id
  RETURN
END

PRO write_tipsy, out_file, data, red 
; dump particle data in tipsy format
; open outputfile
  OPENW, unit, out_file, /get_lun
; write header  
  ds = size(data)
  n_count = ds(2)
  printf, unit, n_count, ' ',0, ' ', 0
; dimensions
  printf, unit, 3
  printf, unit ,red
; part. masses 
  printf, unit, (data(6,*))
  for i=0,5 DO BEGIN
;pos and vel's       
      printf, unit, (data(i,*))
  ENDFOR
  CLOSE, unit
  FREE_LUN, unit

RETURN
END

pro amr_to_tipsy, hier_file, data_d,out_file
; first read in all necessary data
  red = 1
  read_grid_info, hier_file, grid_info, REDSHIFT= red
  
; do not use all grids !
;  grid_info = grid_info(2:*)

  n_total = LONG(TOTAL(grid_info(*).num_particle))
  print, 'total number of particles in simulation:',n_total
  dat = FLTARR(8,n_total)
  N_count = 0
  FOR i=0,N_elements(grid_info)-1 DO BEGIN
      IF grid_info(i).num_particle gt 0 then BEGIN
          read_particle_data,data_d+grid_info(i).particle_file, $
            0, N_p, data, data_label
          dat(*,N_count:N_count+N_p-1) = data(*,*)
          print, 'read ',N_p,' particles from ',data_d+grid_info(i).particle_file
          N_count = N_count+N_p
      ENDIF
  ENDFOR
  print, 'number check:',N_count, N_total

; sort and check uniquness of indices
  ind = (SORT(dat(7,*)))
  print, 'there are ', N_elements(ind), ' where ', N_elements(uniq(ind)),$
    ' are unique !'
  
  write_tipsy, out_file, dat, red
  print, 'written :', out_file

RETURN
END


num = 7
dnum = STRCOMPRESS(STRING(num,FORMAT='(i4.4)'), /REMOVE_ALL)
;dat_d = './f128L12.2/time'+dnum+'/'
dat_d = '/scratch-modi4/tabel/NC1/'
file = 'RedshiftOutput'+dnum+'.grid0001'
hier_file = dat_d + 'RedshiftOutput'+dnum+'.hierarchy'
out_file = dat_d+'t'+dnum+'p.ascii'

; get info on an hdf dump:
;hdf_info, dat_d+file
;for gn=1,15 do hdf_info, dat_d+'128_halo1_00'+dnum+'.grid0'+STRCOMPRESS(STRING(gn, FORMAT='(I3.3)'), /REMOVE_ALL)


;sds_index = 63 ; starting index to search for particle data
; read_particle_data,dat_d+file , sds_index, N_p, data, data_name

amr_to_tipsy, hier_file, dat_d, out_file

END

; run this routine with 
;.run amrtotipsy.com
