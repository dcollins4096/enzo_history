pro READ_GRID_INFO, data_dir, file_name, grid_info, $
        list_str, mass_per_particle, divide_by_density, REDSHIFT=red
; extract interesting information from hierarchy file 

; 3d data:
  data_dim  = 3
  max_num_of_grids = 9000
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
  print, 'reading grid info from ', file_name
  i = 0
  h_particle_file = ""
  h_baryon_file   = ""
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
                  'GridLeftEdge'  : str_to_dbl_arr, value, h_Left_edge
                  'GridRightEdge' : str_to_dbl_arr, value, h_Right_edge
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

; remove .hierarchy
dot_is_at = RStrPos(file_name, '.')
file_name = STRMID(file_name, 0, dot_is_at)

slash_is_at =  RStrPos(file_name, '/')
file_name =  STRMID(file_name, slash_is_at+1,100)
;May be the hierarchy file is also in parent directory:
exist = 0
dummy = FindFile((data_dir+'/'+file_name + '.dir'), COUNT = exist)
IF (exist GT 0) THEN data_dir =  data_dir+'/'+file_name+'.dir/'
file_name = data_dir + file_name

red =  -1.
; if asked for get redshift
if ARG_PRESENT(red) THEN BEGIN
   OPENR, unit, file_name, /get_lun
   WHILE (NOT eof(unit)) DO begin
      readf, unit, help_string
      str_parts  = STR_SEP(help_string, ' = ')
; remove leading and trailing zeros'
      key_word   = STRTRIM(str_parts(0),2)
      if (key_word eq 'CosmologyCurrentRedshift') then $
       red  = FLOAT(STRTRIM(str_parts(1),2))
   ENDWHILE
   CLOSE, unit
   FREE_LUN, unit
   print, 'REDSHIFT:', red, ' from file:', file_name
ENDIF

; get datalabels, etc..
fnum = grid_info(1).num_baryon_fields + 5
mass_per_particle = FLTARR(fnum)
divide_by_density = BYTARR(fnum)
list_str = strarr(fnum+1)
dummy = FindFile((file_name + '.'+'grid0010'), COUNT = exist)
if exist eq 0 then file_name = file_name + '.grid0001' else $
  file_name = file_name + '.grid0010'
; read data labels from hdf data file 
Print, 'Extracting field names, etc from: ', file_name
;DFSD_SETINFO, /RESTART
;DFSD_GETINFO, file_name, NSDS=NumSDS
sd_id = HDF_SD_START(file_name, /READ)
HDF_SD_FILEINFO, sd_id, NumSDS, attributes
i=0
FOR ifake=0,NumSDS-9 DO BEGIN
    r=0.
    l= '' 
    sds_id = HDF_SD_SELECT(sd_id, ifake)
    IF (NOT HDF_SD_ISCOORDVAR(sds_id)) THEN BEGIN
        i=i+1
        HDF_SD_GETINFO, sds_id, label = l, dims=d
        list_str(i-1) = l
        print,i, list_str(i-1)
        mass_per_particle(i-1) = 1.0
        divide_by_density(i-1) = 0 
        IF (STRPOS(l, '_Density') gt 0) THEN BEGIN
            print, 'contains _Density'
            list_str(i-1) = STRMID(list_str(i-1),0,STRPOS(l, '_De'))+'_fraction'
            divide_by_density(i-1) = 1
            IF (STRPOS(l, 'HeI') ge 0) THEN BEGIN 
                mass_per_particle(i-1) = 4.0
            ENDIF
            IF ((strpos(l, 'HI') ge 0) OR (strpos(l,'HM') ge 0)) THEN BEGIN 
                mass_per_particle(i-1) = 1.0
            ENDIF
            IF ((strpos(l, 'H2I') ge 0)) THEN BEGIN 
                mass_per_particle(i-1) = 2.0
            ENDIF
        ENDIF
        print, list_str(i-1), mass_per_particle(i-1), divide_by_density(i-1)
    ENDIF
    HDF_SD_ENDACCESS, sds_id
ENDFOR
; add derived quantities:
print, 'read_grid_info: added angular momentum to list'
list_str(fnum-3) =  'Angular Momentum'
print, 'read_grid_info: added Entropy to list'
list_str(fnum-2) =  'Entropy'
print, 'read_grid_info: added Pressure to list'
list_str(fnum-1) =  'Pressure'

HDF_SD_END, sd_id


order_grids, grid_info

RETURN
END

;.compile TOOLS/read_grid_info
