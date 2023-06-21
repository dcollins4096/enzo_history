;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/16/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + uses IDL's built-in HDF5 routines.
;; + NOTE: there exists a compatibility problem between HDF5 v.1.4.x
;;   and v.1.6.x. Versions of IDL less than v6.1 use a shared library
;;   compiled with HDF5 1.4.3. If your ENZO output dumps were written
;;   with HDF5 v.1.6.x, then you will need to use IDL v.6.1.
;; + reads all datasets in grid0001 and prompts user to choose a a
;;   set. This was done because some of my dumps have more datasets
;;   than others, i.e. DM particle info.
;; + changed all ints to longs. (IDL ints are 2-byte!)
;; + upped maximum number of grids to 50000l
;; + use base_name instead of file_name as argument. Construct
;;   hierarchy filename and grid filenames from base_name.




pro read_grid_info, data_dir, $
                    base_name, $ ; input
                    grid_info, $           ; output
                    list_str, $            ; output
                    mass_per_particle, $   ; output
                    divide_by_density, $   ; output
                    redshift=red           ; input

common packAMR

; extract interesting information from hierarchy file 

; 3d data:
  data_dim  = 3
  max_num_of_grids = 50000l
  grid_num_max = 0l
  grid_level = 0l
  info = {grid, $
          num:          0l,   $
          hier_file:   '',   $
          level:        0l,   $
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
  i = 0l
  h_particle_file = ""
  h_baryon_file   = ""


  hierarchy_file_name = base_name+'.hierarchy'
  message,/info,'reading grid info from '+hierarchy_file_name

  openr, unit, hierarchy_file_name, /get_lun
  while (not eof(unit)) do begin
      i = i+ 1l
;      print, i, found_flag, help_string
      readf, unit, help_string
      str_parts  = str_sep(help_string, ' = ')
; remove leading and trailing zeros'
      key_word   = strtrim(str_parts(0),2)
      if (key_word eq 'Grid') then begin
; extract grid number
          value      = strtrim(str_parts(1),2)
          grid_num   = long(value)
          grid_num_max = max(grid_num, grid_num_max)
;          print, 'Grid number ', grid_num
          while ((strlen(help_string) ge 1) and(not eof(unit))) do begin
; digest AMR terminology
              i = i+ 1l
              str_parts  = str_sep(help_string, ' = ')
; remove leading and trailing zeros'
              key_word   = strtrim(str_parts(0),2)
              value      = strtrim(str_parts(1),2)
              case key_word of 
;                  'GridRank'      : Grid_rank  = LONG(value) 
                  'GridDimension' : str_to_int_arr, value, h_dim
                  'GridStartIndex': str_to_int_arr, value, h_Start_index
                  'GridEndIndex'  : str_to_int_arr, value, h_End_index
                  'GridLeftEdge'  : str_to_dbl_arr, value, h_Left_edge
                  'GridRightEdge' : str_to_dbl_arr, value, h_Right_edge
                  'NumberOfBaryonFields': h_num_baryon_fields = value
                  'BaryonFileName': h_baryon_file = data_dir + value
                  'NumberOfParticles' : h_num_particle  = LONG(value)
                  'ParticleFileName'  : h_particle_file = data_dir + value
                  else: 
              endcase

              readf, unit, help_string
          endwhile


; store grid parameters 
          temp_grid_info[Grid_num-1] =                    $
            {grid, grid_num, hierarchy_file_name,  grid_level, h_dim,           $
             h_Start_index, h_End_index,                  $
             h_Left_edge, h_Right_edge,                   $
             h_num_baryon_fields, h_num_particle,         $
             h_baryon_file, h_particle_file}
      endif      
  endwhile
  close, unit
  free_lun, unit
; return grid parameters
  Grid_info = temp_grid_info(0:grid_num_max-1)


  red =  -1.
; if asked for get redshift
  if arg_present(red) then begin
      openr, unit, base_name, /get_lun
      while (not eof(unit)) do begin
          readf, unit, help_string
          str_parts  = STR_SEP(help_string, ' = ')
; remove leading and trailing zeros'
          key_word   = STRTRIM(str_parts(0),2)
          if (key_word eq 'CosmologyCurrentRedshift') then $
            red  = DOUBLE(STRTRIM(str_parts(1),2))
      endwhile
      close, unit
      free_lun, unit
      print, 'redshift:', red, ' from file:', base_name
  endif
  

;; get datalabels, etc.

  grid_file_name=base_name+'.grid0001'
  dummy=findfile(grid_file_name,count=exist)
  if not exist then begin
      print,grid_file_name+' does not exist!'
      grid_file_name=base_name+'.cpu0000'
      dummy=findfile(grid_file_name,count=exist)
      if not exist then begin
          message,grid_file_name+' does not exist, either!'
      endif else begin
          print,grid_file_name+' does exist!'
          print,'Using packed amr.'
          amr_is_packed = 1
      endelse
  endif else begin
      print, 'Using unpacked amr.'
      amr_is_packed = 0
  endelse
  
; read data labels from hdf data file 
  print, 'Extracting field names, etc from: ',grid_file_name
  
  file_id=h5f_open(grid_file_name)
  Nf=h5g_get_nmembers(file_id,'/')

  ; test if this is a packed amr file
  checkPack = h5g_get_member_name(file_id,'/',0)
  checkPackPrefix = strmid(checkPack,0,4)

  if (checkPackPrefix eq 'Grid') then begin
      ; this is packed amr
      firstGroup = '/' + checkPack
      Nf = h5g_get_nmembers(file_id,firstGroup)
  endif else begin
      ; this is old-style grid file
      firstGroup = '/'
  endelse

  all_labels=strarr(Nf)
  for i=0l,Nf-1 do begin
      all_labels[i]=h5g_get_member_name(file_id,firstGroup,i)
  endfor
  h5f_close,file_id
  
  
;; remove all particle fields (beginning with 'PARTICLE_')
  p1=strmid(all_labels,0,9)
  ind=where(p1 ne 'particle_')
  if ind[0] ne -1 then all_labels=all_labels[ind]
  
  
  Nf=n_elements(all_labels)
  list_str=strarr(Nf+6)
  list_str[0:Nf-1]=all_labels
  list_str[Nf+2]='Angular Momentum'
  list_str[Nf+3]='Entropy'
  list_str[Nf+4]='Pressure'
  mass_per_particle=dblarr(Nf+6)+1.0
  divide_by_density=bytarr(Nf+6)
  

 ;; convert species densities in fractions

  for i=0,Nf-1 do begin
      if (strpos(list_str[i],'_Density') gt 0) then begin
;          print, 'contains _density'

;; PROBLEM: cannot substitute '_fraction' for '_Density' at this
;; point, since the dataset is read out by using the list_str name.

;          list_str[i]=strmid(list_str[i],0,strpos(list_str[i],'_De'))+'_fraction'

          divide_by_density[i]=1
          if(strpos(list_str[i],'HeI') ge 0) then $
            mass_per_particle[i]=4.0
          if((strpos(list_str[i],'HI') ge 0) or $
             (strpos(list_str[i],'HM') ge 0)) then $
            mass_per_particle(i-1) = 1.0
          if(strpos(list_str[i],'H2I') ge 0)then $
            mass_per_particle(i-1) = 2.0
      endif
  endfor

; print,list_str
; print,mass_per_particle
; print,divide_by_density
  
order_grids, grid_info

return
end

;.compile TOOLS/read_grid_info
