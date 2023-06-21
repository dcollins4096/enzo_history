;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/17/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + uses IDL's built-in HDF5 routines
;; + NOTE: there exists a compatibility problem between HDF5 v.1.4.x
;;   and v.1.6.x. Versions of IDL less than v6.1 use a shared library
;;   compiled with HDF5 1.4.3. If your ENZO output dumps were written
;;   with HDF5 v.1.6.x, then you will need to use IDL v.6.1.
;; + read grid by dataset_name, not sds_num
;; + first look in grid pointer array, to see whether this grid has
;;   already been loaded.

pro read_grid,file_name,data,$
             byname=dataset_name
; read hdf file with all data matching field_names of each grid 
; and return it in data_structure
;	See if there is anything there to read


  common gridstorage
  common options


  if not keyword_set(dataset_name) then $
    message,'read_grid: must set keyword BYNAME!'

;; check to see that file exists
  dummy=findfile(file_name,count=exist)  
  if not exist then $
    message,'read_grid: the file "'+file_name+'" does not exist!'



;; check to see whether grid has been loaded

  ; first get gridnumber from name
  pos=strpos(file_name,'.grid')
  gridindex=long(strmid(file_name,pos+5))-1l
  
  message,/info,'reading field "'+dataset_name+'" in grid #'+string(gridindex+1,format='(i4.4)')+'...'

  ; pointer defined? (i.e. has grid been read before?)
  if ptr_valid(allgrids[gridindex]) then begin

      ; correct dataset?
      if ((*allgrids[gridindex]).fieldname eq dataset_name) then begin

          if verbose then $
            message,/info,'Grid '+string(gridindex,format='(i0)')+' is still in memory.'
          data=(*allgrids[gridindex]).grid
          return
      endif else begin 
          ptr_free,allgrids[gridindex]
      endelse
  endif

  if verbose then $
    message,/info,'Grid '+string(gridindex,format='(i0)')+' is not in memory.'

  file_id=h5f_open(file_name)
  dataset_id=h5d_open(file_id,dataset_name)

  data=h5d_read(dataset_id)

  h5d_close,dataset_id
  h5f_close,file_id

  ;; check if found dataset
  if n_elements(data) eq 0 then $
    message,"read_grid: couldn't find dataset!"

  
;; store grid in allgrids array
  
  ; first create structure
  thisgrid={fieldname:dataset_name,$
            grid:data}
  ; assign pointer
  allgrids[gridindex]=ptr_new(thisgrid,/no_copy)

  return
end

