PRO read_peaks, file_name, peaks, N_peaks, max_peak_nr

; Open the file to read
  OPENR, unit, file_name, /GET_LUN, ERROR=error

; If an error occurred
; Notify the User that an error occurred
  IF error LT 0 THEN BEGIN
      PRINT, ' Cannot open the file ' + file_name
      PRINT, ' I will keep the old peaks !'
      RETURN
  ENDIF

; read file
  i = 0
  nst = strarr(100)
  while (not EOF(unit)) DO BEGIN
      line = ''
      readf,unit,line
; lines that contain a #  are comment lines !
      if (strpos(line,'#') eq -1) then BEGIN
          nst(i) = strtrim(strcompress(line),2)
          i = i + 1
      ENDIF ELSE PRINT, 'line',i,' is a comment line'
  END
  close, unit
; use only available peaks  
  nst = nst(0:((i-1)>0))
; the first three columns are x,y, and z positions
  column_num = 3
  peak = FLTARR(column_num, i)   ;
  k = 0 
  FOR j=0,i-1 DO BEGIN
      line_f    = FLOAT(str_sep(nst(j),' '))
      if N_elements(line_f) lt 3 then problem = 1 else BEGIN
          peak(*,j) = line_f(0:2)          
          problem  = ((max(peak(*,j)) gt 1.) or min(peak(*,j) lt 0.))
      ENDELSE
      if problem THEN BEGIN
          k = k + 1 
          PRINT, 'Problem with peak#',j,'in ', file_name
          PRINT, 'peak:', peak(*,j)
          PRINT, 'will not use it !'
          peak(*,j) = -1.
      ENDIF
  ENDFOR

  if (k lt i) THEN BEGIN
      peaks = FLTARR(3,i-k)
      hj = 0
      FOR j=0,i-1 DO BEGIN
          if (min(peak(0:2,j)) gt -1. ) THEN BEGIN
              peaks(0:2,hj) = peak(0:2, j) 
              hj = hj + 1 
          ENDIF
      ENDFOR
  END ELSE BEGIN 
      PRINT, ' no valid peaks'
      PRINT, ' will continue with old ones'
  ENDELSE 
  N_peaks = hj
  IF N_PEAKS ge max_peak_nr THEN BEGIN
      PRINT, ' maximally ',max_peak_nr, 'are supported' 
      PRINT, 'I will only use the first',max_peak_nr,' from ', file_name
      PRINT, 'Use multiple peak files to get around this'
      N_peaks = max_peak_nr
      peaks  = peaks(0:2,0:N_peaks-1)
  END 


  RETURN
END
