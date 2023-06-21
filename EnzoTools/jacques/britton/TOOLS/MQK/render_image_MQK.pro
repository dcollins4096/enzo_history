;; originally written by: Tom Abel
;; modifications by:      Michael Kuhlen
;;
;; begun:                 02/17/2004
;; last modified:         12/18/2004
;;
;; Description:
;;
;; + this routine used to be in the file
;;   construct_interpolated_slice_data.pro, but has been moved to its
;;   own file to allow for automatic compilation.

pro render_image, result, secmin=secmin, missing_value=missed
@TOOLS/common_blocks.inc
  interpolate =  interpolate_i
; interpolate all values
; get triangles:

  if (not keyword_set(secmin)) then secmin = min(z_stuetz)
  
  if verbose then  print, 'max before gridding:',max(z_stuetz), min(z_stuetz)
  
;  print, 'b:',b
  result = dblarr(xy_sl_size,xy_sl_size)
  if n_elements(missed) gt 0 then begin
      if missed then begin 
          result[*] = 0. 
          if verbose then print, 'render_image: set default value to 0.'
      end else begin
          result[*] = secmin
          if verbose then print, 'render_image: set default value to ', secmin
      endelse
  endif
  
;   gs =  double([(slice_size[2]-slice_size[0])/double(xy_sl_size), $
;          (slice_size[3]-slice_size[1])/double(xy_sl_size)])
  
  
  limits =  double([0., 0., 1., 1.])
  if (interpolate ne 0) then begin
      gs =  double([1./double(xy_sl_size), 1./double(xy_sl_size)])
      tx_stuetz =     double((x_stuetz-slice_size[0])/(slice_size[2]-slice_size[0])) 
      ty_stuetz =     double((y_stuetz-slice_size[1])/(slice_size[3]-slice_size[1])) 
      
      ind =  where((tx_stuetz ge 0. ) and (tx_stuetz le 1. ) and (ty_stuetz ge 0. ) and (ty_stuetz le 1. ))
      tx_stuetz =    tx_stuetz[ind]
      ty_stuetz =    ty_stuetz[ind]
      tz_stuetz  =    z_stuetz[ind]
      triangulate, tx_stuetz, ty_stuetz, triangles, b
  endif
  
  if (interpolate eq 2) then result = $
    trigrid(tx_stuetz,ty_stuetz, double(tz_stuetz), triangles, input=result, $
            gs, limits, /quintic, min_value=secmin) <  max(z_stuetz)
  if (interpolate eq 1) then result = $
    trigrid(tx_stuetz,ty_stuetz, double(tz_stuetz), triangles, input=result, $
            gs, limits, min_value=secmin)  <  max(z_stuetz)
  if (interpolate eq 0) then begin 
; no interpolation:  
      dxts =  double(xy_sl_size)/(double(slice_size[2])-double(slice_size[0]))
      dyts =  double(xy_sl_size)/(double(slice_size[3])-double(slice_size[1]))
      
;     print,'dxts: ',dxts
;     print,'dyts: ',dyts
;     print,'n_stuetz: ',n_stuetz
;     print,'slice_size: ',slice_size
;     print,'xy_sl_size: ',xy_sl_size
      
      for i=0L,n_stuetz-1L do begin
          xmi =  round((x_stuetz[i]-0.5D*dx_stuetz[i]-double(slice_size[0]))*dxts) >  0  < (xy_sl_size-1)
          xma =  round((x_stuetz[i]+0.5D*dx_stuetz[i]-double(slice_size[0]))*dxts) < (xy_sl_size-1) > 0
          ymi =  round((y_stuetz[i]-0.5D*dx_stuetz[i]-double(slice_size[1]))*dyts) > 0 < (xy_sl_size-1)
          yma =  round((y_stuetz[i]+0.5D*dx_stuetz[i]-double(slice_size[1]))*dyts) < (xy_sl_size-1) > 0
          
;        print,xmi,xma,ymi,yma
          if ((ymi ge yma) or (xmi ge xma)) then begin
              result(xmi,ymi) = z_stuetz[i]
          endif else begin
              hh =  dblarr((xma-xmi+1),(yma-ymi+1))
              hh =  double(z_stuetz[i])
              
              result(xmi:xma,ymi:yma) = hh
          endelse
      endfor
      
  endif
  
  if verbose then print, 'max after griding:',max(result), min(result)
  
  return
end

