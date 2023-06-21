pro show_peaks
; show peaks projected on the 2d slice you are viewing
@TOOLS/common_blocks.inc

  axes       = (slice_ori gt 1.e-30)
  const_sub  = where(axes gt 1e-10)
  other_subs = where(axes eq 0)

  FOR i=0, N_peaks-1 DO BEGIN
      R_p = peaks(0:2,i)
      print, 'r_p:',r_p, slice_size, xy_sl_size
      R_p = ROUND((R_p(other_subs)-slice_size(0:1))*$
                  DOUBLE(xy_sl_size)/(slice_size(2)-slice_size(0)))
      print, 'pix:', R_p
      xyouts, R_p(0), R_p(1), STRING(i), $
        CHARSIZE = 0.8, /DEVICE, ALIGNMENT = 1.0
      PLOTS, R_p(0)+[-4,4], R_p(1)+[-4,4], /DEVICE
      PLOTS, R_p(0)+[4,-4], R_p(1)+[-4,4], /DEVICE
  ENDFOR
  
RETURN
END
