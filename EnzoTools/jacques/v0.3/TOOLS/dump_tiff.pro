PRO dump_tiff, draw_id, draw_cm_id, un_sensitive_id, suggest_title

filename = pickfile(get_path = path, TITLE = 'select name of TIFF file', FILTER = '*.tif')

;widget_control, un_sensitive_id, sensitive =0 
WIDGET_CONTROL, /HOURGLASS
print, 'dumpe_tiff:copying from screen. This might take some time'
WIDGET_CONTROL, draw_id, GET_VALUE=win
WSET, win   
backing = tvrd(TRUE = 1)
WIDGET_CONTROL, draw_cm_id, GET_VALUE=win2
WSET, win2
cm_backing = tvrd(TRUE = 1)
wset, win
 
back_s    =  size(backing)   
cm_back_s =  size(cm_backing)   
bs1 = back_s(2)
bs2 = back_s(3)
cs1 = cm_back_s(2)
cs2 = cm_back_s(3)
print, bs1, bs2, cs1, cs2
comb_im   = bytarr(3, bs1, bs2+cs2)
comb_im(0:2,0:bs1-1, cs2:(bs2+cs2-1)) = backing
comb_im(0:2,0:bs1-1,       0:(cs2-1)) = cm_backing
comb_im_s = size(comb_im)
print, size(comb_im)

;tvlct, r,g,b, /get
write_tiff, filename,comb_im ; , red = r, green = g, blue = b
print, 'done dumping tiff file', filename
RETURN
END

;comb =  fltarr(3, 200,100)
;comb(*,0:99,*) =  a
;comb(*,100:*,*) =  a
; .compile dump_tiff
