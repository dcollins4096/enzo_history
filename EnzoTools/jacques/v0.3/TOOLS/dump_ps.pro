PRO dump_ps, draw_id, draw_cm_id, un_sensitive_id, suggest_title, GIF=gif

filename = pickfile(get_path = path, TITLE = 'select name of PS file', FILTER = '*')
;filename = './test.ps'

;widget_control, un_sensitive_id, sensitive =0 
widget_control, /hourglass


i_date = 0
i_eps  = 0
IF (not keyword_set(GIF)) then BEGIN
   tt = 'Dump PostScript Image'
   base = widget_base(title=tt, /COLUMN)   
   
   junk =  WIDGET_LABEL(base, VALUE='Plot Title:'  $
;    ,FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1' $
   )
   title_text = WIDGET_TEXT(base, /EDITABLE, UVALUE='PLOT_TITLE', $
                            VALUE = suggest_title, YSIZE=2)
   
   
   opt_but = [  'Date', 'Encapsulated PS' ]
   options = CW_BGROUP(base, /COLUMN,  /NO_REL, /NONEXCLUSIVE, /RETURN_NAME, $
                       opt_but, UVALUE='OPTIONS',    LABEL_TOP='Options:', $
;    FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1',$
   SET_VALUE = [1,0],  FRAME = 4)
   
   
   PRINT_BUTTON = WIDGET_BUTTON( BASE, $
;      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
   UVALUE='PRINT_BUTTON', $
    VALUE='PRINT')
   
   CANCEL_BUTTON = WIDGET_BUTTON( BASE, $
;      FONT='-adobe-helvetica-bold-r-normal--14-100-100-100-p-82-iso8859-1', $
   UVALUE='CANCEL_BUTTON', $
    VALUE='CANCEL')
ENDIF

WIDGET_CONTROL, draw_id, GET_VALUE=win
WSET, win   
;WINDOW, /PIXMAP, /FREE, xs = !d.x_size, ys=!d.y_size  ;Save window
;wset, !d.window
;message,/info,'copying main window into pixmap...'
;DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, win]  ;Save it
WIDGET_CONTROL, /HOURGLASS
message,/info,'tvrd() main window...'
device,retain=2
backing = tvrd(/true)
backing_num = !d.window
WIDGET_CONTROL, draw_cm_id, GET_VALUE=win2
WSET, win2
WINDOW, /PIXMAP, /FREE, xs = !d.x_size, ys=!d.y_size, retain=2  ;Save CM window
wset, !d.window
message,/info,'copying colormap window into pixmap...'
DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, win2]  ;Save it
message,/info,'tvrd() colormap window...'
device,retain=2
cm_backing = tvrd(/true)
cm_backing_num = !d.window
wset, win
 
back_s = size(backing)   
cm_back_s = size(cm_backing)   

; comb_im   = bytarr(back_s(1), back_s(2)+cm_back_s(2)+80)
; comb_im(0:(back_s(1)-1), cm_back_s(2):(cm_back_s(2)+back_s(2)-1)) = backing
; comb_im(0:(back_s(1)-1), 0:(cm_back_s(2)-1)) = cm_backing
; help,comb_im


help,backing
help,cm_backing

comb_im   = bytarr(3,back_s[2], back_s[3]+cm_back_s[3]+80)
help,comb_im
for i=0,2 do begin
    comb_im[i,0:back_s[2]-1,cm_back_s[3]:cm_back_s[3]+back_s[3]-1]=backing[i,*,*]
    comb_im[i,0:back_s[2]-1,0:cm_back_s[3]-1]=cm_backing[i,*,*]
endfor
help,comb_im


IF (not keyword_set(GIF)) then  WIDGET_CONTROL, base, /REALIZE ELSE GOTO, doit

WHILE 1 DO BEGIN                ;Internal event loop   
    ev = WIDGET_EVENT([base])
    WIDGET_CONTROL,Ev.id,GET_UVALUE=Evu

    CASE Evu OF   
        'CANCEL_BUTTON':  BEGIN   
            print, 'event for cancel button'
            goto, all_done
        ENDCASE    
        
        'PRINT_BUTTON': BEGIN   
doit:
           IF (keyword_set(GIF)) THEN BEGIN 
              tvlct,r,g,b,/get
              save,comb_im,backing,cm_backing,file='junk.dat'
              write_png,filename,comb_im
              print, 'done' 
              RETURN 
           END ELSE BEGIN
; read options
              WIDGET_CONTROL, options, get_value = buts
              WIDGET_CONTROL, title_text, get_value = title
           
            tvlct, r,g,b, /get
            set_plot, 'PS'
            !P.MULTI=0
            DEVICE, file = filename, /COLOR, BITS_PER_PIXEL=8,/PORTRAIT,$
              XSIZE=16.,YSIZE=22, $
              XOFFSET=0,YOFFSET=0
            IF (buts(i_eps) gt 0) THEN BEGIN
                print, 'writing encapsulated PS file:', filename
                DEVICE, /ENCAPSULATED
            ENDIF $
            ELSE print, 'write standard PS file:', filename
            
; switch black and white      
;            white = where(comb_im EQ 255B)
;            black = where(comb_im EQ 0B)
;            print, size(white)
;            comb_im(white(*)) = 0B
;            comb_im(black(*)) = 255B

            tvlct, r,g,b
            TV, comb_im
            PRINT, 'also write:',filename+'.gif',filename+'.jpeg'
	    XYOUTS, 50,FIX(!d.y_size*0.95), title(0), /DEVICE, CHARSIZE = 1.3, COLOR=100
            XYOUTS, 50,FIX(!d.y_size*0.92), title(1), /DEVICE, CHARSIZE = 0.7, COLOR=100
            if (buts(i_date) gt 0) THEN BEGIN
                print, 'd.x_size:', !d.y_size, systime()
                XYOUTS, 0,0, systime(), $
                  /DEVICE, CHARSIZE = 0.8, COLOR=100, ALIGNMENT=0.
            END
            DEVICE, /CLOSE
            set_plot, 'x'

	ENDELSE
all_done:
            WIDGET_CONTROL, base, /DESTROY   
;            widget_control, un_sensitive_id, sensitive = 1 
            return
        ENDCASE    
        ELSE: print, 'unspecified event', Evu
    ENDCASE   
    
ENDWHILE                        ;Event loop


END
