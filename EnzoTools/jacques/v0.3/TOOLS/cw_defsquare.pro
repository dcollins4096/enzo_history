pro CW_DEFSQUARE_DRAW, s, i, FILL = fill
; Draw the outline of square
; of the region 
; Use the XOR drawing mode.

n = s.npts
if n lt 1 then return
   
WSET, s.win   
DEVICE, SET_GRAPHICS=6          ;Xor drawing mode   
;col = 1
;while col lt !d.table_size do col = col + col

; mqk
col=255

WIDGET_CONTROL, s.xy_pnts, GET_UVALUE=xy, /NO_COPY  ;Get ROI
  xsave = !x.s & ysave = !y.s       ;Set scaling to pixel coords
  p = float([!d.x_size, !d.y_size])
  f = s.offset / p
  q = s.zoom / p
  !x.s = [f(0), q(0)]
  !y.s = [f(1), q(1)]

  if i lt 0 then plots, COLOR=col, xy(*, 0:n-1)+.5 $ ;All of it?
  else begin 
      plots, COLOR=col, xy(*, i:i+1)+.5 ;One segment
      print, 'xy:', xy
  END            
  IF KEYWORD_SET(FILL) then POLYFILL, xy(*,0:n-1), COLOR=col

!x.s = xsave & !y.s = ysave
WIDGET_CONTROL, s.xy_pnts, SET_UVALUE=xy, /NO_COPY  ;Set ROI   
DEVICE, SET_GRAPHICS=3          ;Copy mode   
end

PRO CW_DEFSQUARE_event, ev, s,  xy_pixels, side_length
; This routine is only called from the CW_DEFSQUARE event loop.
; ev = event structure, s = state structure.

s.button = s.button or ev.press xor ev.release  ;New button state   
n = s.npts
x = (ev.x - s.offset(0)) / s.zoom(0)    ;Pixel coordinates
y = (ev.y - s.offset(1)) / s.zoom(1)   


if (x lt 0) or (y lt 0) or $            ;Within region?
    (x ge s.image_size(0)) or (y ge s.image_size(1)) then return
if ev.press ne 0 then s.drag = [x,y]    ;Start of drag operation

if s.button ne 0 then begin     ;Drag
    if n gt 0 then CW_DEFSQUARE_draw, s, -1 ;Remove old
    t = s.drag
    n = 5
;    th = ([x,y]-t)/abs([x,y]-t) * min(abs([x,y]-t))
;    x = t(0)+th(0)
;    y = t(1)+th(1)
;    xy = [[t], [x, t(1)], [x, y], [t(0), y], [t]]
    th = ([x,y]-t)/abs([x,y]-t) * min(abs([x,y]-t))
    hs = min(th)
    x = t(0)-hs
    y = t(1)-hs
    xy = [[x, y], [x+2*hs,y], [x+2*hs, y+2*hs], [x, y+2*hs], [x,y]]
; update text
WIDGET_CONTROL, s.pos_w, $
;    SET_VALUE=string(x, y0, format='("Position: ",i,", ",i)')
    SET_VALUE=string(side_length*2*hs/xy_pixels, $
                     format='("side_length: ",f)')

    WIDGET_CONTROL, s.xy_pnts, SET_UVALUE=xy, /NO_COPY ;Restore UVALUE
    s.npts = n
    CW_DEFSQUARE_draw, s, -1
ENDIF                           ;DRAG
return
   
done0: WIDGET_CONTROL, s.xy_pnts, SET_UVALUE=xy, /NO_COPY   
end   


function CW_DEFSQUARE, draw, xy_pixels, side_length, $
           ZOOM = zoom, IMAGE_SIZE = image_size, $   
           OFFSET = offset, RESTORE = restore, ORDER = order

base = widget_base(title='Region of Interest', /COLUMN)   
xy_pnts = WIDGET_TEXT(base, YSIZE=6, /FRAME, UVALUE=0, $
    value=['Open a square on the drawing area:', $
           'Add with left button: drag', $   
           '  right  button: Cancel',   $
           '  middle button: Accept & Return', $
          'Note that you have to hit draw,',  $
          'to actually zoom in.' ])   
Options = CW_BGROUP(base, /ROW, /NO_RELEASE, /RETURN_NAME,  $
    ['Clear All', 'Cancel'])   
junk = CW_BGROUP(base, /ROW, /EXCLUSIVE, /NO_REL, /RETURN_NAME, $
    ['Square'], SET_VALUE=0)
junk = CW_BGROUP(base, /ROW, /NO_REL, /RETURN_NAME, ['Done'])
pos_w = WIDGET_TEXT(base, YSIZE=1, XSIZE=18, /FRAME, $
    VALUE='side length:    0.')

WIDGET_CONTROL, draw, GET_VALUE=win
WSET, win   
   
if n_elements(zoom) le 0 then zoom = [1,1]
if n_elements(image_size) le 0 then image_size = [!d.x_size, !d.y_size] / zoom
if n_elements(offset) le 0 then offset = [0,0]   
p  = offset + image_size /2   
if (!version.os NE 'MacOS') THEN TVCRS, p(0), p(1), /DEVICE ELSE TVCRS, 1
   
WINDOW, /PIXMAP, /FREE, xs = !d.x_size, ys=!d.y_size  ;Save window
backing = !d.window
DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, win]  ;Save it

s = { CW_DEFSQUARE_STRUCT, $       ;Structure containing state
    base: base, $       ;Main base widget
    xy_pnts: xy_pnts, $ ;Current roi vertex list
    npts : 0L, $        ;# of points in current roi
    subs : pos_w, $     ;Widget holding prev subscripts
    pos_w : pos_w, $    ;Position text widget
    mode: 0, $          ;major mode
    amode: 0, $         ;0 for add, 1 for remove
    draw: draw, $       ;draw widget id
    win:  win, $        ;draw widget window #
    button: 0, $        ;button state
    image_size : long(image_size), $   ;Image array size
    backing: backing, $ ;Pixmap for backing store
    offset: fix(offset), $   ;offset of array within window
    zoom : fix(zoom), $ ;zoom factor
    order : KEYWORD_SET(order), $  ;Image order
    drag: [0,0]}        ;Beginning of drag motion
   
WIDGET_CONTROL, base, /REALIZE
;WSHOW, win
   
WHILE 1 DO BEGIN                ;Internal event loop   

    ev = WIDGET_EVENT([base, draw])
    IF ev.id eq draw THEN BEGIN
        IF N_elements(ev.press) gt 0 THEN BEGIN 
            if ev.press eq 2 then goto, all_done $
            else $
              if ev.press eq 4 then BEGIN 
                xy = -1 
                goto, get_me_out
            ENDIF
        ENDIF
    ENDIF
    n = s.npts
    if ev.id eq draw then CW_DEFSQUARE_EVENT, ev, s, xy_pixels, side_length $
    else case ev.value of   
'Clear All':  BEGIN
    WIDGET_CONTROL, s.subs, GET_UVALUE=t, /NO_COPY  ;Clr list of subscripts
    t = 0
    WSET, win
    DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, backing]  ;Restore it
    s.npts = 0
    ENDCASE

'Cancel':  BEGIN   
    xy = -1
    goto, get_me_out
    ENDCASE    

'Square' :  

'Done': BEGIN  
all_done: 
    WIDGET_CONTROL, s.subs, GET_UVALUE=t, /NO_COPY  ;List of subscripts
    if n_elements(t) gt 0 then xy(t) = 1
    print, 's_npts:',s.npts
    print, s.xy_pnts
    WIDGET_CONTROL, s.xy_pnts, GET_UVALUE=xy, /NO_COPY ;Get ROI
    IF N_ELEMENTS(xy) gt 1 THEN $
      if (abs(xy(0,0)-xy(0,1)) le 2) THEN xy = -1 ; if too small don't return
get_me_out:
    IF KEYWORD_SET(restore) then begin ;Undo damage?
        WSET, win
        DEVICE, copy = [0,0, !d.x_size, !d.y_size, 0, 0, backing] ;Restore it
    ENDIF
    DEVICE, SET_GRAPHICS=3
    WDELETE, backing
    WIDGET_CONTROL, base, /DESTROY   
    return, xy
    ENDCASE    
ENDCASE   
ENDWHILE            ;Event loop
END   



