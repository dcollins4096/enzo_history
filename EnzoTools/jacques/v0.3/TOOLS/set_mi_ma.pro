PRO set_mi_ma, group, wid_id, MI_MA=mi_ma
; provides little widget to set minimum & max. values for
; colormap
; Tom Abel: 2/00
   on_error, 2                  ;Return to caller if an error occurs.
; set defaults if noe value were provided
   if (n_elements(mi_ma) eq 0) then mi_ma = [0.,1.]
   if (n_elements(fixed) eq 0) then fixed = 0 ; if undefined let it float
   state =  {mi_ma:mi_ma, fixed:fixed}

; create my own base
   base =  widget_base(title="Min and Max", /COLUMN, UVALUE=state, $
                      PRO_SET_VALUE='set_mi_ma_values', GROUP_LEADER=group)
   wid_id =  base

   min_value = cw_field(base, $
                        uvalue = 'MINIMUM', $
                        title='Min', $
                        value=mi_ma(0), /floating, /all_events)

   MAX_VALUE = CW_FIELD(BASE,  $
                        UVALUE = 'MAXIMUM', $
                        TITLE='Max', $
                        VALUE=mi_ma(1), /FLOATING, /ALL_EVENTS)
   opt_but = ['dynamic','static']
   options = CW_BGROUP(base, /ROW,   /EXCLUSIVE, /RETURN_NAME, $
                       opt_but, UVALUE='OPTIONS', $
                       SET_VALUE = [state.fixed])
   DONE_BUTTON = WIDGET_BUTTON( BASE, $
                                UVALUE='DONE', $
                                VALUE='DONE')
   
   WIDGET_CONTROL, base, /REALIZE
   XMANAGER, 'set_mi_ma', BASE, EVENT_HANDLER='SET_MI_MA_EVENT', /NO_BLOCK
   
   RETURN
END

PRO set_mi_ma_values, id, mi_ma

RETURN
END

PRO set_mi_ma_event, event
   
   WIDGET_CONTROL, event.top, GET_UVALUE = state
   fixed =  state.fixed
   min   =  state.mi_ma(0)
   max   =  state.mi_ma(1)

   WIDGET_CONTROL, event.id, GET_UVALUE = eventval
   CASE eventval OF
      "MINIMUM": BEGIN
         widget_control, event.id, GET_VALUE = min
      END
      "MAXIMUM": BEGIN
         widget_control, event.id, GET_VALUE = max
      END
      "OPTIONS": BEGIN
         widget_control, event.id, GET_VALUE = fixed
      END
      ELSE:
   ENDCASE

   state.fixed    = fixed 
   state.mi_ma[0] = min   
   state.mi_ma[1] = max   

;   print,'state.mi_ma: ',state.mi_ma

   WIDGET_CONTROL, event.top, SET_UVALUE = state

   IF eventval eq "DONE" THEN WIDGET_CONTROL, event.top, /DESTROY
  
RETURN
END         

