; $Id: mycalc.pro,v 1.1 2006/12/14 06:20:09 bwoshea Exp $

; Copyright (c) 1991-1995, Research Systems, Inc.  All rights reserved.
;   Unauthorized reproduction prohibited.
;+
; NAME:
;   Scicalc
;
; PURPOSE:
;   This routine emulates a scientific calculator.
;
; CATEGORY:
;   Widgets, math.
;
; CALLING SEQUENCE:
;   Scicalc
;
; INPUTS:
;   None.
;
; KEYWORD PARAMETERS:
;   GROUP:   The widget ID of the widget that calls Scicalc.  When this
;      ID is specified, a death of the caller results in a death of
;      Scicalc.
;
;   FONT:   A string containing the name of the X-Windows font to be
;      used for the calculator buttons and display.  If no font is
;      specified, the first available 20-point font is used.
;      On many systems, you can see the names of available fonts
;      by entering the command "xlsfonts" from the Unix command line.
;
; OUTPUTS:
;   None.
;
; OPTIONAL OUTPUT PARAMETERS:
;   None.
;
; COMMON BLOCKS:
;   SCICALCBLOCK, WTRANSBLOCK
;
; SIDE EFFECTS:
;   Initiates the XManager if it is not already running.
;
; RESTRICTIONS:
;   Math error trapping varies depending upon system.
;
; PROCEDURE:
;   Create and register the widget, allow computations, and then exit.
;
; MODIFICATION HISTORY:
;   WIDGET CALCULATOR by Keith R Crosley, RSI, October 1991
;   Created from a template written by: Steve Richards, January, 1991
;       Modified by: WSO January, 1995
;-

;------------------------------------------------------------------------------
;       procedure UpdateDisplay
;------------------------------------------------------------------------------
PRO UpdateDisplay, number

  COMMON cSciCalcBlock, cDisplayID, cDisplayString, cCurrent, cPrevious, $
            cMemory, cLastOperation

   cDisplayString = cDisplayString + number

   WIDGET_CONTROL, cDisplayID, SET_VALUE=' '

   cCurrent = FLOAT(cDisplayString)

   WIDGET_CONTROL, cDisplayID, SET_VALUE=cDisplayString

END


;------------------------------------------------------------------------------
;       procedure ResetDisplay
;
;
;------------------------------------------------------------------------------
PRO ResetDisplay, NOPREV=noprev

  COMMON cSciCalcBlock

   IF KEYWORD_SET(noprev) EQ 0 THEN $
      cPrevious = cCurrent

   cDisplayString = ''

END

;------------------------------------------------------------------------------
;       procedure SetDisplay
;------------------------------------------------------------------------------
PRO SetDisplay, string

  COMMON cSciCalcBlock

   WIDGET_CONTROL, cDisplayID, SET_VALUE=string

END


;------------------------------------------------------------------------------
;       procedure equals
;------------------------------------------------------------------------------
PRO equals, NOSET=noset

  COMMON cSciCalcBlock

   CASE cLastOperation OF

     '^': BEGIN
          cCurrent = cPrevious^cCurrent
          cPrevious = cCurrent
          IF KEYWORD_SET(noset) EQ 0 THEN ResetDisplay
          SetDisplay, STRING(cPrevious)
          cLastOperation='NOOP'
          END

     '*': BEGIN
          cCurrent = cCurrent*cPrevious
          cPrevious = cCurrent
          IF KEYWORD_SET(noset) EQ 0 THEN ResetDisplay
          SetDisplay, STRING(cPrevious)
          cLastOperation='NOOP'
          END

     '/': BEGIN
          cCurrent = cPrevious/cCurrent
          cPrevious = cCurrent
          IF KEYWORD_SET(noset) EQ 0 THEN ResetDisplay
          SetDisplay, STRING(cPrevious)
          cLastOperation='NOOP'
          END

     '-': BEGIN
          cCurrent = cPrevious-cCurrent
          cPrevious = cCurrent
          IF KEYWORD_SET(noset) EQ 0 THEN ResetDisplay
          SetDisplay, STRING(cPrevious)
          cLastOperation='NOOP'
          END

     '+': BEGIN
          cCurrent = cCurrent+cPrevious
          cPrevious = cCurrent
          IF KEYWORD_SET(noset) EQ 0 THEN ResetDisplay
          SetDisplay, STRING(cPrevious)
          cLastOperation='NOOP'
          END
     'NOOP': BEGIN
          ResetDisplay
          END

   ENDCASE
END



;------------------------------------------------------------------------------
;   procedure Scicalc_ev
;------------------------------------------------------------------------------
; This procedure processes the events being sent by the XManager.
;*** This is the event handling routine for the Scicalc widget.  It is 
;*** responsible for dealing with the widget events such as mouse clicks on
;*** buttons in the Scicalc widget.  The tool menu choice routines are 
;*** already installed.  This routine is required for the Scicalc widget to
;*** work properly with the XManager.
;------------------------------------------------------------------------------
PRO Scicalc_ev, event

  COMMON cSciCalcBlock
  COMMON wtransblock, conv1, conv2

   WIDGET_CONTROL, event.id, GET_UVALUE = eventval      ;find the user value
                        ;of the widget where
                        ;the event occured
   
     ; If a user typed an input instead of clicking the calculator's widget keypad
     ; we need to do some preprocessing
   IF eventval EQ 'DISPLAY' THEN BEGIN
        ; if a character has been typed in from the keyboard
      IF event.type EQ 0 THEN BEGIN
           ; Force "Escape" or "Clear" key to clear the last entry
         IF (event.ch EQ 27) THEN $
            eventval = 'CE' $
         ELSE $
              ; Force "Return" or "Enter" key to perform equal operation
            IF (event.ch EQ 10) THEN $
               eventval = '=' $
            ELSE $
               eventval = STRING(event.ch) 
      ENDIF
   ENDIF
   
   CASE eventval OF
   
   ;*** here is where you would add the actions for your events.  Each widget
   ;*** you add should have a unique string for its user value.  Here you add
   ;*** a case for each of your widgets that return events and take the
   ;*** appropriate action.
   
     'DISPLAY': ; ignore any typed characters that werren't trapped above
      
     '1': UpdateDisplay, '1'
   
     '2': UpdateDisplay, '2'

     '3': UpdateDisplay, '3'
   
     '4': UpdateDisplay, '4'
   
     '5': UpdateDisplay, '5'
   
     '6': UpdateDisplay, '6'
   
     '7': UpdateDisplay, '7'
   
     '8': UpdateDisplay, '8'
   
     '9': UpdateDisplay, '9'
   
     '0': UpdateDisplay, '0'
   
     '.': UpdateDisplay, '.'
   
     '=': EQUALS
   
     '^': $
         CASE cLastOperation OF
           'NOOP': BEGIN
               cLastOperation = '^'
               ResetDisplay
               END

           '^': BEGIN
               cLastOperation='^'
               cCurrent = cPrevious^cCurrent
               cPrevious = cCurrent
               ResetDisplay
               SetDisplay, STRING(cPrevious)
               END

            ELSE: BEGIN
               EQUALS, /NOSET
               cLastOperation = '^'
               ResetDisplay
               END

         ENDCASE

     '*': $
         CASE cLastOperation OF
           'NOOP': BEGIN
               cLastOperation = '*'
               ResetDisplay
               END
           '*': BEGIN
               cLastOperation='*'
               cCurrent = cCurrent*cPrevious
               cPrevious = cCurrent
               ResetDisplay
               SetDisplay, STRING(cPrevious)
               END
            ELSE: BEGIN
               EQUALS, /NOSET
               cLastOperation = '*'
               ResetDisplay
               END
         ENDCASE

     '/': $
         CASE cLastOperation OF
           'NOOP':   BEGIN
               cLastOperation = '/'
               ResetDisplay
               END
           '/': BEGIN
               cLastOperation='/'
               cCurrent = cPrevious/cCurrent
               cPrevious = cCurrent
               ResetDisplay
               SetDisplay, STRING(cPrevious)
               END
            ELSE: BEGIN
               EQUALS, /NOSET
               cLastOperation = '/'
               ResetDisplay
               END
         ENDCASE

     '-': $
         CASE cLastOperation OF
           'NOOP': BEGIN
               cLastOperation = '-' 
               ResetDisplay
               END
           '-': BEGIN
               cLastOperation='-'
               cCurrent = cPrevious-cCurrent
               cPrevious = cCurrent
               ResetDisplay
               SetDisplay, STRING(cPrevious)
               END
            ELSE: BEGIN
               EQUALS, /NOSET
               cLastOperation = '-'
               ResetDisplay
               END
         ENDCASE

     '+': $
         CASE cLastOperation OF
           'NOOP': BEGIN
               cLastOperation = '+'
               ResetDisplay
               END
           '+': BEGIN
               cLastOperation='+'
               cCurrent = cCurrent+cPrevious
               cPrevious = cCurrent
               ResetDisplay
               SetDisplay, STRING(cPrevious)
               END
            ELSE: BEGIN
               EQUALS, /NOSET
               cLastOperation = '+'
               ResetDisplay
               END
         ENDCASE

     'C': BEGIN
         cCurrent=0
         cPrevious=0
         cDisplayString=''
         cLastOperation = 'NOOP'
         ResetDisplay
         SetDisplay, STRING(cCurrent)
         END

     'CE': BEGIN
         cCurrent = 0
         cDisplayString=''
         SetDisplay, STRING(cCurrent)
         END

     '+/-': BEGIN
         cCurrent = -(cCurrent)
         SetDisplay, STRING(cCurrent)
         END

     'M': BEGIN
         cMemory = cCurrent
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'M+': BEGIN
         cMemory = cMemory+cCurrent
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'MR': BEGIN
         cCurrent = cMemory
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'MC': BEGIN
         cMemory = 0
         END

     'DEGREES': BEGIN
         conv1=!dtor
         conv2=!radeg
         END

     'RADIANS': BEGIN
         conv1 = 1
         conv2 = 1
         END

     'SQRT': BEGIN
         cCurrent = SQRT(cCurrent)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV   
         END

     'X^2': BEGIN
         cCurrent = cCurrent^2
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'X^3': BEGIN
         cCurrent = cCurrent^3
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     '1/X': BEGIN
         IF cCurrent NE 0 THEN $
            cCurrent = 1/cCurrent $ 
         ELSE cCurrent = !values.f_infinity
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'ALOG': BEGIN
         cCurrent = ALOG(cCurrent)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'ALOG10': BEGIN
         cCurrent = ALOG10(cCurrent)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'EXP': BEGIN
         cCurrent = EXP(cCurrent)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'SIN': BEGIN
         cCurrent = SIN(cCurrent*conv1)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'COS': BEGIN
         cCurrent = COS(cCurrent*conv1)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'TAN': BEGIN
         cCurrent = TAN(cCurrent*conv1)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'PI': BEGIN
         cCurrent = !pi
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'ASIN': BEGIN
         IF cCurrent GE -1 AND cCurrent LE 1 THEN $
            cCurrent = ASIN(cCurrent)*conv2 $
         ELSE cCurrent = !values.f_nan
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'ACOS': BEGIN
         IF cCurrent GE -1 AND cCurrent LE 1 THEN $
         cCurrent = ACOS(cCurrent)*conv2 $
         ELSE cCurrent = !values.f_nan
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'ATAN': BEGIN
         cCurrent = ATAN(cCurrent)*conv2
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'SINH': BEGIN
         cCurrent = SINH(cCurrent*conv1)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'COSH': BEGIN
         cCurrent = COSH(cCurrent*conv1)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     'TANH': BEGIN
         cCurrent = TANH(cCurrent*conv1)
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

     '!': BEGIN
         IF cCurrent LT 0 THEN $
            RETURN
         IF cCurrent EQ 0 THEN $
            cCurrent = 1 $
         ELSE BEGIN
            cCurrent = FIX(cCurrent)
            temp = 1.
            FOR i=1,cCurrent DO $
               temp = temp * i
            cCurrent = FLOAT(temp)
         ENDELSE
         SetDisplay, STRING(cCurrent)
         ResetDisplay, /NOPREV
         END

       ; When an event occurs in a widget that has no user value in this
       ; case statement, an error message is shown
      ELSE: MESSAGE, "Event User Value Not Found"

   ENDCASE

END ;============= end of Scicalc event handling routine task =============



;------------------------------------------------------------------------------
;   procedure Scicalc
;------------------------------------------------------------------------------
; This routine creates the widget and registers it with the XManager.
;*** This is the main routine for the Scicalc widget.  It creates the
;*** widget and then registers it with the XManager which keeps track of the 
;*** currently active widgets.
;------------------------------------------------------------------------------
PRO Mycalc, GROUP = GROUP, FONT = font

  COMMON cSciCalcBlock
  COMMON wtransblock, conv1, conv2

   ;*** If Scicalc can have multiple copies running, then delete the following
   ;*** line and the comment for it.  Often a common block is used that
   ;*** prohibits multiple copies of the widget application from running.  In
   ;*** this case, leave the following line intact.
   
   IF(XRegistered("Scicalc") NE 0) THEN RETURN      ;only one instance of
                        ;the Scicalc widget
                        ;is allowed.  If it is
                        ;already managed, do
                        ;nothing and return
   
   ;*** Next the main base is created.  You will probably want to specify either
   ;*** a ROW or COLUMN base with keywords to arrange the widget visually.
   
     ; create the main base
   wSciCalcBase = WIDGET_BASE(TITLE = "IDL Calculator", /COLUMN)
   
   ;*** Here some default controls are built in a menu.  The descriptions of
   ;*** these procedures can be found in the Scicalc_ev routine above.  If you
   ;*** would like to add other routines or remove any of these, remove them
   ;*** both below and in the Scicalc_ev routine.
   
   IF KEYWORD_SET(font) EQ 0 THEN font='*20'
   
   
   row1 = WIDGET_BASE(wSciCalcBase, /ROW)
   
   cDisplayID = WIDGET_TEXT(row1, VALUE='0                              ', $
        FONT=font, /FRAME, /ALL_EVENTS, UVALUE='DISPLAY')
   
   special = WIDGET_BASE(wSciCalcBase, /ROW, /FRAME)
   slcol   = WIDGET_BASE(special, /COLUMN)
   sl2col  = WIDGET_BASE(special, /COLUMN)
   smcol   = WIDGET_BASE(special, /COLUMN)
   sr4col   = WIDGET_BASE(special, /COLUMN, SPACE=0, YPAD=0, XPAD=0)
   sr4row   = WIDGET_BASE(sr4col, /ROW, YPAD=0, XPAD=0)
   srcol   = WIDGET_BASE(sr4row, /COLUMN)
   s4col   = WIDGET_BASE(sr4row, /COLUMN)
   
   ; COMMONLY USED FN'S:
   sqrt =    WIDGET_BUTTON(slcol, VALUE='SQRT', UVALUE='SQRT')
   x2   =  WIDGET_BUTTON(slcol, VALUE='X^2', UVALUE='X^2')
   x3   =  WIDGET_BUTTON(slcol, VALUE='X^3', UVALUE='X^3')
   x1   =  WIDGET_BUTTON(slcol, VALUE='1/X', UVALUE='1/X')
   
   
   ;LOGS:
   alog =    WIDGET_BUTTON(sl2col, VALUE='LN', UVALUE='ALOG')
   alog10 = WIDGET_BUTTON(sl2col, VALUE='LOG', UVALUE='ALOG10')
   exp  =   WIDGET_BUTTON(sl2col, VALUE='EXP', UVALUE='EXP')
   
   
   ; TRANSCENDENTALS:
   sin =   WIDGET_BUTTON(smcol, VALUE='SIN', UVALUE='SIN')
   cos =    WIDGET_BUTTON(smcol, VALUE='COS', UVALUE='COS')
   tan =   WIDGET_BUTTON(smcol, VALUE='TAN', UVALUE='TAN')
   pi  =   WIDGET_BUTTON(smcol, VALUE='PI',  UVALUE='PI')
   
   asin =   WIDGET_BUTTON(srcol, VALUE='ASIN', UVALUE='ASIN')
   acos =   WIDGET_BUTTON(srcol, VALUE='ACOS', UVALUE='ACOS')
   atan =   WIDGET_BUTTON(srcol, VALUE='ATAN', UVALUE='ATAN')
   
   sinh =   WIDGET_BUTTON(s4col, VALUE='SINH', UVALUE='SINH')
   cosh =   WIDGET_BUTTON(s4col, VALUE='COSH', UVALUE='COSH')
   tanh =   WIDGET_BUTTON(s4col, VALUE='TANH', UVALUE='TANH')
   
   ;OTHER FN'S:
   fact = WIDGET_BUTTON(sl2col, VALUE='!', UVALUE='!')
   
   ;THE DEGREE/RADIAN TOGGLE:
   togglebase = WIDGET_BASE(sr4col, /ROW, /FRAME, /EXCLUSIVE)
   
   degree = WIDGET_BUTTON(togglebase, VALUE='Degrees', $
                UVALUE = 'DEGREES', /NO_RELEASE)
   
   radian = WIDGET_BUTTON(togglebase, VALUE='Radians', $
                          UVALUE = 'RADIANS', /NO_RELEASE)
   
   
   keypad = WIDGET_BASE(wSciCalcBase, /ROW, /FRAME)
   lcol   = WIDGET_BASE(keypad, /COLUMN)
   rcol   = WIDGET_BASE(keypad, /COLUMN)
   r2col  = WIDGET_BASE(keypad, /COLUMN)
   r3col  = WIDGET_BASE(keypad, /COLUMN)
   
   row2 = WIDGET_BASE(lcol, /ROW)
   seven = WIDGET_BUTTON(row2, VALUE='7', FONT=font, UVALUE='7')
   eight = WIDGET_BUTTON(row2, VALUE='8', FONT=font, UVALUE='8')
   nine = WIDGET_BUTTON(row2, VALUE='9', FONT=font, UVALUE='9')
   
   row3 = WIDGET_BASE(lcol, /ROW)
   four = WIDGET_BUTTON(row3, VALUE='4', FONT=font, UVALUE='4')
   five = WIDGET_BUTTON(row3, VALUE='5', FONT=font, UVALUE='5')
   six = WIDGET_BUTTON(row3, VALUE='6', FONT=font, UVALUE='6')
   
   row4 = WIDGET_BASE(lcol, /ROW)
   one = WIDGET_BUTTON(row4, VALUE='1', FONT=font, UVALUE='1')
   two = WIDGET_BUTTON(row4, VALUE='2', FONT=font, UVALUE='2')
   three = WIDGET_BUTTON(row4, VALUE='3', FONT=font, UVALUE='3')
   
   row5 = WIDGET_BASE(lcol, /ROW)
   zero = WIDGET_BUTTON(row5, VALUE='0', FONT=font, UVALUE='0')
   point = WIDGET_BUTTON(row5, VALUE=' .', FONT=font, UVALUE='.')
   equals = WIDGET_BUTTON(row5, VALUE=' = ', FONT=font, UVALUE='=')
   
   clear= WIDGET_BUTTON(rcol, VALUE='C', FONT=font, UVALUE='C')
   mult = WIDGET_BUTTON(rcol, VALUE='*', FONT=font, UVALUE='*')
   div  = WIDGET_BUTTON(rcol, VALUE='/', FONT=font, UVALUE='/')
   minus= WIDGET_BUTTON(rcol, VALUE='-', FONT=font, UVALUE='-')
   plus = WIDGET_BUTTON(rcol, VALUE='+', FONT=font, UVALUE='+')
   
   ce   = WIDGET_BUTTON(r2col, VALUE='CE', FONT=font, UVALUE='CE')
   power= WIDGET_BUTTON(r2col, VALUE='^', FONT=font, UVALUE='^')
   plusmin=WIDGET_BUTTON(r2col, VALUE='+/-', FONT=font, UVALUE='+/-')
   
   cMemory  = WIDGET_BUTTON(r3col, VALUE='M', FONT=font, UVALUE='M')
   mempl= WIDGET_BUTTON(r3col, VALUE='M+', FONT=font, UVALUE='M+')
   memr = WIDGET_BUTTON(r3col, VALUE='MR', FONT=font, UVALUE='MR')
   memc = WIDGET_BUTTON(r3col, VALUE='MC', FONT=font, UVALUE='MC')
   
   
   ; INITIALIZE:
   
   cCurrent = 0.
   cPrevious = 0.
   cMemory  = 0.
   cDisplayString=''
   cLastOperation='NOOP'
   conv1 = !dtor
   conv2 = !radeg
   
     ; Create and display the widgets that are defined
   WIDGET_CONTROL, wSciCalcBase, /REALIZE
   
   WIDGET_CONTROL, degree, /SET_BUTTON

     ; register the widgets with the XManager and pass through the group leader
   XManager, "Scicalc", wSciCalcBase, $
         EVENT_HANDLER = "Scicalc_ev", $
         GROUP_LEADER = GROUP

END ;==================== end of Scicalc main routine =======================









