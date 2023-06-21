; Copyright (c) 1995, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.
;+
; NAME: Editor
;
; PURPOSE: Display an ASCII text file using widgets.
;
; MAJOR TOPICS: Text manipulation.
;
; CALLING SEQUENCE: Editor [, filename]
;
; INPUTS:
;     Filename:	A scalar string that contains the filename of the file
;		to display.  The filename can include a path to that file.
;               If the filename is omitted, the user will be prompted for
;               a filename, via the system function, PICKFILE().
;
; KEYWORD PARAMETERS:
;	FONT:   The name of the font to use.  If omitted use the default
;		font.
;
;	HEIGHT:	The number of text lines that the widget should display at one
;		time.  If this keyword is not specified, 24 lines is the 
;		default.
;
;	WIDTH:	The number of characters wide the widget should be.  If this
;		keyword is not specified, 80 characters is the default.
;
; PROCEDURE: Editor reads, writes and manipulates text strings ...
;
; MAJOR FUNCTIONS and PROCEDURES:
;
; COMMON BLOCKS and STRUCTURES:
;
; SIDE EFFECTS:
;	Triggers the XMANAGER if it is not already in use.
;
; MODIFICATION HISTORY:  Written by:  WSO, RSI, January 1995
;                        Modified by: ...
;-


FUNCTION FindDelimiter, text, delimiters, column, row, delimiterIndex, $
             endOfLineMarks, regionMarks, START_OF_LINE=startOfLine, $
             END_OF_LINE=endOfLine

   IF (KEYWORD_SET(startOfLine) NE 0) AND (column EQ 0) THEN $
      RETURN, 1

   IF (N_ELEMENTS(endOfLineMarks) EQ 0) THEN $
      endOfLineMarks = ""

   IF (N_ELEMENTS(regionMarks) EQ 0) THEN $
      regionMarks = ""

   delimiters = [delimiters, endOfLineMarks, regionMarks]

   WHILE 1 DO BEGIN

      lastColumn = STRLEN(text(row))-1

        ; Find the delimiter
      WHILE (column LE lastColumn) DO BEGIN

         delimiterIndex = STRPOS(delimiters, STRMID(text(row), column, 1))

         delimiterGroup = (WHERE (delimiterIndex NE -1))(0)

         CASE delimiterGroup OF

              ; If delimiter found is a region marker (quote, double quote, etc)
            2: BEGIN ; Skip to next region marker
               currentColumn = column + 1
               column = STRPOS(STRMID(text(row), currentColumn, lastColumn-column), $
                       STRMID(delimiters(delimiterGroup), delimiterIndex(delimiterGroup), 1))
               IF column EQ -1 THEN $
                  column = lastColumn + 1 $
               ELSE $
                  column = column + currentColumn + 1
            ENDCASE

              ; end of line mark found - force completion of for loop
            1: column = lastColumn+1

              ; normal delimiter found
            0: RETURN, 1

              ; Go to next column
            ELSE: $
               column = column + 1

         ENDCASE

      ENDWHILE

        ; We've reached the end of the current line

        ; If the keyword is set then return success
      IF (KEYWORD_SET(endOfLine) NE 0) THEN $
         RETURN, 1

        ; Go to the beginning of the next row of text
      column = 0
      row = row + 1

        ; Are we past the end of the last line of the text
      IF ((SIZE(text))(0) EQ 0) OR (row GE (SIZE(text))(1)) THEN $
         RETURN, 0

      IF (KEYWORD_SET(startOfLine) NE 0) THEN $
         RETURN, 1

   ENDWHILE

     ; If the delimiters were not found in the string
   RETURN, 0

END


PRO GetNextWord, wTextEdit, SKIP_NUMBERS=skipNumbers

     ; Define a list of characters that make words in the text
     ; Notice single quotes to avoid conflict with constant designaters
   symbols = '0123456789'
   alphabet = "_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

   IF KEYWORD_SET(skipNumbers) NE 0 THEN $
      wordCharacters = alphabet $
   ELSE $
      wordCharacters = alphabet + symbols

     ; Define a list of delimiters that separates words in the text
   wordDelimiters = " (){}/@#&!?*%^~\|-+=<>.[]:,	"

     ; Define the end of line if it's not the carriage return
   endOfLineMarks = ";$"

     ; Define region delimiters to ignore - in our case strings
   regionMarks = "'" + '"'

     ; Define all delimiters that could spearate words
   allDelimiters = [wordDelimiters, endOfLineMarks, regionMarks]

     ; Get the text from the text widget
   WIDGET_CONTROL, wTextEdit, GET_VALUE=text

     ; Get the cursor position as a character offset from start of text
   cursorOffsetPos = WIDGET_INFO(wTextEdit, /TEXT_SELECT)

     ; This is the position to start the word search
   startPos = cursorOffsetPos(0) + cursorOffsetPos(1)

     ; initialize wordlength to zero length - in case no word was found
   wordLength = 0

     ; Convert the character offset to a column/row position
   columnRow = WIDGET_INFO(wTextEdit, TEXT_OFFSET_TO_XY=startPos)

   row = columnRow(1)
   column = columnRow(0)

        ; Find the first delimiter
   found = FindDelimiter(text, wordDelimiters, column, row, delimiterIndex, $
             endOfLineMarks, regionMarks, /START_OF_LINE)

   IF found THEN BEGIN
        ; Find the first alphanumeric character
      found = FindDelimiter(text, wordCharacters, column, row, delimiterIndex, $
             endOfLineMarks, regionMarks)

      IF found THEN BEGIN

         firstCharPos = [column, row]

           ; Find the second delimiter
         found = FindDelimiter(text, wordDelimiters, column, row, delimiterIndex, $
             endOfLineMarks, regionMarks, /END_OF_LINE)

         wordLength = column - firstCharPos(0)

         startPos = WIDGET_INFO(wTextEdit, TEXT_XY_TO_OFFSET=firstCharPos)

      ENDIF

   ENDIF

   WIDGET_CONTROL, wTextEdit, SET_TEXT_SELECT=[startPos, wordLength]

   WIDGET_CONTROL, wTextEdit, /INPUT_FOCUS

END


PRO FindText, state, token

     ; Get the text from the text widget
   WIDGET_CONTROL, state.wTextEdit, GET_VALUE=text

     ; First get the current cursor position and see if any tokens were
     ; found after that position.  If not circle back to the beginning
     ; of the text and continue searching from there.

     ; Get the cursor position as a character offset from start of text
   cursorOffsetPos = WIDGET_INFO(state.wTextEdit, /TEXT_SELECT)

     ; Convert the character offset to a [column,row] position
   cursorColumnRow = WIDGET_INFO(state.wTextEdit, $
                 TEXT_OFFSET_TO_XY=cursorOffsetPos(0)+cursorOffsetPos(1))

     ; For readability only convert the cursorColumnRow array to a structure
   current = {column:cursorColumnRow(0), row:cursorColumnRow(1)}

     ; Search for a token in the current row after the current column
   column = STRPOS(text(current.row), token, current.column)

     ; If any tokens were found on the current line - 
   IF column NE -1 THEN $
        ; Save this [column,row] position as the location where the token was found
      location = [column, current.row] $

   ELSE BEGIN

        ; The token wasn't found in the current row - therefore
        ; search for the token from the next row to the last row
      nextRow = current.row+1
      lastRow = (SIZE(text))(1)-1

        ; If at the end of the text widget
      IF nextRow GT lastRow THEN $
         location = -1 $
      ELSE $
         location = STRPOS(text(nextRow:lastRow), token)

        ; row is equal to the first row indice where a token was found
      row = (WHERE(location GT -1))(0)

        ; If any tokens were found - 
      IF row NE -1 THEN $
           ; Save the token location [column, row]
         location = [location(row), row+current.row+1] $

      ELSE BEGIN

           ; The token wasn't found in the current row to the last row - therefore
           ; search for the token from the first row to the current row
         location = STRPOS(text(0:current.row), token)

           ; row is equal to the first row indice where a token was found
         row = (WHERE(location GT -1))(0)

           ; If any tokens were found - 
         IF row NE -1 THEN $
              ; Save the token location [column, row]
            location = [location(row), row] $
         ELSE $
            location = -1

      ENDELSE
   ENDELSE
 
     ; If location(column) is valid - the token was found
   IF location(0) GT -1 THEN BEGIN

        ; Get the character offset from the column, row position
      startPos = WIDGET_INFO(state.wTextEdit, TEXT_XY_TO_OFFSET=location)

        ; Select the text that was found to match the search token
      WIDGET_CONTROL, state.wTextEdit, SET_TEXT_SELECT=[startPos, STRLEN(token)]

        ; Enable the any controls that require text to be selected
      TextSelected, state, 1

        ; Give the keyboard input focus back to the text widget
      WIDGET_CONTROL, state.wTextEdit, /INPUT_FOCUS

   ENDIF ELSE $
        ; Notify the user that the token was not found
      button = WIDGET_MESSAGE(/INFO, 'The text "' + token + '" was not found.')

END


FUNCTION CountWords, wTextEdit

   wordCount = 0

   GetNextWord, wTextEdit, /SKIP_NUMBERS

     ; Get selection to see if word was found 
   cursorOffsetPos = WIDGET_INFO(wTextEdit, /TEXT_SELECT)

     ; While the word length is greater than zero - continue
   WHILE cursorOffsetPos(1) GT 0 DO BEGIN

      wordCount = wordCount + 1

      GetNextWord, wTextEdit, /SKIP_NUMBERS

        ; Get selection to see if word was found 
      cursorOffsetPos = WIDGET_INFO(wTextEdit, /TEXT_SELECT)

   ENDWHILE

   RETURN, wordCount

END


PRO DisplayEditorInfo, parent

   infoText = [ $
     "   Editor is an IDL example of a simple text "+ $
     "editing application.  It demonstrates the ability to read and "+ $
     "write text from a file and display it in a text widget.  It "+ $
     "also provides a text search capability. ", "", $
     '   Files can be opened and saved via commands in the "File" '+ $
     "menu on the window's menubar.  The New command allows you "+ $
     "to create a new text file. The Open menu command displays a "+ $
     "standard file dialog to help you select a file to edit.  The "+ $
     "Save and Save As menu commands write the text files back out "+ $
     "to the disk.  The Exit Editor command will cause the Editor "+ $
     "example to close and return to IDL.", "", $
     "   To search for text in the text editor, enter the desired "+ $
     "text into the Find Text field and click the Find button."]

   ShowInfo, TITLE="Editor Example Information", GROUP=parent, INFOTEXT=infoText, $
     HEIGHT=15, WIDTH=80

END


PRO EditorKilled, widgetID

   WIDGET_CONTROL, GET_UVALUE = state, widgetID, /NO_COPY

   ; If the Editor could be closed via the window manager, you would
   ; need to add a request here for the user to save any changes they
   ; might have made.  Currently the only way to close this editor 
   ; is through the "Exit Editor" menu command.

   ; No need to reset the user value since it will no longer be used
END


PRO RequestToSave, state

     ; If the user editted (changed) the text in the text editor...
   IF state.textEditted NE 0 THEN BEGIN

        ; See if they want to save any changes they made
      buttonPushed = WIDGET_MESSAGE(/QUESTION, $
                       ["Save changes to ", '"'+state.fileName+'"'])

        ; If the responded by clicking the "Yes" button ...
      IF buttonPushed EQ "Yes" THEN BEGIN

           ; Save the text properly
         IF state.fileName EQ "Untitled" THEN $
            SaveAsText, state $
         ELSE $
            SaveText, state

      ENDIF
   ENDIF
END


PRO GetText, state, FILENAME=fileName

     ; If any current text changed in the text editor -
     ; request that the user save their changes first
   RequestToSave, state

     ; If there isn't a filename defined - request one from the user
   IF (NOT KEYWORD_SET(fileName)) THEN $
      fileName = PICKFILE(GROUP=state.wEditorWindow)

     ; If the user didn't cancel the standard file dialog or
     ; the filename already existed
   IF STRLEN(fileName) GT 0 THEN BEGIN

        ; Save it in the state structure
      state.fileName = fileName

        ; Set the title of the window to include the new file name
      WIDGET_CONTROL, state.wEditorWindow, $
        TLB_SET_TITLE="Editor - " + state.fileName

        ; Display the wait cursor
      WIDGET_CONTROL, state.wEditorWindow, /HOURGLASS

        ; Open the file to read
      OPENR, unit, state.filename, /GET_LUN, ERROR=error

        ; If an error occurred
      IF error LT 0 THEN $
           ; Notify the User that an error occurred
         buttonPushed = WIDGET_MESSAGE( [!err_string, $
                          ' Cannot open the file ' + state.filename] ) $
      ELSE BEGIN
         maxLines = 1000
         lineIncrement = 250
         text = STRARR(maxLines)	; Maximum # of lines
         lineIndex = 0
         lineOfText = ""
         WHILE NOT EOF(unit) DO BEGIN
            READF, unit, lineOfText
            text(lineIndex) = lineOfText
            lineIndex = lineIndex + 1
              ; If the maximum number of lines is hit,
              ; increase the array size by lineIncrement
            IF (lineIndex EQ maxLines) THEN BEGIN
               text = [text, STRARR(lineIncrement)]
               maxLines = maxLines + lineIncrement
            ENDIF

         ENDWHILE

         text = text(0:(lineIndex-1)>0)

         FREE_LUN, unit            ; Free and close the file unit.

           ; Initialize the flag to show the text is initially unchanged
         state.textEditted = 0

           ; Insert the text into the text widget
         WIDGET_CONTROL, state.wTextEdit, SET_VALUE=text

      ENDELSE

   ENDIF

END


PRO SaveAsText, state

     ; Request a filename by which the file is to be saved
   fileName = PICKFILE(/WRITE, GROUP=state.wEditorWindow)

     ; If the user didn't cancel the standard file dialog
   IF STRLEN(fileName) GT 0 THEN BEGIN

        ; Update the window title to reflect the new filename
      WIDGET_CONTROL, state.wEditorWindow, TLB_SET_TITLE="IDL Editor - " + fileName

        ; Store the filename in the application's state structure
      state.fileName = fileName

        ; Write the text to disk
      SaveText, state

   ENDIF
END


PRO SaveText, state

     ; Display the wait cursor
   WIDGET_CONTROL, state.wEditorWindow, /HOURGLASS

     ; Open the file for writing
   OPENW, unit, state.filename, /GET_LUN, ERROR=error

     ; If an error occurred - display an error message
   IF error LT 0 THEN $
      buttonPushed = WIDGET_MESSAGE([!err_string, $
                       'Can not display the file:', state.filename] ) $
   ELSE BEGIN
        ; Get the text to save from the text widget
      WIDGET_CONTROL, state.wTextEdit, GET_VALUE=text

        ; Count how many lines of text to save
      lineCount = (SIZE(text))(1)

      FOR lineIndex = 0, lineCount-1 DO $

           ; Write each line to the file
         PRINTF, unit, text(lineIndex)

      FREE_LUN, unit            ; Free and close the file unit.

        ; Reset the flag to show the text is saved
      state.textEditted = 0

   ENDELSE

END


PRO TextSelected, state, enabled

     ; Update these widgets when text is selected or deselected
   WIDGET_CONTROL, state.wUppercase, SENSITIVE=enabled
   WIDGET_CONTROL, state.wLowercase, SENSITIVE=enabled

END


PRO EditorEventHdlr, event

     ; Get the state structure stored in the user value of the window
   WIDGET_CONTROL, GET_UVALUE = state, event.top, /NO_COPY

     ; Determine in which widget the event occurred
   CASE event.id OF
      
      state.wEditorWindow: BEGIN ; The window has been sized

           ; Get the new size of the window
         WIDGET_CONTROL, state.wEditorWindow, TLB_GET_SIZE=windowSize
        
           ; Determine the change in the window size
         deltaX = windowSize(0) - state.windowSize(0)
         deltaY = windowSize(1) - state.windowSize(1)

           ; Store the new size in the state structure for later comparisons
         state.windowSize = windowSize

           ; Get the pixel size of the text widget
         textEditGeometry = WIDGET_INFO(state.wTextEdit, /GEOMETRY)

           ; Determine the new size based on the amount the window grew
         newTextEditXSize = textEditGeometry.scr_xsize + deltaX
         newTextEditYSize = textEditGeometry.scr_ysize + deltaY

           ; Resize the text widget accordingly
         WIDGET_CONTROL, state.wTextEdit, SCR_XSIZE=newTextEditXSize, $
           SCR_YSIZE=newTextEditYSize

      ENDCASE

      state.wNewButton: BEGIN

           ; If any text changed in the text editor -
           ; request that the user save the changes first
         RequestToSave, state

           ; Set the default filename to "Untitled"
         state.fileName = "Untitled"

           ; Initialize the flag to show the text is initially unchanged
         state.textEditted = 0

           ; Clear any text currently in the text widget
         WIDGET_CONTROL, state.wTextEdit, SET_VALUE=""

           ; Set the title of the window
         WIDGET_CONTROL, event.top, $
           TLB_SET_TITLE="IDL Editor - " + state.fileName

      ENDCASE
         
      state.wOpenButton: $

           ; Get the text from the file and insert it into the text widget
         GetText, state
         
      state.wSaveButton: BEGIN

           ; If a new file was created - first request a filename and save
         IF state.fileName EQ "Untitled" THEN $
            SaveAsText, state $
         ELSE $
              ; Save the text to disk
            SaveText, state

      ENDCASE

      state.wSaveAsButton: $
           ; Request a filename and save the file to disk
         SaveAsText, state
         
      state.wExitButton: BEGIN
           ; If any text changed in the text editor -
           ; request that the user save the changes first
         RequestToSave, state

           ; Restore the state value before the widget app is destroyed
           ; so the KILL_NOTIFY procedure can still use it
         WIDGET_CONTROL, SET_UVALUE = state, event.top, /NO_COPY

           ; Exit the IDL Editor widget application
         WIDGET_CONTROL, event.top, /DESTROY

         RETURN

      ENDCASE
         
      state.wTextEdit: BEGIN

           ; If text has been selected enable the appropriate widgets
           ; initially assume that no text is selected
         selected = 0

           ; If the user selected some text
         IF event.type EQ 3 THEN BEGIN

            IF event.length GT 0 THEN $
               selected = 1  ; Enable the upper and lower case pushbuttons

         ENDIF ELSE $
              ; The user entered or deleted text in the text widget
              ; therefore flag the text as changed to force SAVE request
            state.textEditted = 1

         TextSelected, state, selected

      ENDCASE

      state.wFindTextID: BEGIN

           ; Get the text to search for from the Find Text widget
         WIDGET_CONTROL, event.id, GET_VALUE=findText

           ; If the user has entered text into the Find Text field -
         IF STRLEN(findText(0)) GT 0 THEN $
            WIDGET_CONTROL, state.wFind, /SENSITIVE $  ; enable the Find button
         ELSE $
            WIDGET_CONTROL, state.wFind, SENSITIVE=0  ; disable the Find button

      ENDCASE

      state.wFind: BEGIN

           ; Get the text to search for from the Find Text widget
         WIDGET_CONTROL, state.wFindTextID, GET_VALUE=findText

           ; If there is at least one character ...
         IF STRLEN(findText(0)) GT 0 THEN $
              ; Search for the text in the editor
            FindText, state, findText(0)

      ENDCASE
         
      state.wUppercase: BEGIN

           ; Get the selected text in the editor
         WIDGET_CONTROL, state.wTextEdit, GET_VALUE=selection, /USE_TEXT_SELECT

           ; If there is at least one character selected...
         IF STRLEN(selection(0)) GT 0 THEN BEGIN

              ; Convert the selection to uppercase 
            selection = STRUPCASE(selection)

              ; Set the selection to the uppercase conversion
            WIDGET_CONTROL, state.wTextEdit, SET_VALUE=selection, $
              /USE_TEXT_SELECT, /NO_NEWLINE

         ENDIF
      ENDCASE
         
      state.wLowercase: BEGIN

           ; Get the selected text in the editor
         WIDGET_CONTROL, state.wTextEdit, GET_VALUE=selection, /USE_TEXT_SELECT

           ; If there is at least one character selected...
         IF STRLEN(selection(0)) GT 0 THEN BEGIN

              ; Convert the selection to lowercase 
            selection = STRLOWCASE(selection)

              ; Set the selection to the lowercase conversion
            WIDGET_CONTROL, state.wTextEdit, SET_VALUE=selection, $
              /USE_TEXT_SELECT, /NO_NEWLINE

         ENDIF
      ENDCASE
         
      state.wNextWord: BEGIN
         GetNextWord, state.wTextEdit, /SKIP_NUMBERS

           ; Get selection to see if word was found 
         cursorOffsetPos = WIDGET_INFO(state.wTextEdit, /TEXT_SELECT)

         enabled = (cursorOffsetPos(1) GT 0)

         WIDGET_CONTROL, state.wUppercase, SENSITIVE=enabled
         WIDGET_CONTROL, state.wLowercase, SENSITIVE=enabled

      ENDCASE

      state.wCountWords: BEGIN
         start = SYSTIME(1)
         wordCount = CountWords( state.wTextEdit )
         displayText = 'The number of words in "' + state.fileName + $
                        '" is ' + STRING(wordCount, FORMAT='(I0)')
         buttonPushed = WIDGET_MESSAGE(/INFORMATION, displayText)
      ENDCASE

      state.wInfo: $ ; Display information about the IDL Editor
         DisplayEditorInfo, event.top

      ELSE: $ ; We erroneously received an event for a widget we weren't expecting
         buttonPushed = WIDGET_MESSAGE("An event occurred for a non-existent widget")

      ENDCASE

     ; Reset the windows user value to the updated state structure
   WIDGET_CONTROL, SET_UVALUE = state, event.top, /NO_COPY

END


PRO MyEditor, fileName,  $
              NEW=NEW,WIDTH = WIDTH, HEIGHT = HEIGHT, FONT = font, GROUP = group

     ; If keywords not set - set to defaults
   IF (NOT(KEYWORD_SET(height))) THEN $
      height = 24
   IF(NOT(KEYWORD_SET(width))) THEN $
      width = 80

     ; Create the text editor window with a menu bar
   wEditorWindow = WIDGET_BASE(TITLE = "Editor - Untitled", $
                   MBAR=menuBar, /COLUMN, /TLB_SIZE_EVENTS, $
                   TLB_FRAME_ATTR=8)

     ; Build the menubar
   wFileMenu = WIDGET_BUTTON(menuBar, VALUE="File", /MENU)

     ; Add the following menu items to the window's File menu
   wNewButton = WIDGET_BUTTON(wFileMenu, VALUE="New")
   wOpenButton = WIDGET_BUTTON(wFileMenu, VALUE="Open...")
   wSaveButton = WIDGET_BUTTON(wFileMenu, /SEPARATOR, VALUE="Save")
   wSaveAsButton = WIDGET_BUTTON(wFileMenu, VALUE="Save As...")
   wExitButton = WIDGET_BUTTON(wFileMenu, /SEPARATOR, VALUE="Exit Editor")
   
     ; If the FONT keyword was defined with a font
   IF N_ELEMENTS(font) GT 0 THEN $
        ; Create the text editor's text widget with the specified font
      wTextEdit = WIDGET_TEXT(wEditorWindow, XSIZE=WIDTH, YSIZE=HEIGHT, $
                   /SCROLL, FONT=font, /EDITABLE, /ALL_EVENTS) $
   ELSE $
        ; Create the text editor's text widget with the default font
      wTextEdit = WIDGET_TEXT(wEditorWindow, XSIZE=WIDTH, YSIZE=HEIGHT, $
                   /SCROLL, /EDITABLE, /ALL_EVENTS)

     ; Create a base for all the edit controls (button, etc.)
   wControlBase = WIDGET_BASE(wEditorWindow, /COLUMN)
   
   wTopBase = WIDGET_BASE(wControlBase, SPACE=20, /ROW)
   
   wFindBase = WIDGET_BASE(wTopBase, /ROW)
   
     ; Create a label for the search text field
   wFindLabel = WIDGET_LABEL(wFindBase, VALUE="Find Text:")

     ; Create the search text field
   wFindTextID = WIDGET_TEXT(wFindBase, XSIZE=24, YSIZE=1, $
                   /EDITABLE, /ALL_EVENTS)
   
     ; Create the search text pushbutton
   wFind = WIDGET_BUTTON(wFindBase, VALUE = "Find")

     ; Create a pushbutton to display information about this application
   wInfo = WIDGET_BUTTON(wTopBase, VALUE = "Info...")


   wBottomBase = WIDGET_BASE(wControlBase, /ROW, SPACE=60)
   
   wCaseBase = WIDGET_BASE(wBottomBase, /ROW)
   
   wUppercase = WIDGET_BUTTON(wCaseBase, VALUE = "Uppercase")
   
   wLowercase = WIDGET_BUTTON(wCaseBase, VALUE = "Lowercase")

   wWordBase = WIDGET_BASE(wBottomBase, /ROW)
   
   wNextWord = WIDGET_BUTTON(wWordBase, VALUE = "Next Word")

   wCountWords = WIDGET_BUTTON(wWordBase, VALUE = "Count Words")

     ; Initially disable the pushbuttons, until the required user event
   WIDGET_CONTROL, wUppercase, SENSITIVE=0
   WIDGET_CONTROL, wLowercase, SENSITIVE=0
   WIDGET_CONTROL, wFind, SENSITIVE=0

     ; Save the widget ids and other parameters to be accessed throughout 
     ; this widget application.  This state structure will be stored
     ; in the user value of the window and can be retreived through the 
     ; GET_UVALUE keyword of the IDL WIDGET_CONTROL procedure
   state = { $
             wEditorWindow : wEditorWindow, $
             wNewButton : wNewButton, $
             wOpenButton : wOpenButton, $
             wSaveButton : wSaveButton, $
             wSaveAsButton : wSaveAsButton, $
             wExitButton : wExitButton, $
             wTextEdit : wTextEdit, $
             wUppercase : wUppercase, $
             wLowercase : wLowercase, $
             wNextWord : wNextWord, $
             wCountWords : wCountWords, $
             wFindTextID: wFindTextID, $
             wFind : wFind, $
             wInfo : wInfo, $
             windowSize : [0,0], $
             textEditted : 0, $
             fileName : "Untitled" $
           }

     ; Make the window visible
   WIDGET_CONTROL, wEditorWindow, /REALIZE
   
     ; Display the wait cursor
   WIDGET_CONTROL, state.wEditorWindow, /HOURGLASS

     ; Get the text to be displayed in the text widget
   IF  (KEYWORD_SET(NEW)) THEN BEGIN 
       state.fileName = fileName 
; Set the title of the window to include the new file name
       WIDGET_CONTROL, state.wEditorWindow, $
         TLB_SET_TITLE="Editor - " + state.fileName
   ENDIF ELSE GetText, state, FILENAME=fileName 

     ; Get the current window size to be used when the user resizes the window
   WIDGET_CONTROL, wEditorWindow, TLB_GET_SIZE=windowSize

     ; Save it in the state structure
   state.windowSize = windowSize

     ; Save the state structure in the window's user value
   WIDGET_CONTROL, wEditorWindow, SET_UVALUE=state

     ; Register this widget application with the widget manager
   Xmanager, "Editor", wEditorWindow, GROUP_LEADER=group, $
     EVENT_HANDLER="EditorEventHdlr", CLEANUP="EditorKilled"

END  ;--------------------- procedure Editor ----------------------------
