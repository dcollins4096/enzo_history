pro READ_STRING_FILE, file_name, file_in_string
; read a text file into a string
  file_in_string = ''
  help_string    = ''
  i = 0
  openr, unit, file_name, /get_lun
  WHILE (NOT eof(unit)) DO begin
      i = i+ 1
      print, i
      readf, unit, help_string
      file_in_string = file_in_string + help_string
  ENDWHILE

  print,file_in_string
END

