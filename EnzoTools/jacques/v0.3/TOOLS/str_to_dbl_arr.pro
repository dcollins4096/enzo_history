pro STR_TO_DBL_ARR, string, dbl_arr
; extract floating point array from a string that are separated by spaces
dbl_arr = DOUBLE(STR_SEP(strcompress(string),' '))
RETURN
END
