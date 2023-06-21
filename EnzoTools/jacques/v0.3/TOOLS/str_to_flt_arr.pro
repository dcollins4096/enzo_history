pro STR_TO_FLT_ARR, string, flt_arr
; extract floating point array from a string that are separated by spaces
flt_arr = FLOAT(STR_SEP(strcompress(string),' '))
RETURN
END
