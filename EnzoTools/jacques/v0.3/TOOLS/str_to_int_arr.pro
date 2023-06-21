pro STR_TO_INT_ARR, string, int_arr
; extract integers from a string that are separated by spaces

;remove leading and trailing spaces
int_arr = FIX(STR_SEP(strcompress(string), ' '))
END
