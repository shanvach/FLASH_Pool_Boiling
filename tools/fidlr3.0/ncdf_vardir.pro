FUNCTION NCDF_VARDIR, CDFID
; from Gumley book, p.177
; returns a string array containing the names of all variables in an open netCDF file.
; the null string '' is returned if no variables are found.


;- Check arguments
if (n_params() ne 1) then $
  message, 'Usage: RESULT = NCDF_VARDIR(CDFID)'
if (n_elements(cdfid) eq 0) then $
  message, 'Argument CDFID is undefined'

;- Set default return value
varnames = ''

;- Get file information
fileinfo = ncdf_inquire(cdfid)
nvars = fileinfo.nvars

;- If variables were found, get variable names
if (nvars gt 0) then begin
  varnames = strarr(nvars)
  for index = 0L, nvars - 1L do begin
    varinfo = ncdf_varinq(cdfid, index)
    varnames[index] = varinfo.name
  endfor
endif

;- Return the result
return, varnames

END
