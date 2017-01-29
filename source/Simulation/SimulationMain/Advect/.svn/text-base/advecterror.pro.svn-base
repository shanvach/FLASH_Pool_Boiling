;
; advecterror
;  
;   pro advecterror, filename, RHOIN=rhoin, RHOOUT=rhoout, PRES=pres,
;                 VEL=vel, WIDTH=width, XI=xi, dens, compdens, drho, error,avgdx
;
;       Input
;         filename -- name of the FLASH file to read in
;         left  state (rhol, pl, ul)
;         right state (rhor, pr, ur)
;         gamma of fluid (gamma)
;         initial posn (xi)
;       Output
;         dens     -- density field  from analytic soln
;         compdens -- computed field from FLASH output
;         drho     -- difference of the two
;         error    -- integral over area of error
;         avgdx    -- average resolution in domain


pro advecterror, filename, RHOIN=rhoin, RHOOUT=rhoout, PRES=pres, VEL=vel, $
                 XI=xi, dens, compdens, drho, error, avgdx

    ;
    ; default data
    ;
    if (not(keyword_set(rhoin)))  then rhoin = 1.
    if (not(keyword_set(rhoout))) then rhoout = 1.e-5
    if (not(keyword_set(pres)))   then pres = 1.
    if (not(keyword_set(vel)))    then vel =  10.
    if (not(keyword_set(width)))  then width = .1
    if (not(keyword_set(xi)))     then xi = .25


    ; 
    ; open the file, and get bounding information
    ;
    openflashfile, filename, state
    flashgetbb, state, bb, mindx, maxdx, maxlref, minlref
    idens = where(state.unknames eq "dens")

    ;
    ; set up the 1d analytic soln
    ;

    nxpts = (bb[1,0]-bb[0,0])/mindx[0]
    x = (findgen(nxpts)+.5)*mindx[0] + bb[0,0]
    t = state.params.time
    advect, rhoin, rhoout, pres, width, vel, xi, t, x, psoln, rsoln, vsoln

    avgdx = 0.
    nblocks = 0L

	;
    ; copy the whole file info into dens and compdens

    dens = state
	compdens = state

    ;
    ; loop over blocks
    ;
    error = 0.
    flashnextblock, state, block
    while (n_tags(block) gt 0) do begin
        if (block.nodetype eq 1) then begin
        nxb = block.nxb

        dx = block.xznr[0]-block.xznl[0]

        avgdx = avgdx + dx
        nblocks = nblocks + 1L

        coarsendata, rsoln, x, block.xznl, block.xznr, analytic1d

        if (state.params.ndim eq 1) then begin
	        analyticsolndata = dblarr(nxb)
 
           	analyticsolndata[*] = analytic1d[*]

;         	error = error + total(abs(block.solndata[idens,*]-analyticsolndata))*dx
         	error = error + sqrt(total((block.solndata[idens,*]-analyticsolndata)^2))*dx
         	compdens.unk[*,    block.blocknum,*] = 0.
         	compdens.unk[idens,block.blocknum,*] = analyticsolndata[*]
        endif

        if (state.params.ndim eq 2) then begin
            nyb = block.nxb
            dy = block.yznr[0]-block.yznl[0]
	        analyticsolndata = dblarr(nxb,nyb)
 
         	for j=0, nyb-1 do begin 
           	 analyticsolndata[*,j] = analytic1d[*]
         	endfor

;         	error = error + total(abs(block.solndata[idens,*,*]-analyticsolndata))*dx*dy
;         	error = error + sqrt(total((block.solndata[idens,*,*]-analyticsolndata)^2))*dx*dy
         	compdens.unk[*,    block.blocknum,*,*] = 0.
         	compdens.unk[idens,block.blocknum,*,*] = analyticsolndata[*,*]
         endif

         endif
         flashnextblock, state, block
    endwhile
 
    avgdx = avgdx/(nblocks*1.)
	return

end


pro coarsendata,pres, xi, xznl, xznr, coarsepres
	
	coarsepres = xznl * 0.

    for i=0, n_elements(xznl)-1 do begin
		pstart = min(where(xi ge xznl[i]))
		pend   = max(where(xi le xznr[i]))
        coarsepres[i] = mean(pres[pstart:pend])
    endfor


    return
end

