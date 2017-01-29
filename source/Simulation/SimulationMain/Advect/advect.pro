; 
; advect -- evolve a gaussian pulse advection
;
;  inputs:  rhoin, rhoout
;           pressure
;           width of pulse
;           velocity of pulse
;           initial position of pulse center
;                       
;  outputs:  thermodynamic quantities on the mesh at time t:
;              psoln, rsoln, vsoln

;
pro advect, rhoin, rhoout, pres, width, vel, xi, t, xpts, psoln, rsoln, vsoln

      subsample = 50                ; # pts/grid point to average over
                                    ; this should be even!
      npts      = n_elements(xpts)
      nsspts    = npts*subsample+1L

      psoln = xpts*0.
      rsoln = xpts*0.
      vsoln = xpts*0.
  
      rho = findgen(nsspts)*0.
      u   = findgen(nsspts)*0.
      p   = findgen(nsspts)*0.

      dx = xpts[2] - xpts[1]       ; assuming a uniform mesh!!!
      xl = xpts[0]-.5*dx
      xr = xpts[npts-1]+.5*dx
 
      xznl = xpts - xpts[0] + xl
      xznr = xpts - xpts[0] + xl + dx

      xsspts = findgen(nsspts)*(dx/(1.*subsample)) + xl
      
;
;--------interval over which to compute solution
;
     if (xr lt xl) then begin
         print, 'xr must be greater than xl'
      end
;
;
;-----begin solution; impose periodic BC
;
;
      xscale = (xsspts - vel*t)
      toolow =where(xscale lt xl,lowcount)  
      toohigh=where(xscale gt xr,hicount )  

      if (lowcount gt 0) then begin
          xscale[toolow] = xscale[toolow] + (xr-xl)
      endif
      if (hicount gt 0) then begin
          xscale[toohigh] = xscale[toohigh] + (-xr+xl)
       endif

      xscale=(xscale-xi)/width
      wfac = exp(-(xscale^2.))

      rho = rhoin*wfac + rhoout*(1.-wfac)
      p   = p + pres
      u   = u + vel
;

      for i=ulong(0),ulong(npts-1) do begin
          pstart = min(where(xsspts ge xznl[i]))
          pend   = max(where(xsspts le xznr[i]))
          rsoln(i) = mean(rho[pstart:pend])
          psoln(i) = mean(p[pstart:pend])
          vsoln(i) = mean(u[pstart:pend])
      end

end

