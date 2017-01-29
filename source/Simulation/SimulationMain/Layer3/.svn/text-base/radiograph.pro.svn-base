pro radiograph, filename, im

    cuopacity = 1.0
    chopacity = 0.6
    cfopacity = 0.3

    xcenter = .10
    ycenter = .045

    radius  = .043
    imagedx = .0005

    noiseamp = .2

    cu = loaddata(filename, "cu  ", XCOORDS=xvals, YCOORD=yvals)
    ch = loaddata(filename, "ch  ")
    cf = loaddata(filename, "cf  ")

    tmp = size(cu)
    nx = tmp[1]
    ny = tmp[2]

    xvals = xvals - xcenter
    yvals = yvals - ycenter
 
    dx = xvals[1]-xvals[0]
    dy = yvals[1]-yvals[0]

    dist = cu
    intensity = cu

    opacity   = cuopacity*cu + chopacity*ch + cfopacity*cf

    for i=0,nx-1 do begin 
            for j=0,ny-1 do begin 
                dist[i,j] = xvals[i]*xvals[i] + yvals[j]*yvals[j]
            end
    end

    dist = sqrt(dist)/radius

    intensity = intensity * 0.

    intensity[where(dist le 1)] = 1
    intensity[where(dist gt 1 and dist le 1.5)] = (1.5 - dist[where(dist gt 1 and dist le 1.5)])/.5
;    intensity[where(dist gt 1 and dist le 1.3)] = (1.3 - dist[where(dist gt 1 and dist le 1.3)])/.3
    
;    intensity =  exp(-dist*dist/2.)^2

    im = intensity * (1.-opacity )

    noise = randomu(1, nx*dx/imagedx, ny*dy/imagedx)*2.*noiseamp 
    noise = noise - mean(noise)
   
    im = (congrid(im, nx*dx/imagedx, ny*dy/imagedx, /interp)) +noise
    im = congrid(im, nx*dx/dx, ny*dy/dx) 

    im[where(im lt 0)] = 0 

    return

end
