pro measure, filename, t, spikelen, gaplen, bubblelen, IMAGE=im

    t = 1.
    im = loaddata(filename, "ch  ", XCOORDS=xval, YCOORD=yval, TIME=t)

    dims = size(im)
    nx = dims[1]
    ny = dims[2]

    avg = findgen(nx)
 
    for i=0,nx-1 do  avg[i] = mean(im[i,*])
    navg = nx/40
    avg = smooth(avg,navg)
;    avg = smooth(avg,40)

    ispikestart = min(where(avg gt .05))
    spikestart  = xval[ispikestart]

    ispikeend   = min(where(avg gt .9))
    print, 'ispikeend' 
    print, ispikeend
    spikeend    = xval[ispikeend]
    spikelen    = spikeend - spikestart

    ibubstart   = max(where(avg gt .9))
    bubstart    = xval[ibubstart]

    ibubend   = max(where(avg gt .25))
    bubend    = xval[ibubend]
    bubblelen = bubend-bubstart

    gaplen = bubstart - spikeend

    im[ispikestart,*] = 2 
    im[ispikeend  ,*] = 2 
    im[ibubstart  ,*] = 2 
    im[ibubend    ,*] = 2 

    return
end
