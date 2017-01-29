pro early, im, im2

; crop the array to match that of the experiment.

    ncxhr1 = 194
    ncxhr2 = 678
    ncxlr1 = 97
    ncxlr2 = 339
 
    ncx1 = ncxhr1
    ncx2 = ncxhr2

    im2 = im[ncx1:ncx2,*]

; end of routine
    return

end
