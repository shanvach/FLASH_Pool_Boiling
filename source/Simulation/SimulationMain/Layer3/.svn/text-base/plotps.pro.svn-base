pro plotps, base, im, XSIZE=xs, YSIZE=ys

outfile = base + '.ps'
current_device = !d.name
set_plot, 'PS'

xsize = 4.8
ysize = 4.8

imsize = size(im)
nx = imsize[1]
ny = imsize[2]


if (keyword_set(ys)) then begin
	ysize=ys
	if (keyword_set(xs)) then begin
		xsize=xs
	end else begin
		xsize=(1.*nx/(1.*ny))*ys
        end
end else begin
	if (keyword_set(xs)) then begin
		xsize=xs
                ysize=((1.*ny)/(1.*nx))*xs
        end
end


device, file = outfile, xsize=xsize, ysize=ysize, $
                xoff=0.5, yoff=1.75, /inch, /color, bits_per_pixel = 8
device,filename=outfile
tvscl,im
device,/close

set_plot, 'x'

return
end

