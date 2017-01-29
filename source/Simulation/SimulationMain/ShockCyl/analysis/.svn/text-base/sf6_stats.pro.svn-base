; assuming cartesian coordinates, and 2d data, read in a checkpoint 
; file from the shockcylinder run and calculate and print 
; some diagnostic information:
; 
;  checkpoint simulation time,
;  total sf6 mass, 
;  highest zone density of sf6 with coordinates, 
;  position of center of mass of sf6, 

pro sf6_stats, filename

if n_elements(var_name) EQ 0 then var_name = 'SF6'

print, filename, var_name
 
varnames = get_var_list(filename)
index = (where(varnames EQ 'sf6 '))[0]
if index GT 0 then begin
    sf6name = 'sf6 ' ; flash2 run
end else begin
    sf6name = 'SF6 ' ; otherwise, flash3 run
end 


read_amr, filename, VAR_NAME=sf6name, $
          TREE=tree, DATA=sf6, PARAMETERS=params

read_amr, filename, VAR_NAME='dens', DATA=dens


x = dblarr(params.nxb)
y = dblarr(params.nyb)

xcm = 0.
ycm = 0.
mtot = 0.
sf6densmax = 0.
xsf6max = 0.
ysf6max = 0.
for block = 0, params.totBlocks-1 do begin
    if (tree[block].nodeType EQ 1) then begin
        xmin = tree[block].bndBox[0,0]
        ymin = tree[block].bndBox[0,1]
        dx = tree[block].size[0]/params.nxb
        dy = tree[block].size[1]/params.nyb
        x = dx * (dindgen(params.nxb) + 0.5) + xmin
        y = dy * (dindgen(params.nyb) + 0.5) + ymin

        for j = 0, params.nyb - 1 do begin
            for i = 0, params.nxb - 1 do begin
                
                sf6dens = dens[0,block,i,j] * sf6[0,block,i,j]
                dm = dx*dy * sf6dens
             
                if sf6dens GT sf6densmax then begin
                    sf6densmax = sf6dens
                    xsf6max = x[i]
                    ysf6max = y[j]
                end
                mtot = mtot + dm
                xcm = xcm + (dm * x[i]) 
                ycm = ycm + (dm * y[j])
            end 
        end
    end
end

xcm = xcm / mtot
ycm = ycm / mtot
print, 'time', params.time
print, 'total mass', mtot
print, 'max sf6', sf6densmax, xsf6max, ysf6max
print, 'center of mass: ', xcm, ycm
end  
