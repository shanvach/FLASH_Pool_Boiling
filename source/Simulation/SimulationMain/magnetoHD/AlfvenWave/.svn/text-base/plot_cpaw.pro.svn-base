a=atan(2.)
cos_a= cos(a)
sin_a= sin(a)

; Traveling waves

; initial condition
read_hdf5_procedure,'128x64_hdf5_chk_0000','xnor',xn128,params128
read_hdf5_procedure,'128x64_hdf5_chk_0000','btra',bt128_t0,params128

; 128x64 solution at time t = 5
read_hdf5_procedure,'128x64_hdf5_chk_0005','magx',bx128_t5,params128
read_hdf5_procedure,'128x64_hdf5_chk_0005','magy',by128_t5,params128
bt128_t5 = -bx128_t5*sin_a + by128_t5*cos_a

; 64x32 solution at time t = 5
read_hdf5_procedure,'64x32_hdf5_chk_0005','xnor',xn64,params64
read_hdf5_procedure,'64x32_hdf5_chk_0005','magx',bx64_t5,params64
read_hdf5_procedure,'64x32_hdf5_chk_0005','magy',by64_t5,params64
bt64_t5 = -bx64_t5*sin_a + by64_t5*cos_a

; 32x16 solution at time t = 5
read_hdf5_procedure,'32x16_hdf5_chk_0005','xnor',xn32,params32
read_hdf5_procedure,'32x16_hdf5_chk_0005','magx',bx32_t5,params32
read_hdf5_procedure,'32x16_hdf5_chk_0005','magy',by32_t5,params32
bt32_t5 = -bx32_t5*sin_a + by32_t5*cos_a

; 16x8 solution at time t = 5
read_hdf5_procedure,'16x8_hdf5_chk_0005','xnor',xn16,params16
read_hdf5_procedure,'16x8_hdf5_chk_0005','magx',bx16_t5,params16
read_hdf5_procedure,'16x8_hdf5_chk_0005','magy',by16_t5,params16
bt16_t5 = -bx16_t5*sin_a + by16_t5*cos_a

; 8x4 solution at time t = 5
read_hdf5_procedure,'8x4_hdf5_chk_0005','xnor',xn8,params8
read_hdf5_procedure,'8x4_hdf5_chk_0005','magx',bx8_t5,params8
read_hdf5_procedure,'8x4_hdf5_chk_0005','magy',by8_t5,params8
bt8_t5 = -bx8_t5*sin_a + by8_t5*cos_a


loadct,1
winb
;plot,xn128,bt128_t0,psym=1
;oplot,xn128,bt128_t5,psym=2
;oplot,xn64,bt64_t5,psym=3
;oplot,xn32,bt32_t5,psym=4
;oplot,xn16,bt16_t5,psym=5
;oplot,xn8, bt8_t5, psym=6

n=64
plot,xn128(*,0),bt128_t0(*,0),xr=[0,2],lines=1
for i=1,n-1 do begin
	oplot,xn128(*,i),bt128_t0(*,i),lines=1
	oplot,xn128(*,i),bt128_t5(*,i),lines=2
endfor

n=32
for i=0,n-1 do begin
	oplot,xn64(*,i),bt64_t5(*,i),lines=3
endfor

n=16
for i=0,n-1 do begin
	oplot,xn32(*,i),bt32_t5(*,i),lines=4
endfor

n=8
for i=0,n-1 do begin
	oplot,xn16(*,i),bt16_t5(*,i),lines=5
endfor

n=4
for i=0,n-1 do begin
	oplot,xn8(*,i),bt8_t5(*,i),lines=6
endfor



; Standing waves

read_hdf5_procedure,'steady_128x64_hdf5_chk_0000','xnor',xn128s,params128s
read_hdf5_procedure,'steady_128x64_hdf5_chk_0000','btra',bt128_t0s,params128s

; 128x64 solution at time t = 5
read_hdf5_procedure,'steady_128x64_hdf5_chk_0005','magx',bx128_t5s,params128s
read_hdf5_procedure,'steady_128x64_hdf5_chk_0005','magy',by128_t5s,params128s
bt128_t5s = -bx128_t5s*sin_a + by128_t5s*cos_a

; 64x32 solution at time t = 5
read_hdf5_procedure,'steady_64x32_hdf5_chk_0005','xnor',xn64s,params64s
read_hdf5_procedure,'steady_64x32_hdf5_chk_0005','magx',bx64_t5s,params64s
read_hdf5_procedure,'steady_64x32_hdf5_chk_0005','magy',by64_t5s,params64s
bt64_t5s = -bx64_t5s*sin_a + by64_t5s*cos_a

; 32x16 solution at time t = 5
read_hdf5_procedure,'steady_32x16_hdf5_chk_0005','xnor',xn32s,params32s
read_hdf5_procedure,'steady_32x16_hdf5_chk_0005','magx',bx32_t5s,params32s
read_hdf5_procedure,'steady_32x16_hdf5_chk_0005','magy',by32_t5s,params32s
bt32_t5s = -bx32_t5s*sin_a + by32_t5s*cos_a

; 16x8 solution at time t = 5
read_hdf5_procedure,'steady_16x8_hdf5_chk_0005','xnor',xn16s,params16s
read_hdf5_procedure,'steady_16x8_hdf5_chk_0005','magx',bx16_t5s,params16s
read_hdf5_procedure,'steady_16x8_hdf5_chk_0005','magy',by16_t5s,params16s
bt16_t5s = -bx16_t5s*sin_a + by16_t5s*cos_a

; 8x4 solution at time t = 5
read_hdf5_procedure,'steady_8x4_hdf5_chk_0005','xnor',xn8s,params8s
read_hdf5_procedure,'steady_8x4_hdf5_chk_0005','magx',bx8_t5s,params8s
read_hdf5_procedure,'steady_8x4_hdf5_chk_0005','magy',by8_t5s,params8s
bt8_t5s = -bx8_t5s*sin_a + by8_t5s*cos_a


loadct,1
winb,1
;plot,xn128s,bt128_t0s,psym=1
;oplot,xn128s,bt128_t5s,psym=2
;oplot,xn64s,bt64_t5s,psym=3
;oplot,xn32s,bt32_t5s,psym=4
;oplot,xn16s,bt16_t5s,psym=5
;oplot,xn8s, bt8_t5s, psym=6

n=64
plot,xn128s(*,0),bt128_t0s(*,0),xr=[0,2],lines=1
for i=1,n-1 do begin
	oplot,xn128s(*,i),bt128_t0s(*,i),lines=1
	oplot,xn128s(*,i),bt128_t5s(*,i),lines=2
endfor

n=32
for i=0,n-1 do begin
	oplot,xn64s(*,i),bt64_t5s(*,i),lines=3
endfor

n=16
for i=0,n-1 do begin
	oplot,xn32s(*,i),bt32_t5s(*,i),lines=4
endfor

n=8
for i=0,n-1 do begin
	oplot,xn16s(*,i),bt16_t5s(*,i),lines=5
endfor

n=4
for i=0,n-1 do begin
	oplot,xn8s(*,i),bt8_t5s(*,i),lines=6
endfor

end
