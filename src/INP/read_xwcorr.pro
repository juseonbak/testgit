; pro read_xwcorr

nx1=0 & nx2 = 0  & nw1 = 0 & nw2 = 0
datadir = '/home/xliu/OzoneFit/OZBOREAS-OMI/src/INP/'
;datadir = '/home/xliu/OzoneFit/GLBO3/OMIPROGS/
close, 1 & openr, 1, datadir + 'xw_bias_2006m07_hres_8s48lnr.dat'

readf, 1, nx1, nw1
wave1=fltarr(nw1)         ; wavelength in UV 1
corr1=fltarr(nx1, nw1)    ; ratio of observed radiance to simulated radiance (UV1)
temp=fltarr(nx1+1, nw1)
readf, 1, temp
wave1(*) = temp(0, *)
corr1(*) = temp(1:nx1, *)

; read data in UV 2
readf, 1, nx2, nw2
wave2=fltarr(nw2)         ; wavelength in UV 2
corr2=fltarr(nx2, nw2)    ; ratio of observed radiance to simulated radiance (UV2)
temp=fltarr(nx2+1, nw2)
readf, 1, temp
wave2(*) = temp(0, *)
corr2(*) = temp(1:nx2, *)

;nx1=0 & nx2 = 0  & nw1 = 0 & nw2 = 0
datadir = '/home/xliu/OzoneFit/OZBOREAS-OMI/src/INP/'
datadir = '/home/xliu/OzoneFit/GLBO3/OMIPROGS/
close, 1 & openr, 1, datadir + 'xw_bias_2006m07_hres_8s48lnr_glb.dat'

; read data in UV 1
readf, 1, nhx1, nhw1
hwave1=fltarr(nhw1)         ; wavelength in UV 1
hcorr1=fltarr(nhx1, nhw1)    ; ratio of observed radiance to simulated radiance (UV1)
temp=fltarr(nhx1+1, nhw1)
readf, 1, temp
hwave1(*) = temp(0, *)
hcorr1(*) = temp(1:nhx1, *)

; read data in UV 2
readf, 1, nhx2, nhw2
hwave2=fltarr(nhw2)         ; wavelength in UV 2
hcorr2=fltarr(nhx2, nhw2)    ; ratio of observed radiance to simulated radiance (UV2)
temp=fltarr(nhx2+1, nhw2)
readf, 1, temp
hwave2(*) = temp(0, *)
hcorr2(*) = temp(1:nhx2, *)

; plot the differences between these two correction
diffcorr1 = fltarr(nx1, nw1)  & diffcorr2 = fltarr(nx1, 536)
diffcorr1 = hcorr1 - corr1
tempcorr2 = hcorr2(*, 0:535) - corr2(*, 0:535)
acorr2  = fltarr(nx1, 536)
ahcorr2=acorr2

for i = 0, nx1-1 do begin
    ;temp1 = interpol(hcorr2(i*2, *), hwave2, wave2)-corr2(i*2, *)
    ;temp2 = interpol(hcorr2(i*2+1, *), hwave2, wave2) - corr2(i*2+1, *)
    ;diffcorr2(i, *) = ( temp1 + temp2 ) / 2.0

   diffcorr2(i, *) = (tempcorr2(i*2, *) + tempcorr2(i*2+1, *)) / 2.
   acorr2(i, 0:535) = (corr2(i*2, 0:535) + corr2(i*2+1, 0:535)) / 2.
   ahcorr2(i, 0:535) = (hcorr2(i*2, 0:535) + hcorr2(i*2+1, 0:535)) / 2.
endfor

; plot the difference in correction: figure out what causes it
; Could it make a large difference?
; Where do the structures come from? "fitted parameters"
figdir = '/data/dumbo/xliu/OMIFIGS/PAPERS/'
figname = figdir + 'omicorr_hres-effcrs_corrdf.ps'

set_plot, 'ps'
device, file=figname, /portrait, /inches, /color, xsize=7.5, ysize=9.5, xoffset=0.5, yoffset=0.5
!p.thick=3 & !x.thick=3 & !y.thick=3 & !p.charsize=1.25 & !p.multi = [0, 1, 1]
restore,'/home/xliu/TOOLS/colors/idl_coltab25.dat'  & tvlct, r, g, b

plot, wave1, hcorr1(0, *)*100-100., xrange = [270, 350], xstyle=1, yrange=[-10, 15], ystyle=1, $
      title = '!6', color=1, thick=3, charthick=3.0, /nodata, position=[0.1, 0.3, 0.9, 0.7], $
      charsize=1.25, xtitle='!6Wavelength (nm)', ytitle = '!6HRES'
oplot, [270, 350], [0, 0], linestyle=1, color=1, thick=3
oplot, [310., 310.], [-10.0, 15.0], linestyle=1, color=1, thick=3
xyouts, 0.33, 0.67, '!6UV-1', color=1, charsize=1.25, charthick=3.0, /normal
xyouts, 0.60, 0.67, '!6UV-2', color=1, charsize=1.25, charthick=3.0, /normal

dcol = 1  & dcol = 252. / (nx1 - 1)
cols = findgen(nx1) * dcol + 2.
for ix = 0, nx1-1 do begin
    oplot, wave1, hcorr1(ix, *)*100.-100., color=cols(ix), thick=3
    oplot, wave2, ahcorr2(ix, *)*100.-100., color=cols(ix), thick=3
endfor

barlvl = indgen(nx1) + 1
barval = string(indgen(nx1)+1, format='(I2.2)')
da = where(barlvl mod 5 ne 0)
barval(da) = ' '
vert_color_bar, 0.91, 0.3, 0.93, 0.7, barval, cols, FORM='(A3)', thick_own=3.0, size_own=1.25, color_own=1


plot, wave1, corr1(0, *)*100-100., xrange = [270, 350], xstyle=1, yrange=[-10, 15], ystyle=1, $
      title = '!6', color=1, thick=3, charthick=3.0, /nodata, position=[0.1, 0.3, 0.9, 0.7], $
      charsize=1.25, xtitle='!6Wavelength (nm)', ytitle = '!6EFFCRS'
oplot, [270, 350], [0, 0], linestyle=1, color=1, thick=3
oplot, [310., 310.], [-10.0, 15.0], linestyle=1, color=1, thick=3
xyouts, 0.33, 0.67, '!6UV-1', color=1, charsize=1.25, charthick=3.0, /normal
xyouts, 0.60, 0.67, '!6UV-2', color=1, charsize=1.25, charthick=3.0, /normal

dcol = 1  & dcol = 252. / (nx1 - 1)
cols = findgen(nx1) * dcol + 2.
for ix = 0, nx1-1 do begin
    oplot, wave1, corr1(ix, *)*100.-100., color=cols(ix), thick=3
    oplot, wave2, acorr2(ix, *)*100-100., color=cols(ix), thick=3
endfor

barlvl = indgen(nx1) + 1
barval = string(indgen(nx1)+1, format='(I2.2)')
da = where(barlvl mod 5 ne 0)
barval(da) = ' '
vert_color_bar, 0.91, 0.3, 0.93, 0.7, barval, cols, FORM='(A3)', thick_own=3.0, size_own=1.25, color_own=1


plot, wave1, diffcorr1(0, *)*100, xrange = [270, 350], xstyle=1, yrange=[-3, 3], ystyle=1, $
      title = '!6', color=1, thick=3, charthick=3.0, /nodata, position=[0.1, 0.3, 0.9, 0.7], $
      charsize=1.25, xtitle='!6Wavelength (nm)', ytitle = '!6HRES - EFFCRS Corr. Diff.'
oplot, [270, 350], [0, 0], linestyle=1, color=1, thick=3
oplot, [310., 310.], [-3.0, 3.0], linestyle=1, color=1, thick=3
xyouts, 0.33, 0.67, '!6UV-1', color=1, charsize=1.25, charthick=3.0, /normal
xyouts, 0.60, 0.67, '!6UV-2', color=1, charsize=1.25, charthick=3.0, /normal

dcol = 1  & dcol = 252. / (nx1 - 1)
cols = findgen(nx1) * dcol + 2.
for ix = 0, nx1-1 do begin
    oplot, wave1, diffcorr1(ix, *)*100., color=cols(ix), thick=3
    oplot, wave2, diffcorr2(ix, *)*100, color=cols(ix), thick=3
endfor

barlvl = indgen(nx1) + 1
barval = string(indgen(nx1)+1, format='(I2.2)')
da = where(barlvl mod 5 ne 0)
barval(da) = ' '
vert_color_bar, 0.91, 0.3, 0.93, 0.7, barval, cols, FORM='(A3)', thick_own=3.0, size_own=1.25, color_own=1

device, /close
stop

end

