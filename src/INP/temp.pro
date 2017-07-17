files = ['xw_bias_2006m0711_NMLS_o4ctpclr.dat', $
         'xw_bias_2006m0711_NMLS_o4ctpglbclr.dat', $
         'xw_bias_2006m0711_NMLS_o4ctpglbcld.dat']

files = ['xw_bias_2006m0711_NMLS_o4ctpglbclr.dat', $
         'xw_bias_2006m0711_NMLS_o4ctpglbcld.dat']


nfile = n_elements(files)
nw1 = 250  & nw2 = 536 & nw = nw1 + nw2
nx = 30
corr = dblarr(nfile, nx, nw)
wave = dblarr(nfile, nw)
line=''

for ifile = 0, nfile - 1 do begin
    close, 1 & openr, 1, files(ifile)
    readf, 1, line
    temp1 = dblarr(nx+1, 250)
    readf, 1, temp1
    wave(ifile, 0:249) = temp1(0, *)
    corr(ifile, *, 0:249) = temp1(1:nx, *)
   
    readf, 1, line
    temp2 = dblarr(nx*2+1, 536)
    readf, 1, temp2
    wave(ifile, 250:nw-1) = temp2(0, *)
 
    temp = dblarr(nx, 536)
    for ix = 0, nx - 1 do begin
        temp(ix, *) = (temp2(ix*2, *)+ temp2(ix*2+1, *)) / 2.0
    endfor
    corr(ifile, *, 250:nw-1) = temp(*, *)
endfor
wave = reform(wave(0, *))


stop

; write the difference into a separter file
close, 3 & openw, 3, 'xw_bias_2006m0711_NMLS_o4ctpglbdf.dat'
corrdf = reform(corr(1, *, *)-corr(0, *, *))

printf, 3, nx, nw1, ' ', ' UV1', format='(2I5, 2A5)'
for iw = 0, nw1 - 1 do begin
    printf, 3, wave(iw), corrdf(0:nx-1, iw), format='(F10.4, 100F10.6)'
endfor

printf, 3, nx, nw2, ' ', ' UV2', format='(2I5, 2A5)'
for jw = 0, nw2 - 1 do begin
    iw = nw1 + jw
    printf, 3, wave(iw), corrdf(0:nx-1, iw), format='(F10.4, 100F10.6)'
endfor
stop
end
