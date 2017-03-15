pro omi_hcho

path = "D:\OMI\in\"
rad_name = "OMI-Aura_L1-OML1BRUG_2005m0921t0529-o06305_v003-2011m0120t195711-p1.he4"
irr_name = "OMI-Aura_L1-OML1BIRR_2005m0921t2337-o06316_v003-2007m0417t023751.he4"

;---------------------------------------------------------------------------------------------------
fid = hdf_open(path + rad_name)
sd_id = hdf_sd_start(path + rad_name, /read)

vg_id = hdf_vg_lone(fid)
vds = n_elements(vg_id)
for i = 0, 1 do begin
  vg_h = hdf_vg_attach(fid, vg_id[i])
  hdf_vg_getinfo, vg_h, name = name_vg
  if name_vg eq "Earth UV-2 Swath" then vg_h = vg_h
endfor

hdf_vg_gettrs, vg_h, Tags, Refs
for i = 0, n_elements(Refs) - 1 do begin
  vg_h = hdf_vg_attach(fid, Refs[i])
  hdf_vg_getinfo, vg_h, name = name_vg
  if name_vg eq 'Geolocation Fields' then geo_h = vg_h
  if name_vg eq 'Data Fields' then dat_h = vg_h
endfor

hdf_vg_gettrs, dat_h, Tags_data, Refs_data
for i = 0, n_elements(Refs_data) - 1 do begin
  sd = hdf_sd_reftoindex(sd_id, Refs_data[i])
  if sd ne -1 then begin
    id = hdf_sd_select(sd_id, sd)
    hdf_sd_getinfo, id, NAME = name_sd
    if name_sd eq 'RadianceMantissa' then radman_id = id
    if name_sd eq 'RadianceExponent' then radexp_id = id
    if name_sd eq 'WavelengthCoefficient' then radwc_id = id
  endif
endfor

;get radiance data
hdf_sd_getdata, radman_id, radman
hdf_sd_endaccess, radman_id

hdf_sd_getdata, radexp_id, radexp
hdf_sd_endaccess, radexp_id

hdf_sd_getdata, radwc_id, radwc
hdf_sd_endaccess, radwc_id

;read the data of Latitude and Longitude
hdf_vg_gettrs, geo_h, Tags_data, Refs_data
for i = 0, n_elements(Refs_data) - 1 do begin
  sd = hdf_sd_reftoindex(sd_id, Refs_data[i])
  if sd ne -1 then begin
    id = hdf_sd_select(sd_id, sd)
    hdf_sd_getinfo, id, NAME = name_sd
    if name_sd eq 'Latitude' then lat_id = id
    if name_sd eq 'Longitude' then lon_id = id
    if name_sd eq 'SolarZenithAngle' then sza_id = id
    if name_sd eq 'ViewingZenithAngle' then vza_id = id
  endif
endfor

;get the data of Latitude and Longitude
hdf_sd_getdata, lat_id, lat
hdf_sd_endaccess, lat_id

hdf_sd_getdata, lon_id, lon
hdf_sd_endaccess, lon_id

hdf_sd_getdata, sza_id, sza
hdf_sd_endaccess, sza_id

hdf_sd_getdata, vza_id, vza
hdf_sd_endaccess, vza_id

hdf_sd_end, sd_id
hdf_close, fid
;--------------------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------------------
fid = hdf_open(path + irr_name)
sd_id = hdf_sd_start(path + irr_name, /read)

vg_id = hdf_vg_lone(fid)
vds = n_elements(vg_id)
for i = 0, 1 do begin
  vg_h = hdf_vg_attach(fid, vg_id[i])
  hdf_vg_getinfo, vg_h, name = name_vg
  if name_vg eq "Sun Volume UV-2 Swath" then vg_h = vg_h
endfor

hdf_vg_gettrs, vg_h, Tags, Refs
for i = 0, n_elements(Refs) - 1 do begin
  vg_h = hdf_vg_attach(fid, Refs[i])
  hdf_vg_getinfo, vg_h, name = name_vg
  if name_vg eq 'Geolocation Fields' then geo_h = vg_h
  if name_vg eq 'Data Fields' then dat_h = vg_h
endfor

hdf_vg_gettrs, dat_h, Tags_data, Refs_data
for i = 0, n_elements(Refs_data) - 1 do begin
  sd = hdf_sd_reftoindex(sd_id, Refs_data[i])
  if sd ne -1 then begin
    id = hdf_sd_select(sd_id, sd)
    hdf_sd_getinfo, id, NAME = name_sd
    if name_sd eq 'IrradianceMantissa' then irrman_id = id
    if name_sd eq 'IrradianceExponent' then irrexp_id = id
    if name_sd eq 'WavelengthCoefficient' then irrwc_id = id
  endif
endfor

;get irradiance data
hdf_sd_getdata, irrman_id, irrman
hdf_sd_endaccess, irrman_id

hdf_sd_getdata, irrexp_id, irrexp
hdf_sd_endaccess, irrexp_id

hdf_sd_getdata, irrwc_id, irrwc
hdf_sd_endaccess, irrwc_id

hdf_sd_end, sd_id
hdf_close, fid
;--------------------------------------------------------------------------------------------------------
;hcho = make_array(2, file_lines(path + "hcho_clip.txt")) 
xsc = read_ascii(path + "hcho_clip.txt")
hcho = xsc.(0)
xsc = read_ascii(path + "bro_clip.txt")
bro = xsc.(0)
xsc = read_ascii(path + "oclo_clip.txt")
oclo = xsc.(0)
xsc = read_ascii(path + "o3_clip.txt")
o3 = xsc.(0)
xsc = read_ascii(path + "o4_clip.txt")
o4 = double(xsc.(0))
xsc = read_ascii(path + "no2_clip.txt")
no2 = double(xsc.(0))
temp = read_ascii(path + "amf_reshape.txt")
amf = temp.(0)

sza_model = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, $
            21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, $
            41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, $
            61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, $
            81, 82, 83, 84, 85, 86, 87, $
            87.1, 87.2, 87.3, 87.4, 87.5, 87.6, 87.7, 87.8, 87.9, $
            88, 88.1, 88.2, 88.3, 88.4, 88.5, 88.6, 88.7, 88.8, 88.9, $
            89, 89.1, 89.2, 89.3, 89.4, 89.5, 89.6, 89.7, 89.75, 89.8, 89.85, 89.9]
vza_model = [0, 5, 10, 15, 20, 25, 35, 45, 55, 60, 65]
;--------------------------------------------------------------------------------------------------------
dim = size(radman, /dimension)
n_ban = dim[0]
n_col = dim[1]
n_row = dim[2]
vcd = double(make_array(n_col, n_row))

radman = double(radman) & radexp = double(radexp)
irrman = double(irrman) & irrexp = double(irrexp)
rad = double(radman * (10.0 ^ radexp))
irr = double(irrman * (10.0 ^ irrexp))

for i = 0, n_col -1 do begin 
  for j = 0, n_row - 1 do begin 
;    i = 59
;    j = 1639
    rad_ = rad[*, i, j]
    irr_ = irr[*, i]
    rad_wl = cal_wavlen(radwc[*, i, j])
    if (rad_wl[0] lt 300 or rad_wl[-1] gt 390) then begin
      vcd[i, j] = 0
      continue
    endif
    ind = where(rad_wl gt 328.5 and rad_wl lt 356.5)
    rad_wl = rad_wl[ind]
    rad_ = rad_[ind]
    irr_wl = cal_wavlen(irrwc[*, i])
    irr_ = interpol(irr_, irr_wl, rad_wl)
        
    cosza = abs(cos(sza[i, j] /180 * !pi))
    rad_ = smooth(rad_, 9)
    irr_ = smooth(irr_, 9)
    ref = (!pi * rad_) / (cosza * irr_)
    ref[where(ref lt 0.0 or ref gt 1.0)] = -!values.f_nan
    ref_ = -alog(ref)
    pos = where(finite(ref_))
    s = size(pos)
    if (pos[0] eq -1 or s[1] lt 90) then begin
      vcd[i, j] = 0
      continue
    endif
    ref_ = ref_[pos]
    rad_wl = rad_wl[pos]
    
    y = diff_poly(rad_wl, ref_)
        
    x1 = double(interpol(hcho[1, *], hcho[0, *], rad_wl))
    x1 = smooth(x1, 7)
    x1 = diff_poly(rad_wl, x1)
    x2 = double(interpol(bro[1, *], bro[0, *], rad_wl))
    x2 = diff_poly(rad_wl, x2)
    x3 = double(interpol(oclo[1, *], oclo[0, *], rad_wl))
    x3 = diff_poly(rad_wl, x3)
    x4 = double(interpol(o3[1, *], o3[0, *], rad_wl))
    x4 = diff_poly(rad_wl, x4)
    x5 = double(interpol(o4[1, *], o4[0, *], rad_wl))
    x5 = diff_poly(rad_wl, x5)
;    x5 = x5 / 10.0^40
    x6 = double(interpol(no2[1, *], no2[0, *], rad_wl))
    x6 = diff_poly(rad_wl, x6)
    
    x = [transpose(x1), transpose(x2), transpose(x3), transpose(x4)]
    
    res = regress(x, y, CHISQ = goodness, SIGMA=sigma, CONST=const)
    ;111111111111111111111111111111111
;    x = transpose(x)
;    c = double([[10.0^15], [10.0^13], [10.0^14], [10.0^18]])
;    bc = [1.0]
;    contype = [1]
;    xub_ = double([10.0^17, 10.0^15, 10.0^15, 10.0^19])
;    xlb_ = double([-10.0^17, -10.0^15, -10.0^15, -10.0^19])
;    res = IMSL_LINLSQ(y, x, c, bc, bc, contype, Xlb = xlb_, Xub = xub_, Residual = residual)
    ;111111111111111111111111111111111111
    scd = res[0]
    print, scd
    
    sza_ = sza[i, j]
    vza_ = vza[i, j]
    
    s_pos = where((sza_model - sza_) eq min(abs(sza_model - sza_)))
    v_pos = where((vza_model - vza_) eq min(abs(vza_model - vza_)))   
    
    amf_ = amf[v_pos, s_pos]
;    vcd[i, j] = scd/amf_
vcd[i, j] = scd
    
    print, i, j
  endfor
  ;print, i
endfor
;  openW, lun, path + "vcd.txt", /get_lun
;  printf, lun, vcd
;  free_lun, lun

;  vcd[where(vcd lt 0)] = -!values.f_nan
;;  vcd[where(vcd gt 10^17)] = -!values.f_nan
;  for i = 0, n_col - 1 do begin
;    vcd_ = vcd[i, *]
;    vcd[i, *] = diff_poly(indgen(n_row), vcd_)
;  endfor
  
;  temp = read_ascii("D:\OMI\out\ss.txt")
;  vcd = temp.(0)
;  vcd = vcd
;  vcd[where(~finite(vcd))] = - 10 *30

  path1 = "D:\omi\in\"
d_name = "OMI-Aura_L2-OMHCHO_2005m0921t0529-o06305_v003-2014m0620t065847.he5"
fid = h5f_open(path1 + d_name)
id = h5d_open(fid, "/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/AirMassFactor")
data = h5d_read(id)
  h5f_close, fid
  vcd = vcd/data
  vcd[where(~finite(vcd))] = - 10 *30
  
  a = to_tif(vcd, lon, lat, n_col, n_row, path + "a.tif")
print, "ok"
end
;--------------------------------------------------------------------------------------------------------


function cal_wavlen, wavcoef
  wc = wavcoef
  wavrefcol = 281
  wavran = 557
  ;wavlenref = np.linspace(1, wavran, wavran) - wavrefcol
  wavlenref = findgen(wavran) - wavrefcol
  
  wavlen = wc[0] + $
  wc[1] * wavlenref + $
  wc[2] * wavlenref^2 + $
  wc[3] * wavlenref^3 + $
  wc[4] * wavlenref^4
  ;if (wavlen[0] lt 300 or wavlen[-1] gt 390) then wavlen[*] = !values.f_nan
  return, wavlen
end

function diff_poly, x, y
  y_ = y[where(finite(y))]
  x_ = x[where(finite(y))]
  coef = poly_fit(x_, y_, 3)
  diff = y - poly(x, coef)
  return, diff
end



function to_tif, VCD_O3, O3Rad_Long, O3Rad_Lati, col, row, OutFileName
  index_east = where(O3Rad_Long gt 0 and VCD_O3 gt 0, count)  ;pick up the non-zero pixels in the eastern half
  east = max(O3Rad_Long(index_east))
  west = min(O3Rad_Long(index_east))
  north = max(O3Rad_Lati(index_east))
  south = min(O3Rad_Lati(index_east))

  ; ;;OMI 可见波段分辨率为13*24km(along and across track), 因此经度间隔设为0.2，纬度间隔设为0.125
  ; grid_east = round((east - 0.125) / 0.25)
  ; grid_west = round((west - 0.125) / 0.25) & long_leftup = grid_west * 0.25 + 0.125

  ; 定义栅格图像的边界范围，由中心像元个数差确定，定义每个像元的大小为0.125*0.125
  grid_east = round((east - 0.0625) / 0.125)
  grid_west = round((west - 0.0625) / 0.125) & long_leftup = grid_west * 0.125 + 0.0625
  grid_north = round((89.9375 - north) / 0.125) & lati_leftup = 89.9375 - grid_north * 0.125
  grid_south = round((89.9375 - south) / 0.125)
  col_tiff = grid_east - grid_west + 1
  row_tiff = grid_south - grid_north + 1

  ;VCD_GeoGrid = dblarr(1440, 1440) ; found the Lati-Long grid - 90N~90S, 0~180E
  VCD_GeoGrid = dblarr(col_tiff, row_tiff)
  count_GeoGrid = uintarr(col_tiff, row_tiff)
  num_west = 0  ; count of the pixels of the western longitude
  for i = 0, col - 1 do begin
    for j = 0, row - 1 do begin
      ;if O3Rad_Long[i,j] le 0 then num_west = num_west + 1

      if O3Rad_Long[i,j] gt 0 and VCD_O3[i,j] gt 0 then begin
        m = round((lati_leftup - O3Rad_Lati[i,j]) / 0.125)    ;obtain the row number of the current vcd pixel
        n = round((O3Rad_Long[i,j] - long_leftup) / 0.125)    ;obtain the column number of the current vcd pixel
        VCD_GeoGrid[n, m] = VCD_GeoGrid[n, m] + VCD_O3[i,j]
        count_GeoGrid[n, m] = count_GeoGrid[n, m] + 1
        ;       endif
      endif

    endfor
  endfor
  VCD_GeoGrid = temporary(VCD_GeoGrid) / temporary(count_GeoGrid)
  m = size(temporary(O3Rad_Lati)) & m = size(temporary(O3Rad_Long)) & m = size(temporary(VCD_O3_correct))

  ;############################################################################################
  ;interpolate the tiff image
  point_Z = FLTARR(1)
  point_X = FLTARR(1)
  point_Y = FLTARR(1)
  s=0

  index_valid = where(VCD_GeoGrid gt 0, num_valid)
  point_Y = ulong(index_valid) / ulong(col_tiff)      ;整除得到非0像元的行号
  point_X = ulong(index_valid) mod ulong(col_tiff)    ;整除后余数表示非0像元的列号
  point_Z = reform(VCD_GeoGrid[index_valid], 1, num_valid)    ;将非0像元重新排列为一列向量

  TRIANGULATE, point_X, point_Y, tr
  R = GRIDDATA(point_X, point_Y, point_Z,ANISOTROPY=[1,2,0],/DEGREES,DELTA=[1,1],DIMENSION=[col_tiff,row_tiff],START=[0,0],$
    TRIANGLES=tr,SEARCH_ELLIPSE=[30,2,150],MIN_POINTS=10,MISSING=0,SECTORS=1,SMOOTHING=0.5)

  VCD_GeoGrid = temporary(R)

  PIXELSCALE = [0.125 , 0.125, 0.00000000]
  ControlPoint = [0.00000000, 0.00000000, 0.00000000, long_leftup, lati_leftup, 0.00000000];像元左上角坐标

  GEOTIFF = { MODELPIXELSCALETAG:PIXELSCALE, MODELTIEPOINTTAG:ControlPoint, GTMODELTYPEGEOKEY:2, $
    GTRASTERTYPEGEOKEY:2, GTCITATIONGEOKEY:"Geographic", GEOGRAPHICTYPEGEOKEY:4326, $
    GEOGGEODETICDATUMGEOKEY:6326, GEOGANGULARUNITSGEOKEY:9102}
  WRITE_TIFF, OutFileName, VCD_GeoGrid / 1e15, geotiff=geotiff ,/FLOA; / 2.6867e16
  return, 1

end