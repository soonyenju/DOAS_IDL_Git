pro omps_hcho
  path = "D:\OMPS\in\"
  d_name = "OMPS-NPP-TC_SDR_EV_NASA-v1.0-2016m0226t041822-o00001-2016m0226t013855.h5"
  fid = h5f_open(path + d_name)
  id = h5d_open(fid, "/CALIBRATION_DATA/BinScheme1/SolarFlux")
  irr = h5d_read(id)
  
  id = h5d_open(fid, "/CALIBRATION_DATA/BinScheme1/SolarFluxWavelengths")
  wavlen = h5d_read(id)  
  
  id = h5d_open(fid, "/GEOLOCATION_DATA/BinScheme1/Longitude")
  lon = h5d_read(id)  
  
  id = h5d_open(fid, "/GEOLOCATION_DATA/BinScheme1/Latitude")
  lat = h5d_read(id)  
  
  id = h5d_open(fid, "/GEOLOCATION_DATA/BinScheme1/SolarZenithAngle")
  sza = h5d_read(id)
  
  id = h5d_open(fid, "/GEOLOCATION_DATA/BinScheme1/SatelliteZenithAngle")
  vza = h5d_read(id)    

  id = h5d_open(fid, "/SCIENCE_DATA/BinScheme1/Radiance")
  rad = h5d_read(id)

  h5f_close, fid
;--------------------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------------------
  hcho = make_array(2, file_lines(path + "hcho_clip.txt"))
  xsc = read_ascii(path + "hcho_clip.txt")
  hcho = xsc.(0)
  xsc = read_ascii(path + "bro_clip.txt")
  bro = xsc.(0)
  xsc = read_ascii(path + "oclo_clip.txt")
  oclo = xsc.(0)
  xsc = read_ascii(path + "o3_clip.txt")
  o3 = xsc.(0)
  xsc = read_ascii(path + "o4_clip.txt")
  o4 = xsc.(0)
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
  dim = size(rad, /dimension)
  n_ban = dim[0]
  n_col = dim[1]
  n_row = dim[2]
  vcd = double(make_array(n_col, n_row))
  for i = 0, n_col - 1 do begin
    for j = 0, n_row - 1 do begin
      rad_ = rad[*, i, j]
      irr_ = irr[*, i]
      wl = wavlen[*, i]
      ind = where(wl gt 327.7 and wl lt 356.0)
      wl = wl[ind]
      rad_ = rad_[ind]
      irr_ = irr_[ind]
      
      cosza = abs(cos(sza[i, j] /180 * !pi))
      ref = (!pi * rad_) / (cosza * irr_)
      ref[where(ref lt 0.0 or ref gt 1.0)] = -!values.f_nan
      ref_ = -alog(ref)
      pos = where(finite(ref_))
      s = size(pos)
      if (pos[0] eq -1 or s[1] lt 5) then begin
        vcd[i, j] = 0
        continue
      endif
      
      ref_ = ref_[pos]
      wl = wl[pos]
      
      y = diff_poly(wl, ref_)

      x1 = interpol(hcho[1, *], hcho[0, *], wl)
      x1 = diff_poly(wl, x1)
      x2 = interpol(bro[1, *], bro[0, *], wl)
      x2 = diff_poly(wl, x2)
      x3 = interpol(oclo[1, *], oclo[0, *], wl)
      x3 = diff_poly(wl, x3)
      x4 = interpol(o3[1, *], o3[0, *], wl)
      x4 = diff_poly(wl, x4)
      x5 = interpol(o4[1, *], o4[0, *], wl)
      x5 = diff_poly(wl, x5)

      x = [transpose(x1), transpose(x2), transpose(x3), transpose(x4)]
      res = regress(x, y, CHISQ = goodness, SIGMA=sigma, CONST=const)
      print, res
      scd = res[0]
      
      sza_ = sza[i, j]
      vza_ = vza[i, j]

      s_pos = where((sza_model - sza_) eq min(abs(sza_model - sza_)))
      v_pos = where((vza_model - vza_) eq min(abs(vza_model - vza_)))

      amf_ = amf[v_pos, s_pos]
      vcd[i, j] = scd/amf_      
      
      print, i, j
    endfor
  endfor

  a = to_tif(vcd, lon, lat, n_col, n_row, path + "a.tif")
  print, "ok"
end


function diff_poly, x, y
  coef = poly_fit(x, y, 3)
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