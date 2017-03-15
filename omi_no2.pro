pro OMI_NO2


  NO2GFiles = DIALOG_PICKFILE(/READ, FILTER = ['*.he5'],/MULTIPLE_FILES)
  numfile = N_elements(NO2GFiles)



  filepath = file_DirName(NO2GFiles[0])



  for l=0,numfile-1 do begin

    inFileName = NO2GFiles[l]
    FileBaseName = file_basename(inFileName)
    FileDate = strmid(FileBaseName, 18, 14)
    savefile = filepath+'\'+filedate+'NO2.tif'
    print,FileDate
    File_IDN = H5F_OPEN(inFileName)
    R1_IDN=H5G_OPEN(File_IDN,  'HDFEOS')
    R2_IDN=H5G_OPEN(R1_IDN,  'SWATHS')
    R3_IDN=H5G_OPEN(R2_IDN,  'ColumnAmountNO2')
    R4_IDN=H5G_OPEN(R3_IDN,  'Data Fields')
    Dataset_IDN = H5D_OPEN(R4_IDN, 'ColumnAmountNO2Trop')
   ;Dataset_IDN = H5D_OPEN(R4_IDN, 'ColumnAmountNO2')
    DataNO2 = H5D_READ(Dataset_IDN)

    DIMSN=SIZE(DataNO2,/DIMENSIONS)
    col=DIMSN[0]
    row=DIMSN[1]

    ID_invalidValue = H5A_OPEN_NAME(Dataset_IDN, '_FillValue')
    invalidValue = H5A_READ(ID_invalidValue)
    ;print,invalidValue
    ID_missedValue = H5A_OPEN_NAME(Dataset_IDN, 'MissingValue')
    missedValue = H5A_READ(ID_missedValue)
    ;print,missedValue
    H5D_CLOSE, Dataset_IDN


;    SCD_IDN = H5D_OPEN(R4_IDN, 'SlantColumnAmountNO2')
;    DataNO2 = H5D_READ(SCD_IDN)
;
;    DIMSN=SIZE(DATANO2,/DIMENSIONS)
;    col=DIMSN[0]
;    row=DIMSN[1]
;
;    ID_invalidValue = H5A_OPEN_NAME(SCD_IDN, '_FillValue')
;    invalidValue = H5A_READ(ID_invalidValue)
;    ;print,invalidValue
;    ID_missedValue = H5A_OPEN_NAME(SCD_IDN, 'MissingValue')
;    missedValue = H5A_READ(ID_missedValue)
;    ;print,missedValue
;    H5D_CLOSE,SCD_IDN


    index_0 = where(DataNO2 lt 0, count)
    if count gt 0 then begin
      DataNO2[index_0] = 0
    endif
    ;print,DATANO2
    ;Save dat Files
    ; OPENW, 1, inSaveName
    ; WRITEU, 1, DataNO2
    ; CLOSE, 1

    ;print,col,row

    R5_IDN=H5G_OPEN(R3_IDN,  'Geolocation Fields')
    Lati_IDN = H5D_OPEN(R5_IDN, 'Latitude')
    Long_IDN = H5D_OPEN(R5_IDN, 'Longitude')
    LatiNO2 = H5D_READ(Lati_IDN)
    LongNO2 = H5D_READ(Long_IDN)
    H5D_CLOSE, Lati_IDN
    H5D_CLOSE, Long_IDN
    H5F_CLOSE, File_IDN

;    east = 116
;    west = 112
;    north = 24
;    south = 21
      index_east = where(LongNO2 gt 0 and DataNO2 gt 0 and LatiNO2 gt 0, count)  ;pick up the non-zero pixels in the eastern half
      east = max(LongNO2(index_east))
      west = min(LongNO2(index_east))
      north = max(LatiNO2(index_east))
      south = min(LatiNO2(index_east))
    ;

    grid_east = round((east - 0.0625) / 0.125)
    grid_west = round((west - 0.0625) / 0.125) & long_leftup = grid_west * 0.125 + 0.0625
    grid_north = round((89.9375 - north) / 0.125) & lati_leftup = 89.9375 - grid_north * 0.125
    grid_south = round((89.9375 - south) / 0.125)
    col_tiff = grid_east - grid_west + 1
    row_tiff = grid_south - grid_north + 1

    Data_GeoGrid = dblarr(col_tiff, row_tiff)
    count_GeoGrid = uintarr(col_tiff, row_tiff)
    num_west = 0  ; count of the pixels of the western longitude
    for i = 0, col - 1 do begin
      for j = 0, row - 1 do begin
        if LongNO2[i,j] gt west and LongNO2[i,j] lt east and LatiNO2[i,j] gt south and LatiNO2[i,j] lt north and DataNO2[i,j] gt 0 then begin
          m = round((lati_leftup - LatiNO2[i,j]) / 0.125)   ;obtain the row number of the current vcd pixel
          n = round((LongNO2[i,j] - long_leftup) / 0.125)   ;obtain the column number of the current vcd pixel
          Data_GeoGrid[n, m] = Data_GeoGrid[n, m] + DataNO2[i,j]
          count_GeoGrid[n, m] = count_GeoGrid[n, m] + 1
        endif
      endfor
    endfor

    Data_GeoGrid = temporary(Data_GeoGrid) / temporary(count_GeoGrid)
    m = size(temporary(LatiNO2)) & m = size(temporary(LongNO2)) & m = size(temporary(DataNO2))


    ;interpolate the tiff image
    point_Z = FLTARR(1)
    point_X = FLTARR(1)
    point_Y = FLTARR(1)
    s=0


    index_valid = where(Data_GeoGrid gt 0, num_valid)
    if num_valid gt 200 then begin
      point_Y = ulong(index_valid) / ulong(col_tiff)
      point_X = ulong(index_valid) mod ulong(col_tiff)
      point_Z = reform(Data_GeoGrid[index_valid], 1, num_valid)

      TRIANGULATE, point_X, point_Y, tr
      R = GRIDDATA(point_X, point_Y, point_Z,ANISOTROPY=[1,2,0],/DEGREES,DELTA=[1,1],DIMENSION=[col_tiff,row_tiff],START=[0,0],$
        TRIANGLES=tr,SEARCH_ELLIPSE=[30,2,150],MIN_POINTS=10,MISSING=0,SECTORS=1,SMOOTHING=0.5)

      Data_GeoGrid = temporary(R)


      PIXELSCALE = [0.125, 0.125, 0.00000000]
      ControlPoint = [0.00000000, 0.00000000, 0.00000000, long_leftup, lati_leftup, 0.00000000]

      GEOTIFF = { MODELPIXELSCALETAG:PIXELSCALE, MODELTIEPOINTTAG:ControlPoint, GTMODELTYPEGEOKEY:2, $
        GTRASTERTYPEGEOKEY:2, GTCITATIONGEOKEY:"Geographic", GEOGRAPHICTYPEGEOKEY:4326, $
        GEOGGEODETICDATUMGEOKEY:6326, GEOGANGULARUNITSGEOKEY:9102}
      WRITE_TIFF, savefile, Data_GeoGrid / 1e15, geotiff=geotiff ,/FLOAT
      ;print,savefile
    endif
  endfor
  rr = Dialog_message('Mission Complete', /INFORMATION , TITLE = 'NO2')

END