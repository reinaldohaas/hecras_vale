wkt = 'PROJCS[" SIRGAS_2000_UTM_Zone_22S\,GEOGCS[\GCS_SIRGAS_2000\,DATUM[\D_SIRGAS_2000\,SPHEROID[\GRS_1980\,6378137.0,298.257222101]],PRIMEM[\Greenwich\,0.0],UNIT[\Degree\,0.0174532925199433]],PROJECTION[\Transverse_Mercator\],PARAMETER[\False_Easting\,500000.0],PARAMETER[\False_Northing\,10000000.0],PARAMETER[\Central_Meridian\,-51.0],PARAMETER[\Scale_Factor\,0.9996],PARAMETER[\Latitude_Of_Origin\,0.0],UNIT[\Meter\,1.0]]'
for p in ['_anti/taha_ai/taha_ai.g01', 'taha_ai.g01']:
 try:
 with open(p, 'r', encoding='latin-1') as f:
 lines = f.readlines()
 for i in range(min(10, len(lines))):
 if lines[i].startswith('Spatial Reference System='):
 lines[i] = 'Spatial Reference System=' + wkt + '\n'
 with open(p, 'w', encoding='latin-1') as f:
 f.writelines(lines)
 except:
 pass
with open('_anti/taha_ai/SIRGAS2000_UTM22S.prj', 'w', encoding='utf-8') as f:
 f.write(wkt)
print('CRS atualizado!')