;*********************************************************************
;File: generate_mask_stats.pro
;Date: 28 Apr 2020
;Author: Devin A. White, PhD
;Email: Devin.White@gmail.com
;*********************************************************************
;Copyright (C) 2020  Devin A. White
;
;This program is free software: you can redistribute it and/or modify
;it under the terms of the GNU General Public License as published by
;the Free Software Foundation, either version 3 of the License, or
;(at your option) any later version.
;
;This program is distributed in the hope that it will be useful,
;but WITHOUT ANY WARRANTY; without even the implied warranty of
;MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;GNU General Public License for more details.
;
;You should have received a copy of the GNU General Public License
;along with this program.  If not, see https://www.gnu.org/licenses/.
;*********************************************************************
;Do not remove this header when reusing or repurposing this code
;*********************************************************************

pro generate_mask_stats, land_dir, path_dir
  compile_opt idl2
  
  file_gsd = ['005','010','050','100','500']
  real_gsd = [0.005d, 0.010d, 0.050d, 0.100d, 0.500d]
  num_gsd = n_elements(file_gsd)
  
  cutoffs = ['01','05','10','20']
  real_cutoffs = [0.01d, 0.05d, 0.10d, 0.20d]
  raw = ulon64arr(5,4)
  proportions = dblarr(5,4)
  
  for i=0,num_gsd-1 do begin
    search_string = '*.mask.'+file_gsd[i]+'.tif'
    
    land_masks = file_search(land_dir, search_string, count=land_count)
    if land_count eq 0 then continue

    path_masks = file_search(path_dir, search_string, count=path_count)
    if path_count eq 0 then continue
    
    envi_open_file, land_masks[0], r_fid=land_fid
    envi_file_query, land_fid, dims=land_dims
    land_data = envi_get_data(fid=land_fid, dims=land_dims, pos=0)
    land_total = total(land_data, /double)
            
    for j=0,path_count-1 do begin
      if (i eq 0 and j eq 0) then $
        base_split = strsplit(file_basename(path_masks[j]),'.',/extract)
      envi_open_file, path_masks[j], r_fid=path_fid
      envi_file_query, path_fid, dims=path_dims
      path_data = envi_get_data(fid=path_fid, dims=path_dims, pos=0)
      path_total = total(path_data, /double)
      envi_file_mng, id=path_fid, /remove
      raw[i,j] = path_total
      proportions[i,j] = path_total/land_total
    endfor

    envi_file_mng, id=land_fid, /remove

  endfor
  
  out_table = strarr(6,5)
  out_table[0,0] = '   '
  out_table[1:5,0] = strtrim(string(real_gsd,format='(D5.3)'),2)
  out_table[0,1:4] = strtrim(string(real_cutoffs,format='(D4.2)'),2)
  out_table[1:5,1:4] = strtrim(string(proportions,format='(D12.10)'),2)
  out_file = path_dir + base_split[0] + '.proportions.csv'
  openw, lun, out_file, /get_lun
  for i=0,4 do printf, lun, strjoin(out_table[*,i],',')
  free_lun, lun
  
  out_table[1:5,1:4] = strtrim(string(raw,format='(I30)'),2)
  out_file = path_dir + base_split[0] + '.raw_counts.csv'
  openw, lun, out_file, /get_lun
  for i=0,4 do printf, lun, strjoin(out_table[*,i],',')
  free_lun, lun
  
end