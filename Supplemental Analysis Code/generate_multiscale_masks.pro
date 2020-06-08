;*********************************************************************
;File: generate_multiscale_masks.pro
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

;requires the provision of original FETE output and masked versions of thresholded FETE output
pro generate_multiscale_masks, reference_file, output_dir, fete=fete
  compile_opt idl2
  
  envi_open_file, reference_file, r_fid=ref_fid
  if ref_fid eq -1 then return
  map_info = envi_get_map_info(fid=ref_fid)
  envi_file_query, ref_fid, ns=ns, nl=nl
  cornerx = [0.0d,ns,ns,0.0d]
  cornery = [0.0d,0.0d,nl,nl]
  envi_convert_file_coordinates, ref_fid, cornerx, cornery, mapx, mapy, /to_map
  xmin = min(mapx, max=xmax)
  ymin = min(mapy, max=ymax)
  extents = [xmin,ymin,xmax,ymax]
  te = strjoin(strtrim(string(extents,format='(D10.6)'),2),' ')

  xgsd = [map_info.ps[0],0.010d,0.050d,0.100d,0.500d]
  ygsd = [map_info.ps[1],0.010d,0.050d,0.100d,0.500d]  
  string_xgsd = strtrim(string(xgsd,format='(D12.10)'),2)
  string_ygsd = strtrim(string(ygsd,format='(D12.10)'),2)  
  file_gsd = ['005','010','050','100','500']
  num_gsd = n_elements(xgsd)
  
  if keyword_set(fete) then begin
    suffixes = ['.ramp_20.mask.tif','.ramp_10.mask.tif','.ramp_05.mask.tif','.ramp_01.mask.tif']
    data_files = envi_file_strip(reference_file,/back) + suffixes
    file_cutoff = ['20','10','05','01']
    num_files = n_elements(data_files)
  endif else begin
    num_files = 1
    data_files = reference_file
  endelse
  
  for j=0,num_files-1 do begin
    for i=0,num_gsd-1 do begin
      if keyword_set(fete) then begin
        out_file = output_dir + file_basename(reference_file, '.tif') + '.cutoff_' + file_cutoff[j] + '.mask.' + file_gsd[i] + '.tif'
        ;no need to reprocess the first layer, just make a copy of it
        if i eq 0 then begin
          file_copy, data_files[j], out_file, /overwrite
          continue
        endif      
        result = mask_aggregation(data_files[j], [xgsd[i],ygsd[i]], out_file)
      endif else begin
        out_file = output_dir + file_basename(reference_file, '.tif') + '.mask.' + file_gsd[i] + '.tif'
        if i eq 0 then begin
          file_copy, data_files[j], out_file, /overwrite
          continue
        endif
        result = mask_aggregation(data_files[j], [xgsd[i],ygsd[i]], out_file)        
      endelse
    endfor
  endfor
  
  envi_file_mng, id=ref_fid, /remove
  
end