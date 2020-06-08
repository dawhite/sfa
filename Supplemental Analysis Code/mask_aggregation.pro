;*********************************************************************
;File: mask_aggregation.pro
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

function mask_aggregation, in_file, gsd, out_file
  compile_opt idl2
  
  if n_elements(gsd) ne 2 then return, 0
  
  envi_open_file, in_file, r_fid=in_fid
  if in_fid eq -1 then return, 0
  
  in_map_info = envi_get_map_info(fid=in_fid)
  ;bail if the requested GSD is smaller than the original
  if (in_map_info.ps[0] gt gsd[0]) or (in_map_info.ps[1] gt gsd[1]) then return, 0
  
  ratio = gsd / in_map_info.ps
  
  envi_file_query, in_fid, ns=ns, nl=nl, dims=dims
  
  cornerx = [0.0d,ns,ns,0.0d]
  cornery = [0.0d,0.0d,nl,nl]
  envi_convert_file_coordinates, in_fid, cornerx, cornery, mapx, mapy, /to_map
  minx = min(mapx, max=maxx)
  miny = min(mapy, max=maxy)
  
  out_xsize = ceil((mapx[1] - mapx[0])/gsd[0])
  out_ysize = ceil((mapy[0] - mapy[2])/gsd[1])
  out_array = bytarr(out_xsize, out_ysize)
  
  in_data = envi_get_data(fid=in_fid, dims=dims, pos=0)
  
  ;go cell by cell in the output array and fill with corresponding input data
  for j=0,out_ysize-1 do begin
    ystart = floor(j*ratio[1]) > 0
    yend = ceil((j+1)*ratio[1]) < (nl-1)
    for i=0,out_xsize-1 do begin
      xstart = floor(i*ratio[0]) > 0
      xend = ceil((i+1)*ratio[0]) < (ns-1)
      ;if any input cell is valid, mark the output cell as valid (could be more sophisticated, like a weighting)
      out_array[i,j] = total(in_data[xstart:xend,ystart:yend]) gt 0      
    endfor
  endfor
  
  out_map_info = envi_map_info_create(mc=[0.0d,0.0d,mapx[0],mapy[0]], ps=gsd, proj=in_map_info.proj)
  
  envi_enter_data, out_array, map_info=out_map_info, r_fid=out_fid
  envi_file_query, out_fid, dims=out_dims
  envi_output_to_external_format, fid=out_fid, dims=out_dims, pos=0, /tiff, out_name=out_file
  envi_file_mng, id=out_fid, /remove
  
  file_delete, envi_file_strip(out_file,/back) + '.tfw'
  
  envi_file_mng, id=in_fid, /remove
  
  return, 1
  
end