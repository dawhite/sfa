;*********************************************************************
;File: calculate_wbica.pro
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

;weighted approximation of the Bayesian Information Criteria
pro calculate_wBICa
  compile_opt idl2
  
  num_samples = 40
  weighted_params = 4
  
  ;0 = no rivers
  ;1 = Strahler Order 8 and higher rivers
  ;2 = Strahler Order 9 and higher rivers
  include_river = 2
  
  ;total number of origin and destination pairs per model
  num_od_pairs = dblarr(6)
  num_od_pairs[*] = [$
    1473024400,$
    3620462160,$
    10639328952,$
    8898526224,$
    4328726680,$
    12720681796$    
  ]
  max_od_pairs = max(num_od_pairs)

  base_dir = 'V:\Australia\First Arrival\'

  land_mask = 'V:\Australia\First Arrival\Masks\first_arrival_mask.mask.100.tif'

  fete_dirs = [$
    'V:\Australia\First Arrival\female_unweighted_grid_50\', $
    'V:\Australia\First Arrival\female_unweighted_coastline_to_grid_50\', $
    'V:\Australia\First Arrival\female_unweighted_coastline_to_water\', $
    'V:\Australia\First Arrival\female_unweighted_coastline_to_coastline\', $
    'V:\Australia\First Arrival\female_unweighted_water_to_grid_50\', $
    'V:\Australia\First Arrival\female_unweighted_water_to_water\', $
    'V:\Australia\First Arrival\female_weighted_grid_50\', $
    'V:\Australia\First Arrival\female_weighted_coastline_to_grid_50\', $
    'V:\Australia\First Arrival\female_weighted_coastline_to_water\', $
    'V:\Australia\First Arrival\female_weighted_coastline_to_coastline\', $
    'V:\Australia\First Arrival\female_weighted_water_to_grid_50\', $
    'V:\Australia\First Arrival\female_weighted_water_to_water\' $
    ]

  if include_river gt 0 then begin
    if include_river eq 1 then river_dirs = [$
      'V:\Australia\First Arrival\female_river_weighted_grid_50\', $
      'V:\Australia\First Arrival\female_river_weighted_coastline_to_grid_50\', $
      'V:\Australia\First Arrival\female_river_weighted_coastline_to_water\', $
      'V:\Australia\First Arrival\female_river_weighted_coastline_to_coastline\', $
      'V:\Australia\First Arrival\female_river_weighted_water_to_grid_50\', $
      'V:\Australia\First Arrival\female_river_weighted_water_to_water\' $
      ] else if include_river eq 2 then river_dirs = [$
      'V:\Australia\First Arrival\female_river_weighted_grid_50\', $
      'V:\Australia\First Arrival\female_river_weighted_coastline_to_grid_50\', $
      'V:\Australia\First Arrival\female_river_weighted_coastline_to_water\', $
      'V:\Australia\First Arrival\female_river_weighted_coastline_to_coastline5\', $
      'V:\Australia\First Arrival\female_river_weighted_water_to_grid_50\', $
      'V:\Australia\First Arrival\female_river_weighted_water_to_water\' $
      ]
    use_dirs = [fete_dirs,river_dirs]
    ncols = 3
  endif else begin
    use_dirs = fete_dirs
    ncols = 2
  endelse  

  num_fete = n_elements(use_dirs)

  ;initialize the arrays we'll need
  probabilities = replicate(-1.0d,4*ncols,6)
  land_proportions = probabilities
  od_weights = probabilities
  raw_networks = probabilities
  rss = probabilities
  BICa_prob = probabilities
  wBICa_prob = probabilities
  BICa_rss = probabilities
  wBICa_rss = probabilities
  BICa_combined = probabilities
  wBICa_combined = probabilities  
  
  grids = replicate('none',4*ncols,6)
  
  for i=0,num_fete-1 do begin
    
    if ~file_test(use_dirs[i],/directory) then continue
    
    cur_stats_dir = use_dirs[i] + 'Statistics'
    summary_files = file_search(cur_stats_dir, '*.summary.csv', count=summary_count)
    distance_files = file_search(cur_stats_dir, '*.site_distance.csv', count=distance_count)
    cur_mask_dir = use_dirs[i] + 'Masks'
    proportions_files = file_search(cur_mask_dir, '*.proportions.csv', count=proportions_count)
    raw_files = file_search(cur_mask_dir, '*.raw_counts.csv', count=raw_count)
    mask_files = file_search(cur_mask_dir, '*.mask.100.tif', count=mask_count)
    
    cur_column = i / 6
    cur_row = i mod 6
    
    ;reverse the file orders so highest cutoff is first
    summary_files = reverse(summary_files)
    distance_files = reverse(distance_files)
    mask_files = reverse(mask_files)
    
    for j=0,summary_count-1 do begin
      cur_lines = file_lines(summary_files[j])
      openr, lun, summary_files[j], /get_lun
      summary_lines = strarr(cur_lines)
      readf, lun, summary_lines
      free_lun, lun
      summary_split = strsplit(summary_lines[1], ',', /extract)
      probabilities[j*ncols+cur_column,cur_row] = double(summary_split[22])
      ;take care of the OD weights while we're in here
      od_weights[j*ncols+cur_column,cur_row] = max_od_pairs/num_od_pairs[cur_row]
    endfor

    for j=0,distance_count-1 do begin
      cur_lines = file_lines(distance_files[j])
      openr, lun, distance_files[j], /get_lun
      distance_lines = strarr(cur_lines)
      readf, lun, distance_lines
      free_lun, lun
      ;convert meters to kilometers
      cur_rss = total((double(distance_lines)/1000.0d)^2)
      rss[j*ncols+cur_column,cur_row] = cur_rss
    endfor
    
    for j=0,mask_count-1 do grids[j*ncols+cur_column,cur_row] = mask_files[j]
      
    cur_lines = file_lines(proportions_files[0])
    openr, lun, proportions_files[0], /get_lun
    proportions_lines = strarr(cur_lines)
    readf, lun, proportions_lines
    free_lun, lun
    ;reverse line order so highest cutoff is first
    proportions_lines = reverse(proportions_lines[1:4])
    for j=0,3 do begin
      proportions_split = strsplit(proportions_lines[j], ',', /extract)
      land_proportions[j*ncols+cur_column,cur_row] = double(proportions_split[1])
    endfor
    
    cur_lines = file_lines(raw_files[0])
    openr, lun, raw_files[0], /get_lun
    raw_lines = strarr(cur_lines)
    readf, lun, raw_lines
    free_lun, lun
    ;reverse line order so highest cutoff is first
    raw_lines = reverse(raw_lines[1:4])
    for j=0,3 do begin
      raw_split = strsplit(raw_lines[j], ',', /extract)
      raw_networks[j*ncols+cur_column,cur_row] = double(raw_split[1])
    endfor    
    
  endfor
  
  ;create an epsilon value to avoid taking the natural log of 0.0
  epsilon = 1.0d-16
  
  print, probabilities
  print, ' '
  print, land_proportions
  print, ' '
  print, raw_networks
  print, ' '
  print, rss
  print, ' '
  
  ;keep everything with p=0.100 or better
  where_good = where(probabilities ge 0.9d, good_count)
  
  ;normalize the RSS
  norm_rss = replicate(-1.0d,4*ncols,6)
  norm_rss[where_good] = (rss[where_good] / max(rss[where_good]))
  
  max_raw = max(raw_networks[where_good])
  
  ;For BICa, we are penalizing models with larger networks, so invert the network counts
  BICa_prob[where_good] = alog(num_samples) * weighted_params - 2.0d * alog(probabilities[where_good] * (max_raw/raw_networks[where_good]) * od_weights[where_good] + epsilon)
  
  print, BICa_prob

  print, ' '
  
  min_BICa_prob = min(BICa_prob[where_good])
  delta_from_min = BICa_prob[where_good] - min_BICa_prob
  denominator = total(exp(-delta_from_min/2.0d))
  wBICa_prob[where_good] = exp(-delta_from_min/2.0d)/denominator
  
  print, wBICa_prob
  
  ;output results
  if include_river gt 0 then begin
    if include_river eq 1 then begin
      out_probabilities = base_dir + 'probabilities_river_8.csv'
      out_proportions = base_dir + 'land_proportions_river_8.csv'
      out_networks = base_dir + 'network_cell_counts_river_8.csv'
      out_rss = base_dir + 'rss_river_8.csv'
      out_BICa = base_dir + 'BICa_river_8.csv'
      out_wBICa = base_dir + 'wBICa_river_8.csv'
    endif else if include_river eq 2 then begin
      out_probabilities = base_dir + 'probabilities_river_9.csv'
      out_proportions = base_dir + 'land_proportions_river_9.csv'
      out_networks = base_dir + 'network_cell_counts_river_9.csv'
      out_rss = base_dir + 'rss_river_9.csv'
      out_BICa = base_dir + 'BICa_river_9.csv'
      out_wBICa = base_dir + 'wBICa_river_9.csv'      
    endif
  endif else begin
    out_probabilities = base_dir + 'probabilities_no_river.csv'
    out_proportions = base_dir + 'land_proportions_no_river.csv'
    out_networks = base_dir + 'network_cell_counts_no_river.csv'
    out_rss = base_dir + 'rss_no_river.csv' 
    out_BICa = base_dir + 'BICa_no_river.csv'
    out_wBICa = base_dir + 'wBICa_no_river.csv'       
  endelse
  
  openw, lun_probabilities, out_probabilities, /get_lun
  openw, lun_proportions, out_proportions, /get_lun
  openw, lun_networks, out_networks, /get_lun
  openw, lun_rss, out_rss, /get_lun
  openw, lun_BICa, out_BICa, /get_lun
  openw, lun_wBICa, out_wBICa, /get_lun
  
  for i=0,5 do begin
    out_line = strjoin(string(probabilities[*,i],format='(D25.16)'),',')
    printf, lun_probabilities, out_line

    out_line = strjoin(string(land_proportions[*,i],format='(D25.16)'),',')
    printf, lun_proportions, out_line
    
    out_line = strjoin(string(raw_networks[*,i],format='(D25.16)'),',')
    printf, lun_networks, out_line
    
    out_line = strjoin(string(rss[*,i],format='(D25.16)'),',')
    printf, lun_rss, out_line    
        
    out_line = strjoin(string(BICa_prob[*,i],format='(D25.16)'),',')
    printf, lun_BICa, out_line

    out_line = strjoin(string(wBICa_prob[*,i],format='(D25.16)'),',')
    printf, lun_wBICa, out_line  
  endfor
  
  free_lun, lun_probabilities
  free_lun, lun_proportions
  free_lun, lun_networks
  free_lun, lun_rss
  free_lun, lun_BICa
  free_lun, lun_wBICa
  
  envi_open_file, land_mask, r_fid=land_fid
  envi_file_query, land_fid, ns=ns, nl=nl, dims=dims
  land_data = envi_get_data(fid=land_fid, dims=dims, pos=0)
  map_info = envi_get_map_info(fid=land_fid)
  out_map_info = envi_map_info_create(proj=map_info.proj, ps=map_info.ps, mc=map_info.mc)
  out_grid = dblarr(ns,nl)  
  
  for i=0,good_count-1 do begin
    envi_open_file, grids[where_good[i]], r_fid=grid_fid
    grid_data = envi_get_data(fid=grid_fid, dims=dims, pos=0)
    out_grid += (double(grid_data) * wBICa_prob[where_good[i]])
    envi_file_mng, id=grid_fid, /remove
  endfor
  
  envi_file_mng, id=land_fid, /remove
  
  ;set all non-land cells to -1
  where_not_land = where(land_data eq 0)
  out_grid[where_not_land] = -1.0d
  
  envi_enter_data, out_grid, map_info=out_map_info, r_fid=out_fid
  if include_river gt 0 then begin
    if include_river eq 1 then out_file = base_dir + 'wBICa_grid_composite_river_8.tif' else $
      if include_river eq 2 then out_file = base_dir + 'wBICa_grid_composite_river_9.tif'
  endif else begin
    out_file = base_dir + 'wBICa_grid_composite_no_river.tif'    
  endelse

  if file_test(out_file) then file_delete, out_file, /quiet, /allow_nonexistent
  envi_output_to_external_format, fid=out_fid, dims=dims, pos=0, out_name=out_file, /tiff
  envi_file_mng, id=out_fid, /remove
  file_delete, envi_file_strip(out_file, /back) + '.tfw'
    
end