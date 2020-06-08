;*********************************************************************
;File: run_stats_and_masks.pro
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

;requires original FETE output 
pro run_stats_and_masks, input_file
  compile_opt idl2
   
  site_file = "V:\Australia\First Arrival\SitesMorethan35ka_B_RatingandAbove_.evf"
  dem_file = "V:\Australia\First Arrival\sahul_merged_dem_3px_blend_mean_first_arrival_no_data.tif"
  dem_no_data = -32768
  land_dir = "V:\Australia\First Arrival\Masks\"  

  suffixes = ['.ramp_20.mask.distance.tif','.ramp_10.mask.distance.tif','.ramp_05.mask.distance.tif','.ramp_01.mask.distance.tif']
  data_files = envi_file_strip(input_file,/back) + suffixes
  prob_output_dir = envi_file_strip(input_file, /path) + path_sep() + 'Statistics' + path_sep()
  if ~file_test(prob_output_dir, /directory) then file_mkdir, prob_output_dir

  mask_output_dir = envi_file_strip(input_file, /path) + path_sep() + 'Masks' + path_sep()
  if ~file_test(mask_output_dir, /directory) then file_mkdir, mask_output_dir

  radius = 1000.0d
  increment = 5000.0d
  nsteps = 20
  ntrials = 10000

  generate_statistics, site_file, dem_file, dem_no_data, data_files, radius, increment, nsteps, ntrials, prob_output_dir
  
  generate_multiscale_masks, input_file, mask_output_dir, /fete
  
  generate_mask_stats, land_dir, mask_output_dir
   
end