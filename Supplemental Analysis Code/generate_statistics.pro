;*********************************************************************
;File: generate_statistics.pro
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

;Calculates the Kolmogorov distribution function,
;which gives the probability that Kolmogorov's test statistic will exceed
;the value z assuming the null hypothesis. This gives a very powerful
;test for comparing two one-dimensional distributions.
;see, for example, Eadie et al, "statistocal Methods in Experimental
;Physics', pp 269-270).
;
;This function returns the confidence level for the null hypothesis, where:
;  z = dn*sqrt(n), and
;  dn  is the maximum deviation between a hypothetical distribution
;      function and an experimental distribution with
;  n    events
;
;NOTE: To compare two experimental distributions with m and n events,
;      use z = sqrt(m*n/(m+n))*dn
;
;Accuracy: The function is far too accurate for any imaginable application.
;          Probabilities less than 10^-15 are returned as zero.
;          However, remember that the formula is only valid for "large" n.
;Theta function inversion formula is used for z <= 1
;
;This function was translated by Rene Brun from PROBKL in CERNLIB.
;
;Ported from CERN TMath C++ class (TMath::KolmogorovProb):
;http://root.cern.ch/root/html/src/TMath.cxx.html

function ks_prob, z
  compile_opt idl2

  fj = [-2.0d,-8.0d,-18.0d,-32.0d]
  r = dblarr(4)
  w = 2.50662827d
  ;c1 = -pi**2/8, c2 = 9*c1, c3 = 25*c1
  c1 = -1.2337005501361697d
  c2 = -11.103304951225528d
  c3 = -30.842513753404244d

  u = abs(double(z))
  p = 0.0d

  if (u lt 0.2) then begin
    p = 1
  endif else if (u lt 0.755) then begin
    v = 1.0d / (u*u)
    p = 1.0d - w*(exp(c1*v) + exp(c2*v) + exp(c3*v))/u
  endif else if (u lt 6.8116) then begin
    r[1] = 0.0d
    r[2] = 0.0d
    r[3] = 0.0d
    v = u*u
    maxj = fix(1 > round(3.0d/u))
    for j=0,maxj-1 do r[j] = exp(fj[j]*v)
    p = 2.0d * (r[0] - r[1] + r[2] - r[3])
  endif else begin
    p = 0.0d
  endelse

  return, p

end


;The two-sample K-S test checks whether the two data samples come
;from the same distribution. This does not specify what that common
;distribution is (e.g. normal or not normal)
;
;Ported from CERN TMath C++ class (TMath::KolmogorovTest):
;http://root.cern.ch/root/html/src/TMath.cxx.html
;Probably not very efficient for IDL, but it should work
;
;return values:
;-1 -> failure
; 0 -> null hypothesis not rejected (distributions are the same)
; 1 -> null hypothesis rejected (distributions are different)
;function ks_two_sample, samp1, samp2, alpha_choice, distance=distance

;return value is now the p-value for the statistic!

function ks_two_sample, samp1, samp2, distance=distance
  compile_opt idl2

  ;critical values estimated using this implementation:
  ;http://www.soest.hawaii.edu/wessel/courses/gg313/Critical_KS.pdf
  ;for alpha = 0.10,0.05,0.025,0.01,0.005,0.001
  ;matches values from Wikipedia:
  ;https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test

  num_s1 = n_elements(samp1)
  na = double(num_s1)
  num_s2 = n_elements(samp2)
  nb = double(num_s2)

  a = double(samp1[sort(samp1)])
  b = double(samp2[sort(samp2)])

  ia = 0ULL
  ib = 0ULL
  ;cdf fraction for first sample
  sa = 1.0d/na
  ;cdf fraction for second sample
  sb = 1.0d/nb
  ;running difference
  rdiff = 0.0d
  ;maximum difference (supremum)
  rdmax = 0.0d
  ok = 0

  for i=0,num_s1+num_s2-1 do begin
    if (a[ia] lt b[ib]) then begin
      rdiff -= sa
      ia++
      if (ia ge num_s1) then begin
        ok = 1
        break
      endif
    endif else if (a[ia] gt b[ib]) then begin
      rdiff += sb
      ib++
      if (ib ge num_s2) then begin
        ok = 1
        break
      endif
    endif else begin
      x = a[ia]
      while (a[ia] eq x) do begin
        rdiff -= sa
        ia++
        if (ia ge num_s1) then break
      endwhile
      while (b[ib] eq x) do begin
        rdiff += sb
        ib++
        if (ib ge num_s2) then break
      endwhile
      if (ia ge num_s1) then begin
        ok = 1
        break
      endif
      if (ib ge num_s2) then begin
        ok = 1
        break
      endif
    endelse
    rdmax = rdmax > abs(rdiff)
  endfor

  if ok then begin
    distance = rdmax
    z = rdmax * sqrt((na*nb)/(na+nb))
    prob = ks_prob(z)
    return, prob
  endif else begin
    distance = -1.0d
    return, -1.0d
  endelse

end


;As described here:
;https://en.wikipedia.org/wiki/Spatial_descriptive_statistics
;https://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/h-how-multi-distance-spatial-cluster-analysis-ripl.htm
;https://www.seas.upenn.edu/~ese502/NOTEBOOK/Part_I/4_K_Functions.pdf
;In all three cases, one of the sample loops has only one member, i.e., no loop, because we're comparing point locations to a travel network represented as a single point, 
;where distance to the points is represented by the 3D travel distance rasters. The question we're asking is how much are those points "clustered" around the network at a given threshold.
function ripley_l, distances, search_radius, area
  compile_opt idl2
  
  dnels = double(n_elements(distances))
  ;original version from Wikipedia
  ;lambda = dnels / area
  ;return, sqrt(((1.0d / lambda) * total((distances lt search_radius) / dnels)) / !dpi)
  
  ;ESRI version, which they state is commonly used
  return,  sqrt((area * total(distances lt search_radius, /double)) / (!dpi * dnels * (dnels - 1)))
  
end


pro generate_statistics, site_file, dem_file, dem_no_data, data_files, radius, increment, nsteps, ntrials, prob_output_dir
  compile_opt idl2

  ;determine envelope area, in square meters, for the region of analysis
  envi_open_file, dem_file, r_fid=dem_fid
  envi_file_query, dem_fid, ns=dem_ns, nl=dem_nl, dims=dem_dims
  dem_data = envi_get_data(fid=dem_fid, dims=dem_dims, pos=0)
  dem_map_info = envi_get_map_info(fid=dem_fid)
  geo_proj = envi_proj_create(/geographic)
  corner_imgx = [0.0d,dem_ns,dem_ns,0.0d]
  corner_imgy = [0.0d,0.0d,dem_nl,dem_nl]
  envi_convert_file_coordinates, dem_fid, corner_imgx, corner_imgy, corner_mapx, corner_mapy, /to_map
  envi_convert_projection_coordinates, corner_mapx, corner_mapy, dem_map_info.proj, corner_lon, corner_lat, geo_proj
  min_lon = min(corner_lon, max=max_lon)
  min_lat = min(corner_lat, max=max_lat)  
  center_lat = (min_lat + max_lat)/2.0
  utm_proj = envi_proj_create(/utm, datum='WGS-84', zone=fix(31.0 + ((min_lon + max_lon)/2)/6.0), south=min_lat lt 0)
  envi_proj_units_convert_pixel_size, dem_map_info.ps, meter_ps, $
    i_units=dem_map_info.proj.units, i_proj=dem_map_info.proj, o_proj=utm_proj, ref_lat=center_lat
  envelope_area = total(dem_data ne dem_no_data) * meter_ps[0] * meter_ps[1]

  data_titles = ['FETE 20','FETE 10','FETE 05','FETE 01']
  file_titles = ['cutoff_20','cutoff_10','cutoff_05','cutoff_01']
   
  ndata = n_elements(data_titles)
    
  evf_id = envi_evf_open(site_file)
  envi_evf_info, evf_id, num_recs=num_recs, projection=site_proj
  
  xvals = dblarr(num_recs)
  yvals = dblarr(num_recs)
  
  for i=0,num_recs-1 do begin
    record = envi_evf_read_record(evf_id, i)
    xvals[i] = record[0]
    yvals[i] = record[1]
  endfor
  
  envi_evf_close, evf_id

  alphas = ['0.001','0.005','0.01','0.025','0.05','0.10']
  cutoffs = [0.001d,0.005d,0.01d,0.025d,0.05d,0.10d]
  zscores = [3.0902d,2.5758d,2.3263d,1.9600d,1.6449d,1.2816d]
  ncutoffs = n_elements(cutoffs)
  
  ;search out number of steps from supplied increment
  search_radius = (dindgen(nsteps) + 1) * increment
  
  for a=0,ndata-1 do begin

    envi_open_file, data_files[a], r_fid=data_fid
    envi_file_query, data_fid, dims=dims, ns=ns, nl=nl
    data_map_info = envi_get_map_info(fid=data_fid)
    all_data = envi_get_data(fid=data_fid, dims=dims, pos=0)
    where_good = where(all_data ne -1, good_count)   

    envi_convert_projection_coordinates, xvals, yvals, site_proj, dxvals, dyvals, data_map_info.proj
    envi_convert_file_coordinates, data_fid, dximg, dyimg, dxvals, dyvals
    mean_freq = replicate(-1.0d, num_recs)

    kradpix = radius/meter_ps

    data_split = strsplit(file_basename(data_files[a]),'.',/extract)

    for i=0,num_recs-1 do begin
      cur_floor_x = floor(dximg[i]-kradpix[0]) > 0
      cur_ceil_x = ceil(dximg[i]+kradpix[0]) < (ns - 1)
      cur_floor_y = floor(dyimg[i]-kradpix[1]) > 0
      cur_ceil_y = ceil(dyimg[i]+kradpix[1]) < (nl - 1)
      cur_data = all_data[cur_floor_x:cur_ceil_x,cur_floor_y:cur_ceil_y]
      not_bg = where(cur_data ge 0, not_count)    
      if not_count ge 3 then $
        mean_freq[i] = mean(cur_data[not_bg])       
    endfor
    mean_good = where(mean_freq ne -1, mean_count)

    dist_out = prob_output_dir + data_split[0] + '.' + file_titles[a] + '.site_distance.csv'
    openw, lun, dist_out, /get_lun
    for i=0,num_recs-1 do printf, lun, strtrim(mean_freq[i],2)
    free_lun, lun    
          
    ;begin Monte Carlo    
    mcmean = replicate(-1.0d,num_recs,ntrials)   
    
    envi_report_init, base=base, title='Monte Carlo', 'Running simulations....', /interrupt
    envi_report_inc, base, ntrials
  
    ;create a virtual image that is smaller than the original by the amount
    ;we need to make sure we can get same-sized data windows everywhere (2x radius)
    smaller_offset_x = ceil(kradpix[0])
    smaller_offset_y = ceil(kradpix[1])
    smaller_data = all_data[smaller_offset_x:ns-smaller_offset_x-1,smaller_offset_y:nl-smaller_offset_y-1]
    smaller_dims = size(smaller_data, /dimensions)
    smaller_where_good = where(smaller_data ne -1, smaller_good_count)
        
    for j=0,ntrials-1 do begin
      rlocs = smaller_where_good[floor(randomu(seed,num_recs,/double)*smaller_good_count,/l64)]      
      ;add the radius buffer back in to get real image coordinates
      rx = (rlocs mod smaller_dims[0]) + smaller_offset_x
      ry = (rlocs / smaller_dims[0]) + smaller_offset_y
      
      for i=0,num_recs-1 do begin 
        cur_floor_x = floor(rx[i]-kradpix[0]) > 0
        cur_ceil_x = ceil(rx[i]+kradpix[0]) < (ns - 1)
        cur_floor_y = floor(ry[i]-kradpix[1]) > 0
        cur_ceil_y = ceil(ry[i]+kradpix[1]) < (nl - 1)
        cur_data = all_data[cur_floor_x:cur_ceil_x,cur_floor_y:cur_ceil_y]  
        not_bg = where(cur_data ge 0, not_count)
        if not_count ge 3 then $
          mcmean[i,j] = mean(cur_data[not_bg])      
      end
      
      envi_report_stat, base, j+1, ntrials, cancel=cancel
      if (cancel) then begin
        envi_report_init, base=base, /finish
        return        
      endif
    endfor
  
    envi_report_init, base=base, /finish
    
    envi_file_mng, id=data_fid, /remove

    print, data_titles[a]

    ;run K-S test to see what we get
    ;(basic test for difference from random)
    ks_results = dblarr(ntrials)
    ks_testing = dblarr(ncutoffs)
    for j=0,ntrials-1 do begin
      cur_mcmean = mcmean[*,j]
      mcmean_good = where(cur_mcmean ne -1, mcmean_count)
      ks_results[j] = ks_two_sample(mean_freq[mean_good], cur_mcmean[mcmean_good])
      for k=0,ncutoffs-1 do ks_testing[k] += (ks_results[j] lt (1.0d - cutoffs[k]))
    endfor
    ks_testing /= double(ntrials)
    print, 'K-S testing: fraction of trials that passed'
    print, ks_testing
    print, 'Min prob: ' + strtrim(min(ks_results),2), ' Max prob: ' + strtrim(max(ks_results),2)

    ks_out = prob_output_dir + data_split[0] + '.' + file_titles[a] + '.ks_trials.csv'
    openw, lun, ks_out, /get_lun
    for j=0,ntrials-1 do printf, lun, strtrim(ks_results[j],2)
    free_lun, lun    

    ks_test_out = prob_output_dir + data_split[0] + '.' + file_titles[a] + '.ks_test.csv'
    openw, lun, ks_test_out, /get_lun
    printf, lun, strjoin(strtrim(cutoffs, 2),',')
    printf, lun, strjoin(strtrim(ks_testing,2),',')
    free_lun, lun    

    ripley_site_results = dblarr(nsteps)
    ripley_trial_results = dblarr(nsteps,ntrials)
    summary_ripley_test_results = intarr(nsteps,ncutoffs)
    empirical_ripley_test_results = intarr(nsteps,ncutoffs)
    
    ripley_step_mean = dblarr(nsteps)
    ripley_step_std = dblarr(nsteps)
    
    for i=0,nsteps-1 do begin
      ripley_site_results[i] = ripley_l(mean_freq[mean_good], search_radius[i], envelope_area)
      for j=0,ntrials-1 do begin
        cur_mcmean = mcmean[*,j]
        mcmean_good = where(cur_mcmean ne -1, mcmean_count)
        ripley_trial_results[i,j] = ripley_l(cur_mcmean[mcmean_good], search_radius[i], envelope_area)
      endfor
      ripley_step_mean[i] = mean(ripley_trial_results[i,*])
      ripley_step_std[i] = stddev(ripley_trial_results[i,*])
      ;test for statistical significance at several standard cutoffs
      for k=0,ncutoffs-1 do begin
        ;empirical method
        nremove = ceil((ntrials * cutoffs[k]) / 2.0d)
        cur_results = reform(ripley_trial_results[i,*])
        ripley_sort = cur_results[sort(cur_results)]
        ;remove lowest and highest values to create confidence interval
        test_subset = ripley_sort[nremove:ntrials-nremove-1]
        empirical_test_min = min(test_subset, max=empirical_test_max)

        ;summary stats method
        summary_test_prob = gauss_pdf((ripley_site_results[i] - ripley_step_mean[i])/ripley_step_std[i])
        
        ;just testing for clustering, not clustering or dispersion (simple difference from random) 
        empirical_ripley_test_results[i,k] = ripley_site_results[i] gt empirical_test_max
        summary_ripley_test_results[i,k] = summary_test_prob ge (1.0d - cutoffs[k]) 
        
      endfor
    endfor
    
    ;=========================
    ;pool results across steps
    ;=========================
    pooled_ripley_site_mean = mean(ripley_site_results)
    ;empirical method   
    pooled_ripley_empirical_test_results = intarr(ncutoffs)
    for k=0,ncutoffs-1 do begin
      nvals = long64(ntrials)*nsteps
      nremove = ceil((nvals * cutoffs[k]) / 2.0d)
      ripley_reform = reform(ripley_trial_results, nvals)
      ripley_sort = ripley_reform[sort(ripley_reform)]
      ;remove lowest and highest values to create confidence interval
      test_subset = ripley_sort[nremove:nvals-nremove-1]
      empirical_test_min = min(test_subset, max=pooled_empirical_test_max)
      pooled_ripley_empirical_test_results[k] = pooled_ripley_site_mean gt pooled_empirical_test_max
    endfor        
    ;summary stats method
    ;(we can use a simple approach because the number of trials are the same for each step)
    pooled_ripley_step_mean = mean(ripley_step_mean)
    pooled_ripley_step_std = mean(ripley_step_std)
    pooled_ripley_summary_test_results = intarr(ncutoffs)
    pooled_ripley_summary_test_prob = gauss_pdf((pooled_ripley_site_mean - pooled_ripley_step_mean)/pooled_ripley_step_std)
    print, 'Pooled Probability: ' + strtrim(pooled_ripley_summary_test_prob,2)
    print, ' '
    for k=0,ncutoffs-1 do begin
      pooled_ripley_summary_test_results[k] = pooled_ripley_summary_test_prob ge (1.0d - cutoffs[k]) 
    endfor

    print, 'Empirical'
    empirical_results_array = strarr(nsteps+2,ncutoffs+1)
    empirical_results_array[0,0] = '   '
    empirical_results_array[1:nsteps,0] = strtrim(fix(search_radius/1000),2)
    empirical_results_array[nsteps+1,0] = 'Pooled'
    for i=0,ncutoffs-1 do begin
      empirical_results_array[0,i+1] = strtrim(string(cutoffs[i],format='(D5.3)'),2)
      empirical_results_array[1:nsteps,i+1] = strtrim(empirical_ripley_test_results[*,i],2)
      empirical_results_array[nsteps+1,i+1] = strtrim(pooled_ripley_empirical_test_results[i],2)
    endfor
    print, empirical_results_array
    print, ' '
    prob_out = prob_output_dir + data_split[0] + '.' + file_titles[a] + '.empirical.csv'
    openw, lun, prob_out, /get_lun
    for i=0,ncutoffs do printf, lun, strjoin(empirical_results_array[*,i],',')
    free_lun, lun

    print, 'Summary'
    summary_results_array = strarr(nsteps+3,ncutoffs+1)
    summary_results_array[0,0] = '   '
    summary_results_array[1:nsteps,0] = strtrim(fix(search_radius/1000),2)
    summary_results_array[nsteps+1,0] = 'Pooled'
    summary_results_array[nsteps+2,0] = 'Probability'
    for i=0,ncutoffs-1 do begin
      summary_results_array[0,i+1] = strtrim(string(cutoffs[i],format='(D5.3)'),2)
      summary_results_array[1:nsteps,i+1] = strtrim(summary_ripley_test_results[*,i],2)
      summary_results_array[nsteps+1,i+1] = strtrim(pooled_ripley_summary_test_results[i],2)
      ;yes, this gets repeated for every line
      summary_results_array[nsteps+2,i+1] = strtrim(pooled_ripley_summary_test_prob,2)
    endfor
    print, summary_results_array
    print, ' '
    prob_out = prob_output_dir + data_split[0] + '.' + file_titles[a] + '.summary.csv'
    openw, lun, prob_out, /get_lun
    for i=0,ncutoffs do printf, lun, strjoin(summary_results_array[*,i],',')
    free_lun, lun
       
  endfor
  
  envi_file_mng, id=dem_fid, /remove
  
end