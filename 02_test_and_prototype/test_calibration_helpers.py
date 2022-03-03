# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %%
import prms_lib as prms
import prms_calib_helpers as hlp


# %%
def_ranges = {'ssstor_init_frac': {'max': 1.0, 'min': 0.0}, 'tstorm_mo': {'max': 1.0, 'min': 0.0}, 'parent_ssr': {'max': 1000000.0, 'min': 1.0}, 'va_open_exp': {'max': 10.0, 'min': 0.001}, 'nhm_seg': {'max': 9999999.0, 'min': 1.0}, 'jh_coef_hru': {'max': 150.0, 'min': -99.0}, 'pref_flow_den': {'max': 0.5, 'min': 0.0}, 'soil_moist_max': {'max': 10.0, 'min': 1e-05}, 'snarea_curve': {'max': 1.0, 'min': 0.0}, 'covden_sum': {'max': 1.0, 'min': 0.0}, 'hru_percent_imperv': {'max': 0.999, 'min': 0.0}, 'soil2gw_max': {'max': 5.0, 'min': 0.0}, 'dprst_frac': {'max': 0.999, 'min': 0.0}, 'temp_units': {'max': 1.0, 'min': 0.0}, 'epan_coef': {'max': 3.0, 'min': 0.01}, 'elev_units': {'max': 1.0, 'min': 0.0}, 'tmax_allrain_offset': {'max': 20.0, 'min': 0.0}, 'freeh2o_cap': {'max': 0.2, 'min': 0.01}, 'dday_slope': {'max': 0.9, 'min': 0.2}, 'print_type': {'max': 2.0, 'min': 0.0}, 'adjmix_rain': {'max': 1.4, 'min': 0.6}, 'rain_cbh_adj': {'max': 2.0, 'min': 0.5}, 'segment_type': {'max': 1000.0, 'min': 0.0}, 'smidx_coef': {'max': 0.06, 'min': 0.0001}, 'print_freq': {'max': 15.0, 'min': 0.0}, 'settle_const': {'max': 0.5, 'min': 0.01}, 'rad_trncf': {'max': 1.0, 'min': 0.0}, 'precip_units': {'max': 1.0, 'min': 0.0}, 'gwsink_coef': {'max': 0.05, 'min': 0.0}, 'cov_type': {'max': 4.0, 'min': 0.0}, 'parent_hru': {'max': 1000000.0, 'min': 1.0}, 'hru_lat': {'max': 90.0, 'min': -90.0}, 'K_coef': {'max': 24.0, 'min': 0.01}, 'tmax_index': {'max': 110.0, 'min': -10.0}, 'den_init': {'max': 0.5, 'min': 0.01}, 'obsin_segment': {'max': 1.0, 'min': 0.0}, 'op_flow_thres': {'max': 1.0, 'min': 0.75}, 'transp_tmax': {'max': 1000.0, 'min': 0.0}, 'hru_type': {'max': 3.0, 'min': 0.0}, 'dprst_seep_rate_clos': {'max': 0.1, 'min': 0.0001}, 'hru_x': {'max': 10000000.0, 'min': -10000000.0}, 'ssr2gw_exp': {'max': 3.0, 'min': 0.0}, 'smidx_exp': {'max': 0.5, 'min': 0.0001}, 'dprst_frac_open': {'max': 1.0, 'min': 0.0}, 'poi_type': {'max': 1.0, 'min': 1.0}, 'radj_wppt': {'max': 1.0, 'min': 0.0}, 'snow_intcp': {'max': 1.0, 'min': 0.0}, 'hru_lon': {'max': 360.0, 'min': -360.0}, 'radadj_intcp': {'max': 1.0, 'min': 0.0}, 'snowinfil_max': {'max': 20.0, 'min': 0.0}, 'sat_threshold': {'max': 999.0, 'min': 1.0}, 'tmax_allsnow': {'max': 40.0, 'min': -10.0}, 'cecn_coef': {'max': 10.0, 'min': 2.0}, 'dprst_et_coef': {'max': 1.0, 'min': 0.0}, 'hru_segment': {'max': 1.0, 'min': 0.0}, 'ppt_rad_adj': {'max': 0.5, 'min': 0.0}, 'dday_intcp': {'max': 10.0, 'min': -60.0}, 'radmax': {'max': 1.0, 'min': 0.1}, 'soil_rechr_init_frac': {'max': 1.0, 'min': 0.0}, 'slowcoef_lin': {'max': 0.5, 'min': 1e-05}, 'den_max': {'max': 0.8, 'min': 0.1}, 'x_coef': {'max': 0.5, 'min': 0.0}, 'radj_sppt': {'max': 1.0, 'min': 0.0}, 'tosegment': {'max': 1.0, 'min': 0.0}, 'parent_poigages': {'max': 1000000.0, 'min': 1.0}, 'snowpack_init': {'max': 5000.0, 'min': 0.0}, 'dprst_depth_avg': {'max': 500.0, 'min': 0.0}, 'sro_to_dprst_imperv': {'max': 1.0, 'min': 0.0}, 'gwstor_min': {'max': 1.0, 'min': 0.0}, 'hru_elev': {'max': 30000.0, 'min': -1000.0}, 'hru_deplcrv': {'max': 2.0, 'min': 0.0}, 'soil_moist_init_frac': {'max': 1.0, 'min': 0.0}, 'segment_flow_init': {'max': 10000000.0, 'min': 0.0}, 'outlet_sta': {'max': 1.0, 'min': 0.0}, 'fastcoef_sq': {'max': 1.0, 'min': 1e-05}, 'carea_max': {'max': 1.0, 'min': 0.0}, 'gwflow_coef': {'max': 0.5, 'min': 0.0}, 'gwstor_init': {'max': 10.0, 'min': 0.0}, 'nhm_id': {'max': 9999999.0, 'min': 1.0}, 'hru_aspect': {'max': 360.0, 'min': 0.0}, 'tmax_cbh_adj': {'max': 10.0, 'min': -10.0}, 'parent_gw': {'max': 1000000.0, 'min': 1.0}, 'tosegment_nhm': {'max': 9999999.0, 'min': 1.0}, 'dprst_frac_init': {'max': 1.0, 'min': 0.0}, 'runoff_units': {'max': 1.0, 'min': 0.0}, 'tmin_cbh_adj': {'max': 10.0, 'min': -10.0}, 'poi_gage_id': {'max': 99999999.0, 'min': 0.0}, 'covden_win': {'max': 1.0, 'min': 0.0}, 'va_clos_exp': {'max': 10.0, 'min': 0.001}, 'jh_coef': {'max': 0.1, 'min': 0.0}, 'imperv_stor_max': {'max': 0.1, 'min': 0.0}, 'wrain_intcp': {'max': 1.0, 'min': 0.0}, 'sro_to_dprst_perv': {'max': 1.0, 'min': 0.0}, 'melt_force': {'max': 366.0, 'min': 1.0}, 'dprst_flow_coef': {'max': 0.3, 'min': 0.0001}, 'snarea_thresh': {'max': 200.0, 'min': 0.0}, 'potet_sublim': {'max': 0.75, 'min': 0.1}, 'ssr2gw_rate': {'max': 0.8, 'min': 0.05}, 'soil_rechr_max_frac': {'max': 1.0, 'min': 1e-05}, 'slowcoef_sq': {'max': 1.0, 'min': 1e-05}, 'emis_noppt': {'max': 1.0, 'min': 0.757}, 'srain_intcp': {'max': 1.0, 'min': 0.0}, 'radadj_slope': {'max': 1.0, 'min': 0.0}, 'hru_slope': {'max': 10.0, 'min': 0.0}, 'soil_type': {'max': 3.0, 'min': 1.0}, 'albset_sna': {'max': 1.0, 'min': 0.01}, 'parent_segment': {'max': 1000000.0, 'min': 1.0}, 'albset_rna': {'max': 1.0, 'min': 0.5}, 'melt_look': {'max': 366.0, 'min': 1.0}, 'albset_rnm': {'max': 1.0, 'min': 0.4}, 'albset_snm': {'max': 1.0, 'min': 0.1}, 'fastcoef_lin': {'max': 1.0, 'min': 1e-05}, 'snow_cbh_adj': {'max': 2.0, 'min': 0.5}, 'hru_y': {'max': 10000000.0, 'min': -10000000.0}, 'hru_area': {'max': 1000000000.0, 'min': 0.0001}, 'dprst_seep_rate_open': {'max': 0.1, 'min': 0.0001}, 'transp_end': {'max': 13.0, 'min': 1.0}, 'transp_beg': {'max': 12.0, 'min': 1.0}, 'poi_gage_segment': {'max': 1.0, 'min': 0.0}}

sens_file = '/Users/pnorton/Projects/National_Hydrology_Model/GCPO/test_set/hruSens.csv'

# sens_params = {' tmax_cbh_adj': 1, ' snow_cbh_adj': 1, ' cecn_coef': 1, ' potet_sublim': 1, ' freeh2o_cap': 1, ' tmin_cbh_adj': 1, ' emis_noppt': 1, 'rain_cbh_adj': 1, ' radmax': 1}

pfile = '/Users/pnorton/Projects/National_Hydrology_Model/GCPO/test_set/rGCPO_000000/gcpo.params.expanded'

ofile = '/Users/pnorton/Projects/National_Hydrology_Model/GCPO/test_set/param_limits.txt'

# %%
reload(prms)
reload(hlp)

sens_params = hlp.read_sens_params(sens_file)
hlp.adjust_param_ranges(pfile, sens_params, def_ranges, ofile, make_dups=False)

# %%
