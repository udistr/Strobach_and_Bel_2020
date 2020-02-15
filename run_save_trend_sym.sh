#!/bin/bash

st=0
en=33

run_save_EGA=0
run_save_AVG=0
run_plot=1
run_rel=0


if [ ${run_save_EGA} -eq 1 ] ; then

  ncl save_ts_trend_runmean_20.ncl rcp=0 'SLA="EGA"'
  ncl save_ts_trend_runmean_20.ncl rcp=1 'SLA="EGA"'
  ncl save_ts_trend_runmean_20.ncl rcp=2 'SLA="EGA"'
  ncl save_ts_trend_runmean_20.ncl rcp=3 'SLA="EGA"'
  ncl save_ts_trend_historica_runmean_20.ncl rcp=1 'SLA="EGA"'
  ncl save_ts_trend_historica_runmean_20.ncl rcp=1 'SLA="AVG"'
fi


if [ ${run_save_AVG} -eq 1 ] ; then

  ncl save_ts_trend_runmean_gaus_20.ncl rcp=0 'SLA="AVG"'
  ncl save_ts_trend_runmean_gaus_20.ncl rcp=1 'SLA="AVG"'
  ncl save_ts_trend_runmean_gaus_20.ncl rcp=2 'SLA="AVG"'
  ncl save_ts_trend_runmean_gaus_20.ncl rcp=3 'SLA="AVG"'

fi

if [ ${run_plot} -eq 1 ] ; then

  ncl plot_ts_trend_lim_rcp_runmean_20_2.ncl 'SLA="EGA"'
  ncl plot_ts_trend_lim_rcp_runmean_gaus_20.ncl 'SLA="AVG"'
  matlab -r "run plot_cal_range_runmean_20b.m; exit;"
  matlab -r "plot_dist_rec_combine_runmean_20_2.m; exit"
  matlab -r "plot_dist_rec_combine_runmean_20.m; exit;"
fi


if [ ${run_rel} -eq 1 ] ; then

  for i in `seq ${st} ${en}` ; do

    ncl plot_relibility.ncl vmodel=[$i]

  done

fi

