#!/bin/bash

wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/20237715/primaryscreenreplicatecollapsedtreatmentinfo.csv -O source/prism/primary-screen-replicate-collapsed-treatment-info.csv;
wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/20237709/primaryscreenreplicatecollapsedlogfoldchange.csv -O source/prism/primary-screen-replicate-collapsed-logfold-change.csv;
wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/20237718/primaryscreencelllineinfo.csv -O source/prism/primary-screen-cell-line-info.csv;
wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/20237739/secondaryscreendoseresponsecurveparameters.csv  -O source/prism/secondary-screen-dose-response-curve-parameters.csv;
wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/20237757/secondaryscreenreplicatecollapsedlogfoldchange.csv -O source/prism/secondary-screen-replicate-collapsed-logfold-change.csv;
wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/20237763/secondaryscreenreplicatecollapsedtreatmentinfo.csv -O source/prism/secondary-screen-replicate-collapsed-treatment-info.csv
