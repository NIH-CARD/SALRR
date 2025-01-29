#!/usr/bin/env bash

#script for finding the directory paths and the json files and creating a list of json files
#The output  "json_report_paths.txt" will be used as the input  file for the longread-report-parser

find /data/CARD_AUX/LRS_temp/NABEC_RNA/ -name "*.json" > /data/CARD_AUX/LRS_temp/NABEC_RNA/SEQ_REPORTS/json_report_paths.txt


