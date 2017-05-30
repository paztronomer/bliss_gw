# bliss_gw
Codes for Gravitational Waves in BLISS

(Draft!)

usage: create_json.py [-h] [--source] [--out] [--root_csv] [--root_json]
                      [--utc_diff] [--t_step] [--max_airmass] [--count]
                      [--sequence] [--proposal] [--exptype] [--program]
                      [--band] [--exptime] [--tiling] [--comment] [--note]
                      [--wait]
                      objects [objects ...]

Script to calculate observability from CTIO, using Blanco 4m telecope setup,
for a set (or single) of coordinates given the observing date. There are 2
type of optional arguments: those wo serve as setup and those who will be
directly inserted in the JSON files. NOTE: for JSON arguments use quotes if
any space is in the string

positional arguments:

  objects              Set of space-separated filenames for the tables
                       containig the coordinates to be calculated. Format:
                       tables must have RA and DEC in header, despite the
                       other column names; separator is expected to be spaces
                       Unique file containig all the nights for which the
                       observability will be calculated. Format: 2 columns
                       being YYYY-MM-DD {first/second/all}

optional arguments:

  -h, --help           show this help message and exit
  
  --source , -s        Path to folder containing the source tables for objects
                       and dates. Default is current directory
                       
  --out , -o           Path to folder for output files (JSON and CSV). Default
                       is current directory
                       
  --root_csv           Root string of output name(s) for the file(s)
                       containing the coordinates that passed the
                       observability criteria, plus additional information.
                       Default is 'sel_info'
                       
  --root_json          Root string of JSON output files (if objects were
                       found). Default is 'obs'
                       
  --utc_diff , -u      Difference in hours between UTC and the observing
                       location, namely, CTIO (value=UTC-LOCAL). Default: 4
                       
  --t_step , -t        Step (in minutes) at which the night will be sampled to
                       calculate the observability at each interval. Default:
                       10
                       
  --max_airmass , -m   Maximum airmass at which the objects want to be
                       observed. Default: 1.8
                       
  --count              JSON optional.Number of exposures to be taken for each
                       object. Default: 1
                       
  --sequence           JSON required. Sequence ID, eg: 'LIGO event x'
  
  --proposal           JSON required. Proposal ID
  
  --exptype            JSON optional. Exposure type. Default: 'object'
  
  --program            JSON required. Program ID, example: 'BLISS'
  
  --band               JSON optional. Band to be used. Default: i
  
  --exptime            JSON optional. Exposure time in seconds. Default: 90
  
  --tiling             JSON optional. ID of the tiling. Default: 1
  
  --comment            JSON Optional
  
  --note               JSON optional. Note to be added. Default: 'Added to
                       queue by user, not obstac'
                       
  --wait               JSON optional. Whether to wait to proceed for next
                       exposure. Default: 'False'
