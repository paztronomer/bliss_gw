# BLISS Gravitational Waves
## Codes to get the observability for set of RA-DEC-EXPNUM

1. **create_json.py**: based on inputs RA-DEC and exposure number (EXPNUM),
    calculates the observability at the lowest airmass, for a set of input
    nights. The code writes JSON files (as used on CTIO-4m). Many JSON fields
    are able to be modified through the input arguments. Use the `--help` for
    further information. Horizon limits are not a simple circle, but the ones
    CTIO has listed in its webpage.

    A simple call could be:
    ```bash
    python create_json.py event3modif_ccds.txt observing_2017A.txt
    --root_csv info_gw3 --root_json obs_gw3 --sequence 'GW event 3'
    --proposal 2017A-0260 --program BLISS —req_one_band
    ```
1. **footprint_horizon.py**: print the horizon limits, based in the limits
    given by CTIO. It differs from other codes versions, on which a simple
    circle is used. *NEED MORE TEST AND FIX IF ANY PROBLEM*
