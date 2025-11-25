## Data Details

The datasets in this directory have been appoved for public release under
LA-UR-24-25151. Below are descriptions of each file:

### ibex_real.csv

Contains record of data collected on energetic neutral atoms (ENA) from the
Interstellar Boundary Explorer (IBEX) satellite for the years 2009-2022.

- esa: energy step for ENAs (one of six levels) detected by IBEX. Only ESA 4 is
  provided in this dataset
- map: the 6-month map id (e.g., "2013A"). Only "A" maps are provided
- ecliptic_lon: the ecliptic longitude (values range between 0 and 360) at
  which the satellite was pointed
- ecliptic_lat: the ecliptic latitude (values range between -90 and 90) at
  which the satellite was pointed
- ecliptic_lon_center: the ecliptic longitude in a heliosphere "nose centered"
  frame that is used for plotting purposes
- x, y, z: the corresponding sperical coordinates for the ecliptic longitude
  and latitude
- counts: the number of ENAs detected by the satellite at this row's location
- time: the amount of time the satellite was pointed at the given lat lon (in
  seconds)
- background: the estimated background rate (background ENAs per second)

---

### sims.csv

Contains output from two computer models (GDF, ribbon) added together to form
different simulated sky maps based on varying two parameters.

- lat: the ecliptic latitude (values range between -90 and 90)
- lon: the ecliptic longitude (values range between 0 and 360)
- ESA: energy step for the ENAs. Only ESA 4 is provided in this dataset
- ena_rate: rate at which ENAs are generated at the given location (lat, lon)
- parallel_mean_free_path: parameter representing the distance a charge
  particle travels before undergoing an electron exchange and 'turning back'
- ratio: additional model parameter that is the fraction of perpendicular and
  parallel mean free path
- ecliptic_lon_center: the ecliptic longitude in a heliosphere "nose centered"
  frame and is used for plotting purposes
- blurred_ena_rate: simulated ENA rate after passing through the point spread
  function

---

### synth_sat_data.csv

Contains synthetic counts attempting to mimic the data collected by the IBEX
satellite, but based on ENA rates from different simulated sky maps. Quantities
such as latitude, longitude, time, and background are pulled from satellite
data for a specific year/

- parallel_mean_free_path: parameter representing the distance a charge
  particle travels before undergoing an electron exchange and 'turning back'
- ratio: additional model parameter that is the fraction of perpendicular and
  parallel mean free path
- lon: the ecliptic longitude (values range between 0 and 360)
- lat: the ecliptic latitude (values range between -90 and 90)
- esa: energy step for the ENAs. Only ESA 4 is provided in this dataset
- time: the amount of time the satellite was pointed at the given lat lon (in
  seconds)
- blurred_ena_rate: simulated ENA rate after passing through the point spread
  function
- background: the estimated background rate (background ENAs per second)
- sim_counts: random Poisson draw of counts based on a simulated blurred ENA
  rate and exposure time
