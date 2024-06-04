# IBEX_Calibration

Code performing computer model calibration on IBEX satellite data and simulator output

## Data Details

These datasets have been approved for public release under LA-UR-24-25151.

### ibex_real.csv

Contains record of data collected on ENAs from the IBEX satellite in the years 2009-2022.

- esa: energy step for the IBEX data. Only ESA 4 is provided in this dataset
- map: the 6-month map id (e.g., "2013A"). Only "A" maps are provided
- ecliptic_lon: the ecliptic longitude (values range between 0 and 360)
- ecliptic_lat: the ecliptic latitude (values range between -90 and 90)
- ecliptic_lon_center: the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- x: are the sperical coordinates for the ecliptic longitude and latitude
- y
- z
- counts: the number of ENAs detected
- time: the amount of time the satellite was pointed at the given lat lon (in seconds)
- background: the estimated background rate (background ENAs per second)

### sims.csv

Contains output from two computer models (GDF, ribbon) added together to form different simulated skymaps based on varying two parameters.

- lat: the ecliptic latitude (values range between -90 and 90)
- lon: the ecliptic longitude (values range between 0 and 360)
- ESA: energy step for the IBEX data. Only ESA 4 is provided in this dataset
- ena_rate: rate at which ENAs are generated at a given location (lat, lon)
- parallel_mean_free_path: parameter representing the distance a charge particle travels before encountering a particle exchange
- ratio: additional model parameter that is proportional to another value (perpendicular mean free path)
- ecliptic_lon_center: the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- blurred_ena_rate: simulation of an observed ENA rate after passing through the point spread function

### sims_real.csv

Contains simulated counts repesenting the data collected by the IBEX satellite, but based on ENA rates from different simulated skymaps.

- parallel_mean_free_path: parameter representing the distance a charge particle travels before encountering a particle exchange
- ratio: additional model parameter that is proportional to another value (perpendicular mean free path)
- lon: the ecliptic longitude (values range between 0 and 360)
- lat: the ecliptic latitude (values range between -90 and 90)
- esa: energy step for the IBEX data. Only ESA 4 is provided in this dataset
- time: the amount of time the satellite was pointed at the given lat lon (in seconds)
- blurred_ena_rate: simulation of an observed ENA rate after passing through the point spread function
- background: the estimated background rate (background ENAs per second)
- sim_counts: random draw of counts based on a simulated blurred ENA rate

---

Copyright 2023 for **O4627**

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
