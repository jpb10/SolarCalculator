# SolarCalculator Library for Arduino

SolarCalculator is inspired by the [NOAA Solar Calculator](https://gml.noaa.gov/grad/solcalc/). 

This library provides functions to calculate the times of sunrise, sunset, solar noon, twilight (dawn and dusk), Sun's 
apparent position in the sky, equation of time, etc.

Most formulae are taken from Astronomical Algorithms by Jean Meeus and optimized for 8-bit AVR platform.


## Installation

Download and copy SolarCalculator to your local Arduino/libraries directory.

### Time

Date and time inputs are assumed to be in **Universal Coordinated Time** (UTC).

Although not required, it is recommended to use SolarCalculator along with the 
[Time](https://github.com/PaulStoffregen/Time) library or similar.


## Usage

Include SolarCalculator.h in your sketch:
```cpp
#include <SolarCalculator.h>
```

Calculate the times of sunrise, transit (solar noon), and sunset, in hours:
```cpp
double latitude = 45.55;     // Observer's latitude 
double longitude = -73.633;  // Observer's longitude
int time_zone = -5;          // UTC offset
int year = 2022;             // Calendar year (1901-2099)
int month = 1;               // Calendar month (1-12)
int day = 1;                 // Calendar day (1-31)

double transit, sunrise, sunset;
calcSunriseSunset(year, month, day, latitude, longitude, transit, sunrise, sunset);
```

Or, using the Time library (Unix time):
```cpp
time_t utc = now();
calcSunriseSunset(utc, latitude, longitude, transit, sunrise, sunset);
```

Convert to local standard time:
```cpp
double sunrise_local = sunrise + time_zone;
```
* Refer to the example sketches for more on rounding and printing the results.

Similarly, calculate the times of dawn and dusk, in hours:
```cpp
calcCivilDawnDusk(utc, latitude, longitude, transit, c_dawn, c_dusk);
calcNauticalDawnDusk(utc, latitude, longitude, transit, n_dawn, n_dusk);
calcAstronomicalDawnDusk(utc, latitude, longitude, transit, a_dawn, a_dusk);
```

Sun's equatorial coordinates, in degrees and AUs:
```cpp
calcEquatorialCoordinates(utc, rt_ascension, declination, radius_vector);
```

Sun's horizontal coordinates, corrected for atmospheric refraction, in degrees:
```cpp
calcHorizontalCoordinates(utc, latitude, longitude, azimuth, elevation);
```

Equation of time, in minutes of time:
```cpp
calcEquationOfTime(utc, eq);
```
where the results are passed by reference.


## Examples

The following example sketches are included in this library:

* `SunriseSunset`: Calculate the times of sunrise, solar noon, and sunset for a given date and location.

* `SunriseSunsetAltitude`: Calculate the rise and set times at a height above the level of the horizon.

* `SolarCalculatorTimeLib`: Calculate the rise and set times, equation of time, and current solar coordinates.

* `SolarTrackingTimeLib`: Monitor the Sun's position in the sky for any location on Earth.

* `EquationOfTime`: Plot the equation of time for a given year.


## Notes

### Accuracy

Various things to consider:

* The amount of atmospheric refraction changes with air temperature, pressure, and the elevation of the observer. 
Therefore, sunrise and sunset times can only be accurate to the nearest minute (Meeus, 1998).

* Assuming a purely elliptical motion of the Earth, solar coordinates have a "low accuracy" of 0.01° (Meeus, 1998). To
this precision, we ignore nutation, delta T, and higher-order terms in the relevant expressions.

* Arduino's single precision floating numbers have the equivalent of `23 * log10(2)` ≈ 6.92 significant digits. 
Although this is generally not sufficient for mathematical astronomy (Meeus, 1998), it is good enough for our purposes.

### Sunrise and sunset

The algorithm for finding the times of sunrise and sunset implemented in this library is valid for all latitudes between 
the Arctic and Antarctic circles (about ± 66.5°). Outside this range, a more general algorithm should be used but is not
provided at this time.


## References

ESRL Global Monitoring Laboratory (n.d.). *NOAA Solar Calculator*. https://gml.noaa.gov/grad/solcalc/

Meeus, J. (1998). *Astronomical algorithms* (2nd ed.). Willmann-Bell.

