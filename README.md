# SolarCalculator Library for Arduino

SolarCalculator is based on the [NOAA Solar Calculator](https://gml.noaa.gov/grad/solcalc/). 

This library provides functions to calculate the times of sunrise, sunset, solar noon and twilight (civil, nautical, 
astronomical dawn and dusk), solar coordinates, interpolation of coordinates, atmospheric refraction correction, Sun's 
radius vector, Julian day and century, equation of time, delta T, etc.

Most formulae are taken from the textbook Astronomical Algorithms by Jean Meeus or the NOAA Solar Calculator 
[source code](https://gml.noaa.gov/grad/solcalc/main.js). Other sources are cited in the comments.


## Installation

Download and copy SolarCalculator to your local Arduino/libraries directory.

### Time and Timezone

Date and time inputs are assumed to be in **Universal Coordinated Time** (UTC).

Although not required, it is recommended to use SolarCalculator along with the
[Time](https://github.com/PaulStoffregen/Time) and [Timezone](https://github.com/JChristensen/Timezone) libraries.


## Usage

Include SolarCalculator.h in your sketch:
```
#include <SolarCalculator.h>
```

Calculate the times of sunrise, transit (solar noon) and sunset:
```
int year = 2022;
int month = 2;
int day = 11;
double latitude = 45.55;
double longitude = -73.633;
int time_zone = -5;

double sunrise, transit, sunset;

calcSunriseSunset(year, month, day, latitude, longitude, transit, sunrise, sunset);
```
where the results are passed by reference, in hours.

Convert to local standard time:
```
double sunrise_local = sunrise + time_zone;
```
To hours and minutes:
```
int minutes = int(round(sunrise_local * 60));
int sunrise_local_hours = (minutes / 60) % 24;
int sunrise_local_minutes = minutes % 60;
```
* Due to the varying effect of refraction, sunrise and sunset times should be rounded to the nearest minute.

Similarly, calculate the times of civil, nautical and astronomical dawn and dusk, in hours:
```
calcCivilDawnDusk(year, month, day, latitude, longitude, transit, c_dawn, c_dusk);
calcNauticalDawnDusk(year, month, day, latitude, longitude, transit, n_dawn, n_dusk);
calcAstronomicalDawnDusk(year, month, day, latitude, longitude, transit, a_dawn, a_dusk);
```

Calculate the Sun's equatorial coordinates (right ascension and declination), in degrees:
```
calcEquatorialCoordinates(year, month, day, hour, minute, second, rt_ascension, declination);
```

Calculate the Sun's horizontal coordinates (azimuth and elevation), in degrees:
```
calcHorizontalCoordinates(year, month, day, hour, minute, second, latitude, longitude, azimuth, elevation);
```

Calculate the Sun's radius vector (distance from the Earth), in AUs:
```
calcSunRadiusVector(year, month, day, hour, minute, second, radius_vector);
```

Calculate the equation of time, in minutes of times:
```
calcEquationOfTime(year, month, day, hour, minute, second, eq);
```
where the results are passed by reference.


## Examples

The following example sketches are included in this library:

* `SunriseSunset`: Calculate the times of sunrise, sunset and solar noon for a given date and location.

* `SolarCalculatorTimeLib`: Calculate the times of sunrise, sunset, solar noon, and the solar coordinates for a given 
location. The Arduino Time library is used for timekeeping purposes.

* `SolarTrackingTimeLib`: Monitor the Sun's position in the sky for any location on Earth.

* `EquationOfTime`: Plot the equation of time for a given year.


## Notes

### Sunrise and sunset

The algorithm for finding the times of sunrise and sunset implemented in this library is valid for all latitudes between 
the Arctic and Antarctic circles (about ± 66.5°). Outside this range, a more general algorithm should be used but is not
provided at this time.

### Accuracy

Various things to consider:

* The amount of atmospheric refraction changes with air temperature, pressure, and the elevation of the observer. 
Therefore, sunrise and sunset times can only be accurate to the nearest minute (Meeus, 1998).

* Assuming a purely elliptical motion of the Earth, solar coordinates have a "low accuracy" of 0.01° (Meeus, 1998).

* Arduino's single precision floating numbers have the equivalent of `24 * log10(2)` ≈ 7.22 significant digits. 
Although this is generally not sufficient for mathematical astronomy (Meeus, 1998), it is sufficient for our solar 
calculations.


## References

Meeus, J. (1998). *Astronomical algorithms* (2nd ed.). Willmann-Bell. <br />
ESRL Global Monitoring Laboratory (n.d.). *NOAA Solar Calculator*. https://gml.noaa.gov/grad/solcalc/
