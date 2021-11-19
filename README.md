# SolarCalculator Library for Arduino

SolarCalculator is based on the 
[NOAA Solar Calculator](https://www.esrl.noaa.gov/gmd/grad/solcalc/). 

This library provides functions to calculate the Sun's position in the sky, the times of 
sunrise, sunset, twilight and solar noon for any location on earth, as well as the equation 
of time and more.

Most formulae are taken from the textbook Astronomical Algorithms by Jean Meeus 
or the NOAA Solar Calculator [source code](https://gml.noaa.gov/grad/solcalc/main.js). 
Other sources are cited in the comments.

## Installation

Download and copy SolarCalculator to your local Arduino/libraries directory.

### Time and Timezone

Date and time inputs are assumed to be in **Universal Coordinated Time** (UTC).

Although not required, it is recommended to use SolarCalculator along with the
[Time](https://github.com/PaulStoffregen/Time) and 
[Timezone](https://github.com/JChristensen/Timezone) libraries.

## Usage

Include SolarCalculator.h in your sketch:
```
#include <SolarCalculator.h>
```

Calculate the times of sunrise, transit (solar noon) and sunset:
```
int year = 2021;
int month = 4;
int day = 30;
double latitude = 45.55;
double longitude = -73.633;
int time_zone = -4;

double sunrise;
double transit; 
double sunset; 

calcSunriseSunset(year, month, day, latitude, longitude, transit, sunrise, sunset);
```
where the results are passed by reference, in hours.

Convert to local standard time, then to hours and minutes:
```
double sunrise_local = sunrise + time_zone;

int minutes = int(round(sunrise_local * 60));
int sunrise_local_hours = (minutes / 60) % 24;
int sunrise_local_minutes = minutes % 60;
```
* Due to the varying effect of refraction, sunrise and sunset times should only be rounded
to the nearest minute.


Similarly, calculate the times of civil, nautical and astronomical dawn and dusk:
```
calcCivilDawnDusk(year, month, day, latitude, longitude, transit, c_dawn, c_dusk);
calcNauticalDawnDusk(year, month, day, latitude, longitude, transit, n_dawn, n_dusk);
calcAstronomicalDawnDusk(year, month, day, latitude, longitude, transit, a_dawn, a_dusk);
```

Calculate the Sun's equatorial coordinates (right ascension and declination):
```
calcEquatorialCoordinates(year, month, day, hour, minute, second, rt_ascension, declination);
```

Calculate the Sun's horizontal coordinates (azimuth and elevation):
```
calcHorizontalCoordinates(year, month, day, hour, minute, second, latitude, longitude, azimuth, elevation);
```

Calculate the equation of (ephemeris) time:
```
calcEquationOfTime(year, month, day, hour, minute, second, eq);
```
