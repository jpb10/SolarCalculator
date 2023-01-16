//======================================================================================================================
// SolarCalculator Library example for Linux: Main.cpp
//
// Calculate the rise and set times at a height above the level of the horizon.
//
// Tested with Kubuntu 20.04, compiled using: gcc Main.cpp ../../src/SolarCalculator.cpp -I../../src -o calcSolar -lm
//======================================================================================================================

#include "SolarCalculator.h"
#include <math.h>
#include <time.h>
#include <stdio.h>


// Location: Customize these location parameters for your location on Earth
// you can get it from https://gml.noaa.gov/grad/solcalc/
static double latitude = 32.803200; 
static double longitude = -117.13105;
static int utc_offset = -7; // used for correcting UTC to local time for sunrise and sunset display
//
// Set your location elevation above sealevel for correction to the horizon
// From the Explanatory Supplement to the Astronomical Almanac (1992), p. 484
// Sunrise or sunset at a height above the level of the horizon occurs when the Sun's altitude is approximately:
static int height = 10;  // in meters

// Rounded HH:mm format
char * hoursToString(double h, char *str)
{
  int m = int(round(h * 60));
  int hr = (m / 60) % 24;
  int mn = m % 60;

  str[0] = (hr / 10) % 10 + '0';
  str[1] = (hr % 10) + '0';
  str[2] = ':';
  str[3] = (mn / 10) % 10 + '0';
  str[4] = (mn % 10) + '0';
  str[5] = '\0';
  return str;
}

time_t toLocal(time_t utc)
{
  return utc + utc_offset * 3600;
}


int main()
{
  // Date
  time_t now, utc_now;
  struct tm * timeinfo;
  struct tm * timeinfo_utc;
  time(&now);
  timeinfo = localtime(&now);

  int hours,houru; // storing time for displaying local time and using UTC for calculations
  int seconds,secondu;
  int minutes,minuteu;
  int days,dayu;
  int months,monthu;
  int years,yearu;

  years = timeinfo->tm_year + 1900;   //https://mikaelpatel.github.io/Arduino-RTC/d8/d5a/structtm.html
  months = timeinfo->tm_mon + 1;
  days = timeinfo->tm_mday;
  hours = timeinfo->tm_hour;
  minutes = timeinfo->tm_min;
  seconds = timeinfo->tm_sec;

  utc_now = time( NULL );
  timeinfo_utc = gmtime(&utc_now);

  yearu = timeinfo_utc->tm_year + 1900;   //https://mikaelpatel.github.io/Arduino-RTC/d8/d5a/structtm.html
  monthu = timeinfo_utc->tm_mon + 1;
  dayu = timeinfo_utc->tm_mday;
  houru = timeinfo_utc->tm_hour;
  minuteu = timeinfo_utc->tm_min;
  secondu = timeinfo_utc->tm_sec;

  double transit, sunrise, sunset;
  double sun_altitude = SUNRISESET_STD_ALTITUDE - 0.0353 * sqrt(height);
  double sunElevation;
  double sunAzimuth;

  // Calculate the times of sunrise, transit, and sunset, in hours (UTC)
  calcSunriseSunset(yearu, monthu, dayu, latitude, longitude, transit, sunrise, sunset, sun_altitude);

  // Calculate the Sun's azimuth and elevation (corrected for atmospheric refraction), in degrees
  calcHorizontalCoordinates(yearu, monthu, dayu, houru, minuteu, secondu, latitude, longitude, sunAzimuth, sunElevation);

  printf("Day: %d: %d: %d\n",months,days,years);
  printf("Time: %d: %d: %d\n",hours,minutes,seconds);
  printf("UTC Day: %d: %d: %d\n",monthu,dayu,yearu);
  printf("UTC Time: %d: %d: %d\n",houru,minuteu,secondu);

  // Print results
  char str[6];
  printf("SunRise:");
  printf("%s\n",hoursToString(sunrise + utc_offset, str));
  printf("Transit:");
  printf("%s\n",hoursToString(transit + utc_offset, str));
  printf("SunSet:");
  printf("%s\n",hoursToString(sunset + utc_offset, str));

  printf("Sun angles : \n");
  printf("    Elevation : ");
  printf("%f",sunElevation);
  printf("    panel Elevation : ");
  printf("%f\n",90 - sunElevation);
  printf("    Azimuth   : ");
  printf("%f\n",sunAzimuth);
}

