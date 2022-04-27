//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SunriseSunsetAltitude.ino
//
// Calculate the rise and set times at a height above the level of the horizon.
//
// Tested with Arduino IDE 1.8.19 and Arduino Uno
//======================================================================================================================

#include <SolarCalculator.h>

void setup()
{
  Serial.begin(9600);

  // Date
  int year = 2022;
  int month = 1;
  int day = 1;

  // Location
  double latitude = 45.5034;
  double longitude = -73.5869;
  int utc_offset = -5;

  double transit, sunrise, sunset;

  // From the Explanatory Supplement to the Astronomical Almanac (1992), p. 484
  // Sunrise or sunset at a height above the level of the horizon occurs when the Sun's altitude is approximately:

  int height = 200;  // in meters
  double sun_altitude = SUNRISESET_STD_ALTITUDE - 0.0353 * sqrt(height);

  // Calculate the times of sunrise, transit, and sunset, in hours (UTC)
  calcSunriseSunset(year, month, day, latitude, longitude, transit, sunrise, sunset, sun_altitude);

  // Print results
  char str[6];
  Serial.println(hoursToString(sunrise + utc_offset, str));
  Serial.println(hoursToString(transit + utc_offset, str));
  Serial.println(hoursToString(sunset + utc_offset, str));
}

void loop()
{
}

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
