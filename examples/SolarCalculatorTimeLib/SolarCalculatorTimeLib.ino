//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SolarCalculatorTimeLib.ino
//
// Calculate the times of sunrise, sunset, solar noon, and the solar coordinates for a given location.
//
// The Arduino Time library is used for timekeeping. For more accurate timekeeping, system time should be synchronized
// with an RTC module and/or NTP server.
//
// Tested with Arduino IDE 1.8.19 and Arduino Uno
//======================================================================================================================

#include <SolarCalculator.h>
#include <TimeLib.h>

// Location
const double latitude = 45.55;
const double longitude = -73.633;
const int time_zone = -5;

void setup() 
{
  Serial.begin(9600);

  time_t sunrise;        // Sunrise (UTC)      
  time_t transit;        // Solar noon (UTC)
  time_t sunset;         // Sunset (UTC)
  double eq;             // Equation of time, in minutes
  double rt_ascension;   // Sun's right ascension, in degrees
  double declination;    // Sun's declination, in degrees
  double azimuth;        // Sun's azimuth, in degrees
  double elevation;      // Sun's elevation, in degrees

  // Set system time to compile time (UTC)
  setTime(toUtc(compileTime()));
  
  // Or, set time manually (hr, min, sec, day, mo, yr)
  //setTime(12, 0, 0, 11, 2, 2022);

  // Get current time
  time_t utc = now();

  double transit_h, sunrise_h, sunset_h;

  // Calculate the times of sunrise, transit and sunset, in hours (UTC)
  calcSunriseSunset(year(utc), month(utc), day(utc), latitude, longitude, transit_h, sunrise_h, sunset_h);
  
  // Convert hours to time_t
  sunrise = hoursToTime(utc, sunrise_h);
  transit = hoursToTime(utc, transit_h);
  sunset = hoursToTime(utc, sunset_h);

  // Calculate the equation of time
  calcEquationOfTime(year(utc), month(utc), day(utc), hour(utc), minute(utc), second(utc), eq);

  // Calculate the Sun's right ascension and declination
  calcEquatorialCoordinates(year(utc), month(utc), day(utc), hour(utc), minute(utc), second(utc), 
                            rt_ascension, declination);

  // Calculate the Sun's azimuth and elevation
  calcHorizontalCoordinates(year(utc), month(utc), day(utc), hour(utc), minute(utc), second(utc), latitude, longitude, 
                            azimuth, elevation);
  
  // Print results
  Serial.print("Date: ");
  Serial.print(year(utc));
  Serial.print('-');
  printDigits(month(utc));
  Serial.print('-');
  printDigits(day(utc));
  Serial.println();
  
  Serial.print("Latitude: ");
  Serial.print(latitude, 3);
  Serial.print(" Longitude: ");
  Serial.println(longitude, 3);
  Serial.print("UTC offset: ");
  Serial.println(time_zone);
  Serial.println("--"); 

  Serial.print("Sunrise: ");
  printSunTime24h(toLocal(sunrise));
  Serial.print("Transit: ");
  printSunTime24h(toLocal(transit));
  Serial.print("Sunset:  ");
  printSunTime24h(toLocal(sunset));
  Serial.println("--"); 

  Serial.print("Equation of time: ");
  Serial.print(eq);
  Serial.println(" minutes");
  Serial.print("Right ascension: ");
  printCoordHours(rt_ascension);
  Serial.print(" Declination: ");
  printCoordDegrees(declination);
  Serial.println();
  Serial.print("Azimuth: ");
  Serial.print(azimuth);
  Serial.print("° Elevation: ");
  Serial.print(elevation);
  Serial.println("°");
}

void loop() 
{
}

time_t toUtc(time_t local)
{
  return local - time_zone * 3600;
}

time_t toLocal(time_t utc)
{
  return utc + time_zone * 3600;
}

time_t hoursToTime(time_t day, double t)
{
  return previousMidnight(day) + time_t(round(t * 3600));
}

// Code from JChristensen/Timezone Clock.ino example
time_t compileTime()
{
  const uint8_t COMPILE_TIME_DELAY = 5;
  const char *compDate = __DATE__, *compTime = __TIME__, *months = "JanFebMarAprMayJunJulAugSepOctNovDec";
  char chMon[4], *m;
  tmElements_t tm;

  strncpy(chMon, compDate, 3);
  chMon[3] = '\0';
  m = strstr(months, chMon);
  tm.Month = ((m - months) / 3 + 1);

  tm.Day = atoi(compDate + 4);
  tm.Year = atoi(compDate + 7) - 1970;
  tm.Hour = atoi(compTime);
  tm.Minute = atoi(compTime + 3);
  tm.Second = atoi(compTime + 6);
  time_t t = makeTime(tm);
  return t + COMPILE_TIME_DELAY;
}

double degreesToHours(double deg)
{
  return 24 * wrapTo360(deg) / 360;
}

void printCoordHours(double d)
{
  long s = long(round((degreesToHours(d) * 3600)));
  int seconds = int(s % 60);
  int minutes = int(s / 60);
  int hours = minutes / 60; 
  minutes = minutes % 60;
  Serial.print(hours);
  Serial.print("h ");
  Serial.print(minutes);
  Serial.print("m ");
  Serial.print(seconds);
  Serial.print("s");
}

void printCoordDegrees(double d)
{
  long s = long(round(d * 3600));
  int arcsec = int(abs(s) % 60);
  int arcmin = int(s / 60);
  int degrees = arcmin / 60; 
  arcmin = abs(arcmin) % 60;
  Serial.print(degrees);
  Serial.print("° ");
  Serial.print(arcmin);
  Serial.print("' ");
  Serial.print(arcsec);
  Serial.print("\"");
}

void printSunTime24h(time_t t) 
{
  if (second(t) >= 30)
    t += 60;
  if (hour(t) < 10)
    Serial.print(' ');
  Serial.print(hour(t));
  Serial.print(':');
  printDigits(minute(t));
  Serial.println();
}

void printDigits(int digits)
{
  if (digits < 10)
    Serial.print('0');
  Serial.print(digits);
}
