//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SolarCalculatorTimeLib
//  
// Calculate the times of sunrise, sunset, solar noon, twilight and solar position for a given location.
//
// The Arduino Time library is used for timekeeping purposes. In other applications such as a clock, the Timezone 
// library may also be handy for switching back and forth between Daylight saving time and Standard time.
// 
// For more accurate timekeeping, system time should be synchronized with an RTC module and/or NTP server.
//
// Tested with Arduino IDE 1.8.13 and Arduino Uno/Nano
//======================================================================================================================

#include <SolarCalculator.h>
#include <TimeLib.h>

// Location
double latitude = 45.55;
double longitude = -73.633;
int time_zone = -5;

double sunrise;       // Sunrise, in hours (UTC)      
double transit;       // Solar noon, in hours (UTC)
double sunset;        // Sunset, in hours (UTC)
double dawn;          // Civil dawn, in hours (UTC)
double dusk;          // Civil dusk, in hours (UTC)
double eq;            // Equation of Time, in minutes
double rt_ascension;  // Sun's right ascension, in degrees
double declination;   // Sun's declination, in degrees
double azimuth;       // Sun's azimuth, in degrees
double elevation;     // Sun's elevation, in degrees

void setup() 
{
  Serial.begin(9600);

  // Set system time to compile time (UTC)
  setTime(compileTime() - time_zone * 3600);
  
  // Or, set time manually (hr, min, sec, day, mo, yr)
  //setTime(12, 0, 0, 22, 11, 2021);

  // Calculate the times of sunrise, transit and sunset (UTC)
  calcSunriseSunset(year(), month(), day(), latitude, longitude, transit, sunrise, sunset);
  
  // Calculate the times of civil dawn and dusk (UTC)
  calcCivilDawnDusk(year(), month(), day(), latitude, longitude, transit, dawn, dusk);

  // Calculate the equation of time
  calcEquationOfTime(year(), month(), day(), hour(), minute(), second(), eq);

  // Calculate the Sun's right ascension and declination
  calcEquatorialCoordinates(year(), month(), day(), hour(), minute(), second(), rt_ascension, declination);
  
  // To local standard time
  sunrise += time_zone;
  transit += time_zone;
  sunset += time_zone;
  dawn += time_zone;
  dusk += time_zone;  

  // Print results
  Serial.print("Date: ");
  Serial.print(year());
  Serial.print("-");
  Serial.print(month());
  Serial.print("-");
  Serial.println(day());
  
  Serial.print("Latitude: ");
  Serial.print(latitude, 3);
  Serial.print(" Longitude: ");
  Serial.println(longitude, 3);
  Serial.print("UTC offset: ");
  Serial.println(time_zone);
  Serial.println("--"); 

  Serial.print("Sunrise: ");
  printTime24h(sunrise);
  Serial.print("Transit: ");
  printTime24h(transit);
  Serial.print("Sunset:  ");
  printTime24h(sunset);
  Serial.print("Civil dawn: ");
  printTime24h(dawn);
  Serial.print("Civil dusk: ");
  printTime24h(dusk);
  Serial.println("--"); 

  Serial.print("Equation of Time: ");
  Serial.print(eq);
  Serial.println(" minutes");
  Serial.print("Right ascension: ");
  Serial.print(rt_ascension);
  Serial.print(" Declination: ");
  Serial.println(declination);
  Serial.println("--"); 
}

void loop() 
{
  time_t utc = now();
  static time_t last;
  static double last_azimuth;
  static double last_elevation;

  if (utc != last)
  {
    // Calculate the Sun's azimuth and elevation
    calcHorizontalCoordinates(year(utc), month(utc), day(utc), hour(utc), minute(utc), second(utc), latitude, longitude, 
                              azimuth, elevation);

    // Round to two decimal places
    azimuth = double(round(azimuth * 100)) / 100;
    elevation = double(round(elevation * 100)) / 100;

    // Print only if azimuth or elevation has changed
    if (azimuth != last_azimuth || elevation != last_elevation)
    {
      Serial.print("Azimuth: ");
      Serial.print(azimuth);
      Serial.print(" Elevation: ");
      Serial.println(elevation);
      
      last_azimuth = azimuth;
      last_elevation = elevation;
    }
    last = utc;
  }
}

// Code from JChristensen/Timezone Clock example
time_t compileTime()
{
  const uint8_t COMPILE_TIME_DELAY = 0;
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

void printTime24h(double h) 
{
  int m = int(round(h * 60));
  int hours = (m / 60) % 24;
  printDigits(hours);
  Serial.print(":");
  int minutes = m % 60;
  printDigits(minutes);
  Serial.println();
}

void printDigits(int digits)
{
  if (digits < 10)
    Serial.print('0');
  Serial.print(digits);
}
