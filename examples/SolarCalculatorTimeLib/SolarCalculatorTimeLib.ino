//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SolarCalculatorTimeLib.ino
//
// Calculate the rise and set times, equation of time, and current solar coordinates.
//
// Tested with Arduino IDE 1.8.19 and Arduino Uno
//======================================================================================================================

#include <SolarCalculator.h>
#include <TimeLib.h>

// Location
double latitude = 45.55;
double longitude = -73.633;
int utc_offset = -5;

void setup()
{
  Serial.begin(9600);

  double transit, sunrise, sunset;  // Event times, in hours (UTC)
  double eq;                        // Equation of time, in minutes
  double ra, dec, r;                // Equatorial coordinates, in degrees and AUs
  double az, el;                    // Horizontal coordinates, in degrees

  // Set system time to compile time
  setTime(toUtc(compileTime()));

  // Set time manually (hr, min, sec, day, mo, yr)
  //setTime(0, 0, 0, 1, 1, 2022);

  // Get current time
  time_t utc = now();

  calcEquationOfTime(utc, eq);
  calcEquatorialCoordinates(utc, ra, dec, r);
  calcHorizontalCoordinates(utc, latitude, longitude, az, el);
  calcSunriseSunset(utc, latitude, longitude, transit, sunrise, sunset);

  // Print results
  Serial.print(F("Sunrise: "));
  printSunTime24h(sunrise + utc_offset);
  Serial.print(F("Transit: "));
  printSunTime24h(transit + utc_offset);
  Serial.print(F("Sunset:  "));
  printSunTime24h(sunset + utc_offset);
  Serial.print(F("Eq of time: "));
  Serial.print(eq);
  Serial.println(F(" min"));
  Serial.print(F("RA: "));
  Serial.print(degreesToHours(ra), 3);
  Serial.print(F("h  Dec: "));
  Serial.print(dec);
  Serial.print(F("°  R: "));
  Serial.print(r, 6);
  Serial.println(F(" AU"));
  Serial.print(F("Az: "));
  Serial.print(az);
  Serial.print(F("°  El: "));
  Serial.print(el);
  Serial.println(F("°"));
}

void loop()
{
}

time_t toUtc(time_t local)
{
  return local - utc_offset * 3600L;
}

double degreesToHours(double deg)
{
  return deg / 15;
}

// Code from JChristensen/Timezone Clock example
time_t compileTime()
{
  const uint8_t COMPILE_TIME_DELAY = 8;
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

void printSunTime24h(double hours)
{
  int m = int(round(hours * 60));
  int hr = (m / 60) % 24;
  int mn = m % 60;
  printDigits(hr);
  Serial.print(':');
  printDigits(mn);
  Serial.println();
}

void printDigits(int digits)
{
  if (digits < 10)
    Serial.print('0');
  Serial.print(digits);
}
