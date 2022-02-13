//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SolarTrackingTimeLib.ino
//
// Monitor the Sun's position in the sky for any location on Earth.
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

// Refresh interval, in seconds
const int interval = 10;

void setup() 
{
  Serial.begin(9600);

  // Set system time to compile time (UTC)
  setTime(toUtc(compileTime()));
  
  // Or, set time manually (hr, min, sec, day, mo, yr)
  //setTime(12, 0, 0, 11, 2, 2022);
  
  // Print results
  Serial.print("Date: ");
  Serial.print(year());
  Serial.print('-');
  printDigits(month());
  Serial.print('-');
  printDigits(day());
  Serial.println();
  
  Serial.print("Latitude: ");
  Serial.print(latitude, 3);
  Serial.print(" Longitude: ");
  Serial.println(longitude, 3);
  Serial.print("UTC offset: ");
  Serial.println(time_zone);
  Serial.println("--");
}

void loop() 
{
  static unsigned long next_millis = 0;

  // At every interval
  if (millis() > next_millis)
  {
    time_t utc = now();
    double azimuth, elevation;
    
    // Calculate the Sun's azimuth and elevation (corrected for atmospheric refraction), in degrees
    calcHorizontalCoordinates(year(utc), month(utc), day(utc), hour(utc), minute(utc), second(utc), latitude, longitude, 
                              azimuth, elevation);      

    Serial.print('(');
    Serial.print(hour(toLocal(utc)));
    Serial.print(':');
    printDigits(minute(toLocal(utc)));
    Serial.print(':');
    printDigits(second(toLocal(utc)));
    Serial.print(") Azimuth: ");
    Serial.print(azimuth);
    Serial.print("° Elevation: ");
    Serial.print(elevation);
    Serial.println("°");

    next_millis = millis() + interval * 1000;
  }
}

time_t toUtc(time_t local)
{
  return local - time_zone * 3600;
}

time_t toLocal(time_t utc)
{
  return utc + time_zone * 3600;
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

void printDigits(int digits)
{
  if (digits < 10)
    Serial.print('0');
  Serial.print(digits);
}
