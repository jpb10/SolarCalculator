//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SolarTrackingTimeLib.ino
//
// Monitor the Sun's position in the sky for any location on Earth.
//
// Tested with Arduino IDE 1.8.19 and Arduino Uno
//======================================================================================================================

#include <SolarCalculator.h>
#include <TimeLib.h>

// Location
double latitude = 45.55;
double longitude = -73.633;
int utc_offset = -5;

// Refresh interval, in seconds
int interval = 10;

void setup()
{
  Serial.begin(9600);

  // Set system time to compile time
  setTime(toUtc(compileTime()));

  // Set time manually (hr, min, sec, day, mo, yr)
  //setTime(0, 0, 0, 1, 1, 2022);
}

void loop()
{
  static unsigned long next_millis = 0;

  // At every interval
  if (millis() > next_millis)
  {
    time_t utc = now();
    double az, el;

    // Calculate the solar position, in degrees
    calcHorizontalCoordinates(utc, latitude, longitude, az, el);

    // Print results
    Serial.print(F("Az: "));
    Serial.print(az);
    Serial.print(F("°  El: "));
    Serial.print(el);
    Serial.println(F("°"));

    next_millis = millis() + interval * 1000L;
  }
}

time_t toUtc(time_t local)
{
  return local - utc_offset * 3600L;
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
