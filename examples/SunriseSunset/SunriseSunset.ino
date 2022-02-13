//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SunriseSunset.ino
//  
// Calculate the times of sunrise, sunset and solar noon for a given date and location.
//
// Tested with Arduino IDE 1.8.19 and Arduino Uno
//======================================================================================================================

#include <SolarCalculator.h>

// Date
const int year = 2022;
const int month = 2;
const int day = 11;

// Location
const double latitude = 45.55;
const double longitude = -73.633;
const int time_zone = -5;

void setup() 
{
  Serial.begin(9600);

  double sunrise;  // Sunrise, in hours (UTC)      
  double transit;  // Solar noon, in hours (UTC)
  double sunset;   // Sunset, in hours (UTC)

  // Calculate the times of sunrise, transit and sunset (UTC)
  calcSunriseSunset(year, month, day, latitude, longitude, transit, sunrise, sunset);

  // Print result
  Serial.print("Date: ");
  Serial.print(year);
  Serial.print('-');
  printDigits(month);
  Serial.print('-');
  printDigits(day);
  Serial.println();

  Serial.print("Latitude: ");
  Serial.print(latitude, 3);
  Serial.print(" Longitude: ");
  Serial.println(longitude, 3);
  Serial.print("UTC offset: ");
  Serial.print(time_zone);
  Serial.println();
  Serial.println("--"); 

  Serial.print("Sunrise: ");
  printSunTime24h(toLocal(sunrise));
  Serial.print("Transit: ");
  printSunTime24h(toLocal(transit));
  Serial.print("Sunset:  ");
  printSunTime24h(toLocal(sunset));
}

void loop() 
{
}

double toLocal(double utc)
{
  return utc + time_zone;
}

void printSunTime24h(double h) 
{
  int m = int(round(h * 60));
  int hours = (m / 60) % 24;
  int minutes = m % 60;
  printDigits(hours);
  Serial.print(':');
  printDigits(minutes);
  Serial.println();
}

void printDigits(int digits)
{
  if (digits < 10)
    Serial.print('0');
  Serial.print(digits);
}
