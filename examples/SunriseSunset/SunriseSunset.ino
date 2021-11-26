//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: SunriseSunset
//  
// Calculate the times of sunrise, sunset, solar noon and twilight for a given date and location.
//
// Tested with Arduino IDE 1.8.13 and Arduino Uno/Nano
//======================================================================================================================

#include <SolarCalculator.h>

// Date
int year = 2021;
int month = 11;
int day = 22;

// Location
double latitude = 45.55;
double longitude = -73.633;
int time_zone = -5;

double sunrise;  // Sunrise, in hours (UTC)      
double transit;  // Solar noon, in hours (UTC)
double sunset;   // Sunset, in hours (UTC)
double dawn;     // Civil dawn, in hours (UTC)
double dusk;     // Civil dusk, in hours (UTC)

void setup() 
{
  Serial.begin(9600);

  // Calculate the times of sunrise, transit and sunset (UTC)
  calcSunriseSunset(year, month, day, latitude, longitude, transit, sunrise, sunset);

  // Calculate the times of civil dawn and dusk (UTC)
  calcCivilDawnDusk(year, month, day, latitude, longitude, transit, dawn, dusk);

  // To local standard time
  sunrise += time_zone;
  transit += time_zone;
  sunset += time_zone;
  dawn += time_zone;
  dusk += time_zone; 

  // Print result
  Serial.print("Date: ");
  Serial.print(year);
  Serial.print("-");
  Serial.print(month);
  Serial.print("-");
  Serial.println(day);
  
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
}

void loop() 
{
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
