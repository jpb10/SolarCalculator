//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: EquationOfTime.ino
//  
// Plot the equation of time for a given year.
// 
// Tested with Arduino IDE 1.8.19 and Arduino Uno
//======================================================================================================================

#include <SolarCalculator.h>

int year = 2022;

void setup()
{
  Serial.begin(9600);

  // Starting day (January 1)
  double jd = calcJulianDay(year, 1, 1);

  for (int i = 0; i < 365; i++)
  {
    // Calculate the equation of time
    double t = calcJulianCent(jd + i);
    double eq = 4 * equationOfTimeSmart(t);  // convert degrees to minutes of time

    // View with serial plotter (Ctrl+Shift+L)
    Serial.println(eq);
  }
}

void loop() 
{
}
