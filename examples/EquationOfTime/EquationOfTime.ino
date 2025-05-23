//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: EquationOfTime.ino
//  
// Plot the equation of time for a given year.
// 
// Tested with Arduino IDE 2.3.5 and Arduino Uno
//======================================================================================================================

#include <SolarCalculator.h>

int year = 2000;

void setup()
{
  Serial.begin(9600);

  // Starting day (January 1)
  JulianDay day(year, 1, 1);

  for (int i = 0; i < 365; ++i)
  {
    double eq;
    calcEquationOfTime(day, eq);

    // Print and view with serial plotter
    Serial.println(eq);

    // Next day
    ++day.JD;
  }
}

void loop() 
{
}
