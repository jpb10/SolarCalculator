//======================================================================================================================
// SolarCalculator Library for Arduino example sketch: EquationOfTime
//  
// Plot the equation of time for a given year.
// 
//   The equation of time is the difference between apparent time (sundials) and mean time (clocks). This difference is 
//   mainly due to the Earth's elliptic orbit around the Sun and the tilt of the Earth's axis of rotation.
// 
// Julian days are utilized to go through every day of a given year without a calendar.
// 
//   The Julian Day is a continuous count of days from the beginning of the year -4712. By convention, Julian days begin 
//   at noon (whole numbers) and any instant can be expressed as a fraction of a day (decimal numbers).
// 
// Display the curve with the Arduino Serial Plotter (Ctrl+Shift+L).
// 
// Tested with Arduino IDE 1.8.13 and Arduino Uno/Nano
//======================================================================================================================

#include <SolarCalculator.h>

int year = 2021;

// Starting Julian Day (January 1)
double JD = calcJulianDay(year, 1, 1);

void setup() 
{
  Serial.begin(9600);

  for (int i = 0; i < 365; i++) 
  {
    // Calculate the equation of time
    double T = calcJulianCent(JD);
    double E = 4 * equationOfTimeSmart(T);  // Multiply by 4 to convert degrees to minutes of time 

    // Print result
    Serial.println(E);  // Press Ctrl+Shift+L to launch serial plotter

    // Next Julian Day
    ++JD;
  }
}

void loop() 
{
}
