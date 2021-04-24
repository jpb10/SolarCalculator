#ifndef SOLARCALCULATOR_H
#define SOLARCALCULATOR_H

#include "Arduino.h"

const double SUNRISESET_STD_ALTITUDE = -0.8333;
const double CIVIL_DAWNDUSK_STD_ALTITUDE = -6.0;
const double NAUTICAL_DAWNDUSK_STD_ALTITUDE = -12.0;
const double ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE = -18.0;

//
// Intermediate calculations
// Time T is measured in Julian Centuries (36525 ephemeris days from the epoch J2000.0)
//

double between0And1(double n);
double interpolateCoordinates(double n, double y1, double y2, double y3);
double fractionalDay(int hour, int minute, int second);
double calcJulianDay(int year, int month, int day);
double calcJulianCent(double JD);
double calcJulianCentSplit(double JD, double m);
double calcGeomMeanLongSun(double T);
double calcGeomMeanAnomalySun(double T);
double calcEccentricityEarthOrbit(double T);
double calcSunEqOfCenter(double T);
double calcSunTrueLong(double T);
double calcSunTrueAnomaly(double T);
double calcSunRadVector(double T);
double calcSunApparentLong(double T);
double calcMeanObliquityOfEcliptic(double T);
double calcNutationLongitude(double T);
double calcNutationObliquity(double T);
double calcObliquityCorrection1(double T);
double calcObliquityCorrection2(double T);
double calcSunRtAscension(double T);
double calcSunDeclination(double T);
double calcNutationRtAscension(double T);
double calcMeanSiderealTime(double JD);
double calcApparentSiderealTime(double JD);
double calcSiderealTimeInstant(double GAST, double m);
double calcSolarElevation(double ha, double decl, double lat);
double calcSolarAzimuth(double ha, double decl, double lat);
double calcRefractionCorr(double elev);
double equationOfTime1(double T);
double equationOfTime2(double T);
double equationOfTime3(double T);
double calcDeltaT(double year, double month);

//
// Solar calculator
// Results are passed by reference
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//

// Calculate the Equation of (Ephemeris) Time, in minutes of time
//
// sel = 1: Equation by W.M. Smart, Textbook on Spherical Astronomy (1971) (default)
// sel = 2: Equation by D.W. Hughes, The equation of time - NASA/ADS (1989)
// sel = 3: As defined by Jean Meeus, Astronomical Algorithms (1991)
//
void calcEquationOfTime(int year, int month, int day, int hour, int minute, int second,
                        double &E, int sel = 1);

// Calculate the Sun's right ascension and declination, in degrees
//
void calcEquatorialCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double &rt_ascension, double &declination);

// Calculate the Sun's azimuth and elevation (altitude), in degrees
//
void calcHorizontalCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double latitude, double longitude,
                               double &azimuth, double &elevation);

// Calculate the Sun's radius vector (distance), in AUs
//
void calcSunRadiusVector(int year, int month, int day, int hour, int minute, int second,
                         double &radius_vector);

// Calculate the Sun's times of rising, transit and setting, in fraction of days
//
// local = false:   Results are between 0 and 1, Universal Time.
// local = true:    Results can be less than 0 or greater than 1, Universal Time. (default)
//                  *** Use the default option if you intend to convert your results to local standard time. ***
//
void calcSunriseSunset(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset,
                       double altitude = SUNRISESET_STD_ALTITUDE, bool local = true);

// Calculate the times of civil, nautical and astronomical dawn and dusk, in fraction of days
//
void calcCivilDawnDusk(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &dawn, double &dusk);

void calcNauticalDawnDusk(int year, int month, int day, double latitude, double longitude,
                          double &transit, double &dawn, double &dusk);

void calcAstronomicalDawnDusk(int year, int month, int day, double latitude, double longitude,
                              double &transit, double &dawn, double &dusk);

#endif //SOLARCALCULATOR_H
