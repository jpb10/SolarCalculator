//======================================================================================================================
// SolarCalculator Library for Arduino
//
// SolarCalculator is based on the NOAA Solar Calculator: https://gml.noaa.gov/grad/solcalc/
//
// This library provides functions to calculate the Sun's position in the sky, the times of sunrise, sunset, twilight
// and solar noon for any location on earth, as well as the equation of time and more.
//
// Most formulae are taken from the textbook Astronomical Algorithms by Jean Meeus or the NOAA Solar Calculator
// source code. Other sources are cited in the comments.
//======================================================================================================================

#ifndef SOLARCALCULATOR_H
#define SOLARCALCULATOR_H

#include <Arduino.h>

const double SUNRISESET_STD_ALTITUDE = -0.8333;
const double CIVIL_DAWNDUSK_STD_ALTITUDE = -6.0;
const double NAUTICAL_DAWNDUSK_STD_ALTITUDE = -12.0;
const double ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE = -18.0;

//======================================================================================================================
// Intermediate calculations
//
// Time T is measured in Julian Centuries (36525 ephemeris days from the epoch J2000.0)
//======================================================================================================================

double wrapTo360(double angle);
double wrapTo180(double angle);
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
double calcObliquityCorrection(double T);
double calcSunRtAscension(double T);
double calcSunDeclination(double T);
double calcNutationRtAscension(double T);
double calcGrMeanSiderealTime(double JD);
double calcGrApparentSiderealTime(double JD);
double calcGrSiderealTimeInstant(double GAST, double m);
double calcSunAzimuth(double HA, double decl, double lat);
double calcSunElevation(double HA, double decl, double lat);
double calcRefraction(double elev);
double equationOfTimeSmart(double T);
double equationOfTimeHughes(double T);
double equationOfTimeMeeus(double T);
double calcDeltaT(double year, double month);
double calcDeltaTPoly(double year, double month);

//======================================================================================================================
// Solar calculator
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//
// Results are passed by reference
//======================================================================================================================

// Calculate the equation of (ephemeris) time, in minutes of time
//
void calcEquationOfTime(int year, int month, int day, int hour, int minute, int second,
                        double &E);

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

// Calculate the Sun's times of rising, transit and setting, in hours
//
void calcSunriseSunset(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset,
                       double altitude = SUNRISESET_STD_ALTITUDE);

// Calculate the times of civil, nautical and astronomical dawn and dusk, in hours
//
void calcCivilDawnDusk(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &dawn, double &dusk);

void calcNauticalDawnDusk(int year, int month, int day, double latitude, double longitude,
                          double &transit, double &dawn, double &dusk);

void calcAstronomicalDawnDusk(int year, int month, int day, double latitude, double longitude,
                              double &transit, double &dawn, double &dusk);

#endif //SOLARCALCULATOR_H
