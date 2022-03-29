//======================================================================================================================
// SolarCalculator Library for Arduino
//
// This library provides functions to calculate the times of sunrise, sunset, solar noon and twilight (dawn and dusk),
// solar coordinates, interpolation of coordinates, atmospheric refraction correction, equation of time, delta T, etc.
//
// Most formulae are taken from Astronomical Algorithms by Jean Meeus and adapted for 8-bit AVR platform.
//======================================================================================================================

#ifndef SOLARCALCULATOR_H
#define SOLARCALCULATOR_H

#include <Arduino.h>

//namespace solarcalculator {

constexpr double SUNRISESET_STD_ALTITUDE = -0.8333;
constexpr double CIVIL_DAWNDUSK_STD_ALTITUDE = -6.0;
constexpr double NAUTICAL_DAWNDUSK_STD_ALTITUDE = -12.0;
constexpr double ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE = -18.0;

struct JulianDay
{
    double JD;  // Julian day at 0h UT (JD ending in .5)
    double m;   // Fractional day, 0h to 24h (decimal number between 0 and 1)

    explicit JulianDay(unsigned long utc);  // Unix time, i.e. seconds since 0h UT 1 January 1970
    JulianDay(int year, int month, int day, int hour = 0, int minute = 0, int second = 0);  // Calendar date (UTC)
};

//======================================================================================================================
// Intermediate calculations
//
// Time T is measured in Julian Centuries (36525 ephemeris days from the epoch J2000.0)
//======================================================================================================================

// Utilities
double wrapTo360(double angle);
double wrapTo180(double angle);
double between0And1(double n);
double interpolateCoordinates(double n, double y1, double y2, double y3);

// Julian day and century
double fractionalDay(int hour, int minute, int second);
double calcJulianDay(int year, int month, int day);
double calcJulianCent(double JD, double m = 0);

// Solar coordinates
double calcGeomMeanLongSun(double T);
double calcGeomMeanAnomalySun(double T);
double calcEccentricityEarthOrbit(double T);
double calcSunEqOfCenter(double T);
double calcSunTrueLong(double T);
double calcSunTrueAnomaly(double T);
double calcSunRadVector(double T);
double calcMeanObliquityOfEcliptic(double T);
void calcSolarCoordinates(double T, double &ra, double &dec);

// Sidereal time at Greenwich
double calcGrMeanSiderealTime(double JD, double m = 0);

// Sun's position in the sky, solar time, and ΔT
double calcAzimuth(double H, double delta, double lat);
double calcElevation(double H, double delta, double lat);
double calcRefraction(double elev);
double equationOfTimeSmart(double T);
double calcDeltaT(double year);

//======================================================================================================================
// Solar calculator
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

// Calculate the equation of time, in minutes of time
void calcEquationOfTime(JulianDay jd, double &E);

// Calculate the Sun's right ascension, declination and radius vector (distance), in degrees and AUs
void calcEquatorialCoordinates(JulianDay jd, double &rt_ascension, double &declination, double &radius_vector);

// Calculate the Sun's azimuth and elevation (altitude), corrected for atmospheric refraction, in degrees
void calcHorizontalCoordinates(JulianDay jd, double latitude, double longitude, double &azimuth, double &elevation);

// Calculate the times of sunrise, transit and sunset, in hours
void calcSunriseSunset(JulianDay jd, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset,
                       double altitude = SUNRISESET_STD_ALTITUDE, int iterations = 1);

//======================================================================================================================
// Wrapper functions
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

void calcEquationOfTime(unsigned long utc, double &E);
void calcEquationOfTime(int year, int month, int day, int hour, int minute, int second, double &E);

void calcEquatorialCoordinates(unsigned long utc, double &rt_ascension, double &declination, double &radius_vector);
void calcEquatorialCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double &rt_ascension, double &declination, double &radius_vector);

void calcHorizontalCoordinates(unsigned long utc, double latitude, double longitude,
                               double &azimuth, double &elevation);
void calcHorizontalCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double latitude, double longitude, double &azimuth, double &elevation);

void calcSunriseSunset(unsigned long utc, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset,
                       double altitude = SUNRISESET_STD_ALTITUDE, int iterations = 1);
void calcSunriseSunset(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset,
                       double altitude = SUNRISESET_STD_ALTITUDE, int iterations = 1);

void calcCivilDawnDusk(unsigned long utc, double latitude, double longitude,
                       double &transit, double &dawn, double &dusk);
void calcCivilDawnDusk(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &dawn, double &dusk);

void calcNauticalDawnDusk(unsigned long utc, double latitude, double longitude,
                          double &transit, double &dawn, double &dusk);
void calcNauticalDawnDusk(int year, int month, int day, double latitude, double longitude,
                          double &transit, double &dawn, double &dusk);

void calcAstronomicalDawnDusk(unsigned long utc, double latitude, double longitude,
                              double &transit, double &dawn, double &dusk);
void calcAstronomicalDawnDusk(int year, int month, int day, double latitude, double longitude,
                              double &transit, double &dawn, double &dusk);

//}  // namespace
#endif  //SOLARCALCULATOR_H
