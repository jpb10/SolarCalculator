//======================================================================================================================
// SolarCalculator Library for Arduino
//
// This library provides functions to calculate the times of sunrise, sunset, solar noon and twilight (dawn and dusk),
// solar coordinates, interpolation of coordinates, atmospheric refraction correction, equation of time, delta T, etc.
//
// Most formulae are taken from Astronomical Algorithms by Jean Meeus and adapted for 8-bit AVR platform.
//======================================================================================================================

#include "SolarCalculator.h"

//namespace solarcalculator {

JulianDay::JulianDay(int year, int month, int day, int hour, int minute, int second)
{
    JD = calcJulianDay(year, month, day);
    m = fractionalDay(hour, minute, second);
}

JulianDay::JulianDay(unsigned long utc)
{
    JD = static_cast<unsigned long>(utc / 86400) + 2440587.5;
    m = (utc % 86400) / 86400.0;
}

//======================================================================================================================
// Intermediate calculations
//
// Time T is measured in Julian centuries (36525 ephemeris days from the epoch J2000.0)
//======================================================================================================================

double wrapTo360(double angle)
{
    angle = fmod(angle, 360);
    if (angle < 0) angle += 360;
    return angle;  // [0, 360)
}

double wrapTo180(double angle)
{
    angle = fmod(angle + 180, 360);
    if (angle < 0) angle += 360;
    return angle - 180;  // [-180, 180)
}

double between0And1(double n)
{
    while (n < 0) ++n;
    while (n > 1) --n;
    return n;
}

double interpolateCoordinates(double n, double y1, double y2, double y3)
{
    if (fabs(y2 - y1) > 180)  // if coordinate is discontinuous
    {                         // add or subtract 360 degrees
        if (y1 < 0) y1 += 360;
        else if (y2 < 0) y1 -= 360;
    }
    else if (fabs(y3 - y2) > 180)
    {
        if (y3 < 0) y3 += 360;
        else if (y2 < 0) y3 -= 360;
    }
    double a = y2 - y1;
    double b = y3 - y2;
    double c = b - a;
    return y2 + n * (a + b + n * c) / 2;
}

double fractionalDay(int hour, int minute, int second)
{
    return (hour + minute / 60.0 + second / 3600.0) / 24;
}

double calcJulianDay(int year, int month, int day)
{
//    if (month <= 2)
//    {
//        year -= 1;
//        month += 12;
//    }
//    double A = floor(year / 100.0);
//    double B = 2 - A + floor(A / 4);
//    return floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5;

    // Valid only from 1901 to 2099, Van Flandern & Pulkkinen (1979)
    return 367.0 * year - static_cast<int>(7 * (year + (month + 9) / 12) / 4) +
                          static_cast<int>(275 * month / 9) + day + 1721013.5;
}

double calcJulianCent(double JD, double m)
{
    return (JD - 2451545 + m) / 36525;
}

double calcGeomMeanLongSun(double T)
{
    return wrapTo360(280.46646 + T * 36000.76983);  // in degrees
}

double calcGeomMeanAnomalySun(double T)
{
    return wrapTo360(357.52911 + T * 35999.05029);  // in degrees
}

double calcEccentricityEarthOrbit(double T)
{
    return 0.016708634 - T * 0.000042037;  // no units
}

double calcSunEqOfCenter(double T)
{
    double M = calcGeomMeanAnomalySun(T);
    return sin(radians(M)) * 1.914602 + sin(2 * radians(M)) * 0.019993;  // in degrees
}

double calcSunTrueLong(double T)
{
    double L0 = calcGeomMeanLongSun(T);
    double C = calcSunEqOfCenter(T);
    return L0 + C;  // in degrees
}

double calcSunTrueAnomaly(double T)
{
    double M = calcGeomMeanAnomalySun(T);
    double C = calcSunEqOfCenter(T);
    return M + C;  // in degrees
}

double calcSunRadVector(double T)
{
//    double v = calcSunTrueAnomaly(T);
//    double e = calcEccentricityEarthOrbit(T);
//    return 1.000001018 * (1 - e * e) / (1 + e * cos(radians(v)));
    double M = calcGeomMeanAnomalySun(T);
    return 1.00014 - 0.01671 * cos(radians(M)) - 0.00014 * cos(2 * radians(M));  // in AUs
}

double calcMeanObliquityOfEcliptic(double T)
{
    return 23.4392911 - T * 0.0130042;  // in degrees
}

// Accurate to 0.01 degree (valid from 1950 to 2050)
void calcSolarCoordinates(double T, double &ra, double &dec)
{
    double L = calcSunTrueLong(T) - 0.00569;
    double epsilon0 = calcMeanObliquityOfEcliptic(T);
    ra = degrees(atan2(cos(radians(epsilon0)) * sin(radians(L)), cos(radians(L))));  // [-180, 180)
    dec = degrees(asin(sin(radians(epsilon0)) * sin(radians(L))));
}

// Valid only for JD at 0h UT, Greenwich (JD ending in .5)
double calcGrMeanSiderealTime(double JD, double m)
{
    double T = calcJulianCent(JD);
    double GMST = wrapTo360(100.46061837 + T * 36000.770053608);
    return wrapTo360(GMST + 360.985647 * m);  // in degrees
}

double calcAzimuth(double H, double delta, double lat)
{
    return degrees(atan2(sin(radians(H)), cos(radians(H)) * sin(radians(lat)) -
                   tan(radians(delta)) * cos(radians(lat))));  // in degrees
}

double calcElevation(double H, double delta, double lat)
{
    return degrees(asin(sin(radians(lat)) * sin(radians(delta)) +
                   cos(radians(lat)) * cos(radians(delta)) * cos(radians(H))));  // in degrees
}

// Equation of ephemeris time by Smart (1978)
double equationOfTimeSmart(double T)
{
    double L0 = calcGeomMeanLongSun(T);
    double M = calcGeomMeanAnomalySun(T);
    return 2.465 * sin(2 * radians(L0)) - 1.913 * sin(radians(M)) +
           0.165 * sin(radians(M)) * cos(2 * radians(L0)) -
           0.053 * sin(4 * radians(L0)) - 0.02 * sin(2 * radians(M));  // in degrees
}

// Approximate atmospheric refraction correction, in degrees
double calcRefraction(double elev)
{
    if (elev < -0.575)
        return -20.774 / tan(radians(elev)) / 3600;  // Zimmerman (1981)
    else
        return 1.02 / tan(radians(elev + 10.3 / (elev + 5.11))) / 60;  // Sæmundsson (1986)
}

// Simple polynomial expressions for delta T (ΔT), in seconds of time
double calcDeltaT(double year)
{
    if (year > 1997)
    {
        double t = year - 2015;
        return 67.62 + t * (0.3645 + 0.0039755 * t);  // Fred Espenak (2014)
    }
    else  // y > 948
    {
        double u = (year - 2000) / 100;
        return 64.69 + u * (80.59 + 23.604 * u);  // fitted to historical data, very approximate before 1900
    }
}

//======================================================================================================================
// Solar calculator
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

// Calculate the equation of time, in minutes of time
void calcEquationOfTime(JulianDay day, double &E)
{
    double T = calcJulianCent(day.JD, day.m);
    E = 4 * equationOfTimeSmart(T);
}

// Calculate the Sun's right ascension, declination and radius vector (distance), in degrees and AUs
void calcEquatorialCoordinates(JulianDay day, double &rt_ascension, double &declination, double &radius_vector)
{
    double T = calcJulianCent(day.JD, day.m);
    calcSolarCoordinates(T, rt_ascension, declination);  // [-180, 180)
    radius_vector = calcSunRadVector(T);
}

// Calculate the Sun's azimuth and elevation (altitude), in degrees
void calcHorizontalCoordinates(JulianDay day, double latitude, double longitude, double &azimuth, double &elevation)
{
    double JD = day.JD;
    double m = day.m;
    double T = calcJulianCent(JD, m);
    double theta0 = calcGrMeanSiderealTime(JD, m);

    double alpha, delta;
    calcSolarCoordinates(T, alpha, delta);
    double H = theta0 + longitude - alpha;

    azimuth = 180 + calcAzimuth(H, delta, latitude);  // measured from the North
    elevation = calcElevation(H, delta, latitude);
    elevation += calcRefraction(elevation);
}

// Sunrise equation helper
double sunriseEquation(double JD, double &m, double latitude, double longitude, double h0)
{
    double T = calcJulianCent(JD, m);
    double theta0 = calcGrMeanSiderealTime(JD, m);

    double alpha, delta;
    calcSolarCoordinates(T, alpha, delta);

    // Approximate time of transit for instant m
    m = between0And1((alpha - longitude - theta0) / 360 + m);

    // Local hour angle at sunrise or sunset (±NaN if body is circumpolar)
    return degrees(acos((sin(radians(h0)) - sin(radians(latitude)) * sin(radians(delta))) /
                   (cos(radians(latitude)) * cos(radians(delta)))));
}

// Calculate the times of sunrise, transit and sunset, in hours
void calcSunriseSunset(JulianDay day, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset, double altitude, int iterations)
{
    double JD = day.JD;
    double m0 = 0.5 - longitude / 360;

    // Transit
    double H0 = sunriseEquation(JD, m0, latitude, longitude, altitude);

    // Approximate times of sunrise and sunset
    double m1 = m0 - H0 / 360;
    double m2 = m0 + H0 / 360;

    // Correction of sunrise and sunset
    if (iterations)
    {
        for (int i = 0; i < iterations; i++)
        {
            H0 = sunriseEquation(JD, m1, latitude, longitude, altitude);
            m1 = m1 - H0 / 360;

            H0 = sunriseEquation(JD, m2, latitude, longitude, altitude);
            m2 = m2 + H0 / 360;
        }
    }

    transit = m0 * 24;
    sunrise = m1 * 24;
    sunset = m2 * 24;
}

//======================================================================================================================
// Wrapper functions
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

void calcEquationOfTime(unsigned long utc, double &E)
{
    JulianDay JD(utc);
    calcEquationOfTime(JD, E);
}

void calcEquationOfTime(int year, int month, int day, int hour, int minute, int second, double &E)
{
    JulianDay JD(year, month, day, hour, minute, second);
    calcEquationOfTime(JD, E);
}

void calcEquatorialCoordinates(unsigned long utc, double &rt_ascension, double &declination, double &radius_vector)
{
    JulianDay JD(utc);
    calcEquatorialCoordinates(JD, rt_ascension, declination, radius_vector);
}

void calcEquatorialCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double &rt_ascension, double &declination, double &radius_vector)
{
    JulianDay JD(year, month, day, hour, minute, second);
    calcEquatorialCoordinates(JD, rt_ascension, declination, radius_vector);
}

void calcHorizontalCoordinates(unsigned long utc, double latitude, double longitude,
                               double &azimuth, double &elevation)
{
    JulianDay JD(utc);
    calcHorizontalCoordinates(JD, latitude, longitude, azimuth, elevation);
}

void calcHorizontalCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double latitude, double longitude, double &azimuth, double &elevation)
{
    JulianDay JD(year, month, day, hour, minute, second);
    calcHorizontalCoordinates(JD, latitude, longitude, azimuth, elevation);
}

void calcSunriseSunset(unsigned long utc, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset, double altitude, int iterations)
{
    JulianDay JD(utc);
    calcSunriseSunset(JD, latitude, longitude, transit, sunrise, sunset, altitude, iterations);
}

void calcSunriseSunset(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset, double altitude, int iterations)
{
    JulianDay JD(year, month, day);
    calcSunriseSunset(JD, latitude, longitude, transit, sunrise, sunset, altitude, iterations);
}

void calcCivilDawnDusk(unsigned long utc, double latitude, double longitude,
                       double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(utc, latitude, longitude, transit, dawn, dusk, CIVIL_DAWNDUSK_STD_ALTITUDE);
}

void calcCivilDawnDusk(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(year, month, day, latitude, longitude, transit, dawn, dusk, CIVIL_DAWNDUSK_STD_ALTITUDE);
}

void calcNauticalDawnDusk(unsigned long utc, double latitude, double longitude,
                          double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(utc, latitude, longitude, transit, dawn, dusk, NAUTICAL_DAWNDUSK_STD_ALTITUDE);
}

void calcNauticalDawnDusk(int year, int month, int day, double latitude, double longitude,
                          double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(year, month, day, latitude, longitude, transit, dawn, dusk, NAUTICAL_DAWNDUSK_STD_ALTITUDE);
}

void calcAstronomicalDawnDusk(unsigned long utc, double latitude, double longitude,
                              double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(utc, latitude, longitude, transit, dawn, dusk, ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE);
}

void calcAstronomicalDawnDusk(int year, int month, int day, double latitude, double longitude,
                              double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(year, month, day, latitude, longitude, transit, dawn, dusk, ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE);
}

//}  // namespace
