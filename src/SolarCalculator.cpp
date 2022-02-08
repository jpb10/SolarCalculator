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

#include "SolarCalculator.h"

//======================================================================================================================
// Intermediate calculations
//
// Time T is measured in Julian Centuries (36525 ephemeris days from the epoch J2000.0)
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
    if (month <= 2)
    {
        year -= 1;
        month += 12;
    }
    double A = floor(year / 100.0);
    double B = 2 - A + floor(A / 4);
    return floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5;
}

double calcJulianCent(double JD)
{
    return (JD - 2451545) / 36525;
}

double calcJulianCentSplit(double JD, double m)
{
    return calcJulianCent(JD) + m / 36525;
}

double calcGeomMeanLongSun(double T)
{
    return wrapTo360(280.46646 + T * (36000.76983 + 0.0003032 * T));  // in degrees
}

double calcGeomMeanAnomalySun(double T)
{
    return wrapTo360(357.52911 + T * (35999.05029 - 0.0001537 * T));  // in degrees
}

double calcEccentricityEarthOrbit(double T)
{
    return 0.016708634 - T * (0.000042037 + 0.0000001267 * T);  // no units
}

double calcSunEqOfCenter(double T)
{
    double M = radians(calcGeomMeanAnomalySun(T));
    return sin(M) * (1.914602 - T * (0.004817 + 0.000014 * T)) + sin(2 * M) * (0.019993 - 0.000101 * T) +
           sin(3 * M) * 0.000289;  // in degrees
}

double calcSunTrueLong(double T)
{
    return calcGeomMeanLongSun(T) + calcSunEqOfCenter(T);  // in degrees
}

double calcSunTrueAnomaly(double T)
{
    return calcGeomMeanAnomalySun(T) + calcSunEqOfCenter(T);  // in degrees
}

double calcSunRadVector(double T)
{
    double v = radians(calcSunTrueAnomaly(T));
    double e = calcEccentricityEarthOrbit(T);
    return 1.000001018 * (1 - e * e) / (1 + e * cos(v));  // in AUs
}

/*
// These expressions require double-precision arithmetic
//
double calcNutationLongitude(double T)
{
    double omega = radians(125.04452 - 1934.136261 * T);
    double L1 = radians(280.4665 + 36000.7698 * T);
    double L2 = radians(218.3165 + 481267.8813 * T);
    return (-17.2 * sin(omega) - 1.32 * sin(2 * L1) - 0.23 * sin(2 * L2) + 0.21 * sin(2 * omega)) / 3600;  // in degrees
}

double calcNutationObliquity(double T)
{
    double omega = radians(125.04452 - 1934.136261 * T);
    double L1 = radians(280.4665 + 36000.7698 * T);
    double L2 = radians(218.3165 + 481267.8813 * T);
    return (9.20 * cos(omega) + 0.57 * cos(2 * L1) + 0.10 * cos(2 * L2) - 0.09 * cos(2 * omega)) / 3600;  // in degrees
}
*/

double calcNutationLongitude(double T)
{
    double omega = radians(125.04 - 1934.136 * T);
    return -0.00478 * sin(omega);  // in degrees
}

double calcNutationObliquity(double T)
{
    double omega = radians(125.04 - 1934.136 * T);
    return 0.00256 * cos(omega);  // in degrees
}

double calcSunApparentLong(double T)
{
    return calcSunTrueLong(T) - 0.00569 + calcNutationLongitude(T);  // in degrees
}

double calcMeanObliquityOfEcliptic(double T)
{
    return 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - 0.001813 * T))) / 60) / 60;  // in degrees
}

double calcObliquityCorrection(double T)
{
    double epsilon0 = calcMeanObliquityOfEcliptic(T);
    double delta_epsilon = calcNutationObliquity(T);
    return epsilon0 + delta_epsilon;  // in degrees
}

double calcSunRtAscension(double T)
{
    double epsilon = radians(calcObliquityCorrection(T));
    double lambda = radians(calcSunApparentLong(T));
    return degrees(atan2(cos(epsilon) * sin(lambda), cos(lambda)));  // in degrees
}

double calcSunDeclination(double T)
{
    double epsilon = radians(calcObliquityCorrection(T));
    double lambda = radians(calcSunApparentLong(T));
    return degrees(asin(sin(epsilon) * sin(lambda)));  // in degrees
}

double calcNutationRtAscension(double T)
{
    double delta_psi = calcNutationLongitude(T);
    double epsilon = radians(calcObliquityCorrection(T));
    return delta_psi * cos(epsilon);  // in degrees
}

// Valid only at 0h UT, Greenwich (JD ending in .5)
double calcGrMeanSiderealTime(double JD)
{
    double T = calcJulianCent(JD);
    return wrapTo360(100.46061837 + T * (36000.770053608 + T * (0.000387933 - T / 38710000)));  // in degrees
}

// We neglect the very small variation of Δψ during time ΔT
double calcGrApparentSiderealTime(double JD)
{
    double T = calcJulianCent(JD);
    return calcGrMeanSiderealTime(JD) + calcNutationRtAscension(T);  // in degrees
}

// We neglect the small variation of Δψ during time m
double calcGrSiderealTimeInstant(double GAST, double m)
{
    return wrapTo360(GAST + 360.985647 * m);  // in degrees
}

/*
// Alternate implementation
//
double calcGrMeanSiderealTimeInstant(double JD, double m)
{
    double theta0 = calcGrMeanSiderealTime(JD);
    return wrapTo360(theta0 + 360.985647 * m);  // in degrees
}

double calcGrApparentSiderealTimeInstant(double JD, double m)
{
    double T = calcJulianCentSplit(JD, m);
    return calcGrMeanSiderealTimeInstant(JD, m) + calcNutationRtAscension(T);  // in degrees
}
*/

double calcSunAzimuth(double HA, double decl, double lat)
{
    return degrees(atan2(sin(radians(HA)), cos(radians(HA)) * sin(radians(lat)) -
                   tan(radians(decl)) * cos(radians(lat))));  // in degrees
}

double calcSunElevation(double HA, double decl, double lat)
{
    return degrees(asin(sin(radians(lat)) * sin(radians(decl)) +
                   cos(radians(lat)) * cos(radians(decl)) * cos(radians(HA))));  // in degrees
}

// Approximate atmospheric refraction correction
// National Oceanic and Atmospheric Administration (NOAA)
//
double calcRefraction(double elev)
{
    if (elev > 85.0)
        return 0.0;
    else
    {
        double te = tan(radians(elev));
        if (elev > 5.0)
            return (58.1 / te - 0.07 / (te * te * te) + 0.000086 / (te * te * te * te * te)) / 3600;
        else if (elev > -0.575)
            return (1735.0 + elev * (-518.2 + elev * (103.4 + elev * (-12.79 + 0.711 * elev)))) / 3600;
        else
            return (-20.774 / te) / 3600;  // in degrees
    }
}

// The equation of (ephemeris) time by W.M. Smart
// Textbook on Spherical Astronomy (1971)
//
double equationOfTimeSmart(double T)
{
    double e = calcEccentricityEarthOrbit(T);
    double L = radians(calcGeomMeanLongSun(T));
    double M = radians(calcGeomMeanAnomalySun(T));
    double epsilon = radians(calcObliquityCorrection(T));
    double y = tan(epsilon / 2) * tan(epsilon / 2);
    return degrees(y * sin(2 * L) - 2 * e * sin(M) + 4 * e * y * sin(M) * cos(2 * L) - 0.5 * y * y * sin(4 * L) -
                   1.25 * e * e * sin(2 * M));  // in degrees
}

// The equation of (ephemeris) time by D.W. Hughes
// The Equation of Time - NASA/ADS (1989)
//
double equationOfTimeHughes(double T)
{
    double e = calcEccentricityEarthOrbit(T);
    double L = radians(calcGeomMeanLongSun(T));
    double M = radians(calcGeomMeanAnomalySun(T));
    double epsilon = radians(calcObliquityCorrection(T));
    double y = tan(epsilon / 2) * tan(epsilon / 2);
    return degrees(0.00000447 * T + 0.00000149 * T * T - 2 * e * sin(M) - 1.25 * e * e * sin(2 * M) +
                   y * sin(2 * L) - 0.5 * y * y * sin (4 * L) + 4 * e * y * sin(M) * cos(2 * L) +
                   2.5 * e * e * y * sin(2 * M) * cos(2 * L) - 4 * e * y * y * sin(M) * cos(4 * L) -
                   8 * e * e * y * sin(M) * sin(M) * sin(2 * L));  // in degrees
}

// The equation of time as defined by Jean Meeus
// Astronomical Algorithms (1991)
//
double equationOfTimeMeeus(double T)
{
    double t = T / 10;
    double L0 = wrapTo360(280.4664567 + t * (360007.6982779 + t * (0.03032028 + t * (1 / 49931.0 - t * (1 / 15300.0 +
                          t / 2000000)))));
    return wrapTo180(L0 - 0.0057183 - calcSunRtAscension(T) - calcNutationRtAscension(T));  // in degrees
}

// Simple polynomial expressions for delta T (ΔT)
// Long-term parabolas fitted to historical data, very approximate before 1900
//
double calcDeltaT(double year, double month)
{
    double y = year + (month - 0.5) / 12;
    if (y > 1997)
    {
        double t = y - 2015;
        return 67.62 + t * (0.3645 + 0.0039755 * t);
    }
    else if (y > 948)
    {
        double u = (y - 2000) / 100;
        return 64.69 + u * (80.59 + 23.604 * u);
    }
    else  // y < 948
    {
        double u = (y - 2000) / 100;
        return 2177 + u * (497 + 44.1 * u);  // in seconds of time
    }
}

// Polynomial expressions for delta T (ΔT) by Fred Espenak
// Valid from year -1999 to +3000
//
double calcDeltaTPoly(double year, double month)
{
    double y = year + (month - 0.5) / 12;
    if (y > 2015)
    {
        double t = y - 2015;
        return 67.62 + t * (0.3645 + 0.0039755 * t);
    }
    else if (y > 2005)
    {
        double t = y - 2005;
        return 64.69 + 0.2930 * t;
    }
    else if (y > 1986)
    {
        double t = y - 2000;
        return 63.86 + t * (0.3345 + t * (-0.060374 + t * (0.0017275 + t * (0.000651814 + 0.00002373599 * t))));
    }
    else if (y > 1961)
    {
        double t = y - 1975;
        return 45.45 + t * (1.067 - t * (1 / 260.0 + t / 718));
    }
    else if (y > 1941)
    {
        double t = y - 1950;
        return 29.07 + t * (0.407 - t * (1 / 233.0 - t / 2547));
    }
    else if (y > 1920)
    {
        double t = y - 1920;
        return 21.20 + t * (0.84493 + t * (-0.076100 + 0.0020936 * t));
    }
    else if (y > 1900)
    {
        double t = y - 1900;
        return -2.79 + t * (1.494119 + t * (-0.0598939 + t * (0.0061966 - 0.000197 * t)));
    }
    else if (y > 1860)
    {
        double t = y - 1860;
        return 7.62 + t * (0.5737 + t * (-0.251754 + t * (0.01680668 + t * (-0.0004473624 + t / 233174))));
    }
    else if (y > 1800)
    {
        double t = y - 1800;
        return 13.72 + t * (-0.332447 + t * (0.0068612 + t * (0.0041116 + t * (-0.00037436 + t * (0.0000121272 +
               t * (-0.0000001699 + 0.000000000875 * t))))));
    }
    else if (y > 1700)
    {
        double t = y - 1700;
        return 8.83 + t * (0.1603 + t * (-0.0059285 + t * (0.00013336 - t / 1174000)));
    }
    else if (y > 1600)
    {
        double t = y - 1600;
        return 120 + t * (-0.9808 + t * (-0.01532 + t / 7129));
    }
    else if (y > 500)
    {
        double u = (y - 1000) / 100;
        return 1574.2 + u * (-556.01 + u * (71.23472 + u * (0.319781 + u * (-0.8503463 + u * (-0.005050998 +
               0.0083572073 * u)))));
    }
    else if (y > -500)
    {
        double u = y / 100;
        return 10583.6 + u * (-1014.41 + u * (33.78311 + u * (-5.952053 + u * (-0.1798452 + u * (0.022174192 +
               0.0090316521 * u)))));
    }
    else  // y < -500
    {
        double u = (y - 1820) / 100;
        return -20 + 32 * u * u;  // in seconds of time
    }
}

//======================================================================================================================
// Solar calculator
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//
// Results are passed by reference
//======================================================================================================================

// Calculate the equation of (ephemeris) time, in minutes of time
//
void calcEquationOfTime(int year, int month, int day, int hour, int minute, int second, double &E)
{
    // Find Julian day
    double JD = calcJulianDay(year, month, day);
    double m = fractionalDay(hour, minute, second);
    double T = calcJulianCentSplit(JD, m);

    // Write result, in minutes of time
    E = 4 * equationOfTimeSmart(T);
}

// Calculate the Sun's right ascension and declination, in degrees
//
void calcEquatorialCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double &rt_ascension, double &declination)
{
    // Find Julian day
    double JD = calcJulianDay(year, month, day);
    double m = fractionalDay(hour, minute, second);

    // Convert to Dynamical Time (TD)
    double delta_t = calcDeltaT(year, month);
    double T = calcJulianCentSplit(JD, m + delta_t / 86400);

    // Write results, in degrees
    rt_ascension = calcSunRtAscension(T);
    declination = calcSunDeclination(T);
}

// Calculate the Sun's azimuth and elevation (altitude), in degrees
//
void calcHorizontalCoordinates(int year, int month, int day, int hour, int minute, int second,
                               double latitude, double longitude, double &azimuth, double &elevation)
{
    // Find Julian day
    double JD = calcJulianDay(year, month, day);
    double m = fractionalDay(hour, minute, second);

    // Convert to Dynamical Time (TD)
    double delta_t = calcDeltaT(year, month);
    double T = calcJulianCentSplit(JD, m + delta_t / 86400);

    // Find solar coordinates
    double alpha = calcSunRtAscension(T);
    double delta = calcSunDeclination(T);

    // Find apparent sidereal time at Greenwich
    double theta0 = calcGrApparentSiderealTime(JD);
    theta0 = calcGrSiderealTimeInstant(theta0, m);

    // Find local angle hour
    double H = theta0 + longitude - alpha;

    // Write results, in degrees
    azimuth = 180 + calcSunAzimuth(H, delta, latitude);
    elevation = calcSunElevation(H, delta, latitude);
    elevation += calcRefraction(elevation);
}

// Calculate the Sun's radius vector (distance), in AUs
//
void calcSunRadiusVector(int year, int month, int day, int hour, int minute, int second, double &radius_vector)
{
    // Find Julian day
    double JD = calcJulianDay(year, month, day);
    double m = fractionalDay(hour, minute, second);

    // Convert to Dynamical Time (TD)
    double delta_t = calcDeltaT(year, month);
    double T = calcJulianCentSplit(JD, m + delta_t / 86400);

    // Write result, in AUs
    radius_vector = calcSunRadVector(T);
}

// Calculate the Sun's times of rising, transit and setting, in hours
//
void calcSunriseSunset(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &sunrise, double &sunset, double altitude)
{
    // Find Julian days
    double JD = calcJulianDay(year, month, day);
    double T1 = calcJulianCent(JD - 1);
    double T2 = calcJulianCent(JD);
    double T3 = calcJulianCent(JD + 1);

    // Standard altitude of the Sun
    double h0 = altitude;

    // Find solar coordinates
    double alpha1 = calcSunRtAscension(T1);
    double delta1 = calcSunDeclination(T1);
    double alpha2 = calcSunRtAscension(T2);
    double delta2 = calcSunDeclination(T2);
    double alpha3 = calcSunRtAscension(T3);
    double delta3 = calcSunDeclination(T3);

    // Find apparent sidereal time at Greenwich
    double theta0 = calcGrApparentSiderealTime(JD);

    // Local angle hour at sunrise or sunset (NaN if body is circumpolar)
    double H0 = degrees(acos((sin(radians(h0)) - sin(radians(latitude)) * sin(radians(delta2))) /
                        (cos(radians(latitude)) * cos(radians(delta2)))));

    // Approximate times of transit, rising and setting
    double m0 = between0And1((alpha2 - longitude - theta0) / 360);
    double m1 = m0 - H0 / 360;
    double m2 = m0 + H0 / 360;

    // Find ΔT = TD - UT
    double delta_t = calcDeltaT(year, month);

    // Number of iterations (Newton's method)
    int max_iterations = 3;
    for (int i = 0; i < max_iterations; i++)
    {
        double n0 = m0 + delta_t / 86400;
        double transit_alpha = interpolateCoordinates(n0, alpha1, alpha2, alpha3);
        double transit_theta0 = calcGrSiderealTimeInstant(theta0, m0);
        double transit_ha = wrapTo180(transit_theta0 + longitude - transit_alpha);
        double transit_corr = -transit_ha / 360;

        double n1 = m1 + delta_t / 86400;
        double rise_alpha = interpolateCoordinates(n1, alpha1, alpha2, alpha3);
        double rise_delta = interpolateCoordinates(n1, delta1, delta2, delta3);
        double rise_theta0 = calcGrSiderealTimeInstant(theta0, m1);
        double rise_ha = rise_theta0 + longitude - rise_alpha;
        double rise_elev = calcSunElevation(rise_ha, rise_delta, latitude);
        double rise_corr = (rise_elev - h0) /
                           (360.0 * cos(radians(rise_delta)) * cos(radians(latitude)) * sin(radians(rise_ha)));

        double n2 = m2 + delta_t / 86400;
        double set_alpha = interpolateCoordinates(n2, alpha1, alpha2, alpha3);
        double set_delta = interpolateCoordinates(n2, delta1, delta2, delta3);
        double set_theta0 = calcGrSiderealTimeInstant(theta0, m2);
        double set_ha = set_theta0 + longitude - set_alpha;
        double set_elev = calcSunElevation(set_ha, set_delta, latitude);
        double set_corr = (set_elev - h0) /
                          (360.0 * cos(radians(set_delta)) * cos(radians(latitude)) * sin(radians(set_ha)));

        m0 += transit_corr;
        m1 += rise_corr;
        m2 += set_corr;
    }

    // Write final times, in hours
    transit = m0 * 24;
    sunrise = m1 * 24;
    sunset = m2 * 24;
}

// Calculate the times of civil, nautical and astronomical dawn and dusk, in hours
//
void calcCivilDawnDusk(int year, int month, int day, double latitude, double longitude,
                       double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(year, month, day, latitude, longitude,
                      transit, dawn, dusk, CIVIL_DAWNDUSK_STD_ALTITUDE);
}

void calcNauticalDawnDusk(int year, int month, int day, double latitude, double longitude,
                          double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(year, month, day, latitude, longitude,
                      transit, dawn, dusk, NAUTICAL_DAWNDUSK_STD_ALTITUDE);
}

void calcAstronomicalDawnDusk(int year, int month, int day, double latitude, double longitude,
                              double &transit, double &dawn, double &dusk)
{
    calcSunriseSunset(year, month, day, latitude, longitude,
                      transit, dawn, dusk, ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE);
}
