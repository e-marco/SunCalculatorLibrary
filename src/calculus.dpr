{
Updated version July 2019. 
Forked from the original repository https://github.com/ISFH/SunCalculatorLibrary

Copyright (c) 2015 Institute for Solar Energy Research Hamelin
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
 3. Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.
 4. Redistributions of any form whatsoever must retain the following
    acknowledgment: 'This product includes the SunCalculator-Library developed 
	by Marco Ernst and Hendrik Holst at the 
	Institute for Solar Energy Research Hamelin (www.isfh.de)'
	
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
}
library calculus;

{ Important note about DLL memory management: ShareMem must be the
  first unit in your library's USES clause AND your project's (select
  Project-View Source) USES clause if your DLL exports any procedures or
  functions that pass strings as parameters or function results. This
  applies to all strings passed to and from your DLL--even those that
  are nested in records and classes. ShareMem is the interface unit to
  the BORLNDMM.DLL shared memory manager, which must be deployed along
  with your DLL. To avoid using BORLNDMM.DLL, pass string information
  using PChar or ShortString parameters. }

uses
  SysUtils,
  Classes,
  DateUtils,
  Math, coordinateTypes;

{$R *.res}

const
  SUNS_ORBIT_EXCENTRICITY = 0.0167; // excentricity  of sun's orbit

// Calculation of the Julian Day for a given Date+Time
function JulianDay(UniversalTime : TDateTime):double; STDCALL;
var
MON, YEA, DAY, HOUR, MIN, SEC, mSEC : Word;
A,B,H,GGG,JUL,S : double;
begin
  DecodeDateTime(UniversalTime, YEA, MON, DAY, HOUR, MIN, SEC, mSEC);

  H := (HOUR + (MIN / 60) + (SEC/3600))/24;      // partially fractured day
  GGG := 1;
  if (YEA <= 1585) then GGG := 0;
  JUL := -1 * floor(7 * (floor((MON + 9) / 12) + YEA) / 4);
  S := 1;
  if ((MON - 9)<0) then S:=-1;
  A := abs(MON - 9);
  B := floor(YEA + S * floor(A / 7));
  B := -1 * floor((floor(B / 100) + 1) * 3 / 4);
  JUL := JUL + floor(275 * MON / 9) + DAY + (GGG * B);
  JUL := JUL + 1721027 + 2 * GGG + 367 * YEA - 0.5;
  JUL := JUL + H ;
  // Calculation from
  //http://www.physics.sfasu.edu/astro/javascript/julianday.html
  result:=JUL;
end;

// Correction function with the equation of Time
function EquationOfTime(JUL: double):double; STDCALL;
var  numberOfDaysSinceJ2000, L, g, ecliptic, tau, REC, lambda : double;
begin
  numberOfDaysSinceJ2000:=JUL - 2451545; // Time coordinate n, Number of days since the
                                        // Standard equinox J2000.0 12:00
  L := 280.460 + 0.9856474 * numberOfDaysSinceJ2000; // middle ecliptic length L of the sun
  g := 357.528 + 0.9856003 * numberOfDaysSinceJ2000; // middle anomaly g
  ecliptic := 23.439 - 0.0000004 * numberOfDaysSinceJ2000;    // Inclined one of the ecliptic

  while g > 360 do g:=g - 360; // values between 0° - 360°
  while L > 360 do L:=L - 360;
  while ecliptic<0 do ecliptic := ecliptic + 360;

  lambda := L + ( 2*SUNS_ORBIT_EXCENTRICITY*sin(g * (DEGREE_TO_RAD_FACTOR)) + 1.25 *SUNS_ORBIT_EXCENTRICITY*SUNS_ORBIT_EXCENTRICITY
                *sin(2*g * DEGREE_TO_RAD_FACTOR) ) * (RAD_TO_DEGREE_FACTOR);
  // ecliptical length Lambda
  tau := tan((ecliptic/2) * (DEGREE_TO_RAD_FACTOR))*tan((ecliptic/2) * (DEGREE_TO_RAD_FACTOR));
  REC := lambda - (tau * sin(2*lambda * (DEGREE_TO_RAD_FACTOR)) - (Tau*Tau/2)
                * sin(4*lambda*(DEGREE_TO_RAD_FACTOR)) ) * (RAD_TO_DEGREE_FACTOR);
  // Right ascension

  Result:=(L - REC) ;
  Result:=Result*4;   // 4 min per degree

  // accuracy amounts to about 0.1 minutes for the period
  // between 1950 and 2050.

  // The equation is in minutes, we need it in
  // same format as the Julian Day:
  Result:=Result/1440;
end;

// Calculate the day of the year, including the partially fractured day
function DayOfYear(aDate: tDateTime):double; STDCALL;
VAR
    Jahr, Monat, Tag : WORD;
begin
  DecodeDate(aDate, Jahr, Monat, Tag);
  result := aDate -  EncodeDate(Jahr, 1, 1);
end; // Day Of Year


function CoordinatesOfTheSun(julianDate,longitude,latitude : double):TSunsCoordinates; STDCALL;
var
  fullDaysSinceEpocheStart,fullJulianCenturies,siderialTimeInHours,siderialTimeInUT,totalHoursOfCurrentDay,localSiderialTime,numberOfDaysFromJ2000,julianCenturies,meanLongitude,meanAnomaly,equationOfCenter,obliquity,eclipticalLongitude,rightAscension,Lquadrant,RECquadrant,x,y : double;
  result_sun_coordiante : TSunsCoordinates;
const
    daysPerTropicalYear = 365.2422;
    julianDateAtStartOfEpoche_J2000 = 2451545;
    daysPerJulianCentury = 36525;
    hoursPerDay = 24;

begin
  // Number of days from 01.01.2000 12h UT
  fullDaysSinceEpocheStart := int(julianDate - julianDateAtStartOfEpoche_J2000);

  // Number of Julian centuries from J2000
  fullJulianCenturies := fullDaysSinceEpocheStart / daysPerJulianCentury;

  // The input date/time in the JUL function already must be corrected to
  // Universal Time
  totalHoursOfCurrentDay:=((julianDate - julianDateAtStartOfEpoche_J2000) - fullDaysSinceEpocheStart)*hoursPerDay;

  // sidereal time for meridian of Greenwich 00:00h UT in hours
  siderialTimeInHours := 6.697376 + 2400.05134 * fullJulianCenturies ;

  // sidereal time SG of the Greenwich meridian for the universal time tUT
  siderialTimeInUT := siderialTimeInHours + (366.2422 / 365.2422) * totalHoursOfCurrentDay;

  //Local sidereal time S for the geographical longitude LON
  localSiderialTime := siderialTimeInUT*15 + longitude;

  //Number of days d from J2000.0 for the given date and time
  numberOfDaysFromJ2000 := fullDaysSinceEpocheStart + totalHoursOfCurrentDay / hoursPerDay;

  //Number of centuries T from the reference time
  julianCenturies := numberOfDaysFromJ2000 / daysPerJulianCentury;

  // Sun's mean longitude L0 and its mean anomaly M0
  meanLongitude := 280.466 + 36000.770 * julianCenturies;
  meanAnomaly := 357.529 + 35999.050 * julianCenturies;

  // Sun's equation of center C, which accounts for the eccentrisity
  // of Earth's orbit around the Sun
  equationOfCenter := (1.915 - 0.005 * julianCenturies) * sin(meanAnomaly*DEGREE_TO_RAD_FACTOR) + 0.020 * sin(2*meanAnomaly*DEGREE_TO_RAD_FACTOR);

  // True ecliptical longitude of the Sun LS
  eclipticalLongitude := meanLongitude + equationOfCenter;

  // K is the obliquity of the ecliptic
  obliquity := 23.439 - 0.013 * julianCenturies;

  // Sun's right ascension REC
  rightAscension := arctan(tan(eclipticalLongitude*DEGREE_TO_RAD_FACTOR) * cos(obliquity*DEGREE_TO_RAD_FACTOR))*RAD_TO_DEGREE_FACTOR;

  // (right ascension value needs to be in the same quadrant as L)
  Lquadrant  := (floor( eclipticalLongitude/90)) * 90;
  RECquadrant := (floor(rightAscension/90)) * 90;
  rightAscension := rightAscension + (Lquadrant - RECquadrant);

  // Declination DEC
  result_sun_coordiante.declination_angle := arcsin( sin(rightAscension*DEGREE_TO_RAD_FACTOR)*sin(obliquity*DEGREE_TO_RAD_FACTOR) )*RAD_TO_DEGREE_FACTOR;

  // Hour angle of the object HA
  // HRA is 0 if sun at noon, positive in the afternoon
  result_sun_coordiante.hour_angle := (localSiderialTime - rightAscension) + 180;
  // HRA_value := (S - REC);

  // HRA between -180 - 180°
  if result_sun_coordiante.hour_angle>180 then result_sun_coordiante.hour_angle:=result_sun_coordiante.hour_angle-360;

  // Altitude ALT of the Sun

  result_sun_coordiante.altitude_angle := arcsin(sin(latitude*DEGREE_TO_RAD_FACTOR)*sin(result_sun_coordiante.declination_angle*DEGREE_TO_RAD_FACTOR)
    +cos(latitude*DEGREE_TO_RAD_FACTOR)*cos(result_sun_coordiante.declination_angle*DEGREE_TO_RAD_FACTOR)*cos(result_sun_coordiante.hour_angle*DEGREE_TO_RAD_FACTOR))*RAD_TO_DEGREE_FACTOR;

  // Azimuth AZI measured eastward from the North

  y := -cos(result_sun_coordiante.declination_angle*DEGREE_TO_RAD_FACTOR)*cos(latitude*DEGREE_TO_RAD_FACTOR)*sin(result_sun_coordiante.hour_angle*DEGREE_TO_RAD_FACTOR);
  x := sin(result_sun_coordiante.declination_angle*DEGREE_TO_RAD_FACTOR) - sin(latitude*DEGREE_TO_RAD_FACTOR)*sin(result_sun_coordiante.altitude_angle*DEGREE_TO_RAD_FACTOR);

  result_sun_coordiante.azimuthal_angle := arctan (y/x) * RAD_TO_DEGREE_FACTOR;

  if x<0 then
    result_sun_coordiante.azimuthal_angle := result_sun_coordiante.azimuthal_angle + 180;
  if (x>0) and (y<0) then
    result_sun_coordiante.azimuthal_angle := result_sun_coordiante.azimuthal_angle + 360;

  result := result_sun_coordiante;
end;


// Calculates the Universal time from the given local time value
// 13.06.2014: UTC offset input changed from hours to minutes
function UniversalTime(TimeZone : integer; longitude : Double; LocalMeanTime : TDateTime; isZoneTime:boolean):TDateTime; STDCALL;
const
   MINUTES_PER_DAY = 24.0 * 60.0;
begin
  // 13.06.2014 changes in offset input
  if isZoneTime then
    Result := LocalMeanTime - (TimeZone) / MINUTES_PER_DAY
  else
  // Anhand longitude
    Result := LocalMeanTime - ((15*TimeZone/60) - longitude)/(15.0*24.0);
end;

function FaimanMills(DEC_value,HRA_value,LAT_value,ORI_value,TIL_value,SKW_value :double):TPanelCoordinates; STDCALL;
var
  n_sx,n_sy,n_sz,
  n1_sx,n1_sy,n1_sz,
  n2_sx,n2_sy,n2_sz,
  n3_sx,n3_sy,n3_sz,
  n4_sx,n4_sy,n4_sz : double;
  YAN_value, ZAN_value, TTA_value : double;
begin
  n_sx:=cos(DEC_value*DEGREE_TO_RAD_FACTOR)*cos(HRA_value*DEGREE_TO_RAD_FACTOR);
  n_sy:=-cos(DEC_value*DEGREE_TO_RAD_FACTOR)*sin(HRA_value*DEGREE_TO_RAD_FACTOR);
  n_sz:=sin(DEC_value*DEGREE_TO_RAD_FACTOR);

  n1_sx:=n_sz*sin(LAT_value*DEGREE_TO_RAD_FACTOR) + n_sx*cos(LAT_value*DEGREE_TO_RAD_FACTOR);
  n1_sy:=n_sy;
  n1_sz:=n_sz*cos(LAT_value*DEGREE_TO_RAD_FACTOR) - n_sx*sin(LAT_value*DEGREE_TO_RAD_FACTOR);

  n2_sx:=n1_sx;
  n2_sy:=n1_sy*cos(ORI_value*DEGREE_TO_RAD_FACTOR) + n1_sz*sin(ORI_value*DEGREE_TO_RAD_FACTOR);
  n2_sz:=-n1_sy*sin(ORI_value*DEGREE_TO_RAD_FACTOR) + n1_sz*cos(ORI_value*DEGREE_TO_RAD_FACTOR);

  n3_sx:=-n2_sz*sin(TIL_value*DEGREE_TO_RAD_FACTOR) + n2_sx*cos(TIL_value*DEGREE_TO_RAD_FACTOR);
  n3_sy:=n2_sy;
  n3_sz:=n2_sz*cos(TIL_value*DEGREE_TO_RAD_FACTOR) + n2_sx*sin(TIL_value*DEGREE_TO_RAD_FACTOR);

  n4_sx:=n3_sx;
  n4_sy:=n3_sy*cos(SKW_value*DEGREE_TO_RAD_FACTOR) + n3_sz*sin(SKW_value*DEGREE_TO_RAD_FACTOR);
  n4_sz:=-n3_sy*sin(SKW_value*DEGREE_TO_RAD_FACTOR) + n3_sz*cos(SKW_value*DEGREE_TO_RAD_FACTOR);

  if n4_sy<>0 then
    YAN_value := (ARCCOS(ABS(n4_sx)/sqrt(n4_sx*n4_sx+n4_sy*n4_sy))
      *(n4_sy/ABS(n4_sy)))*RAD_TO_DEGREE_FACTOR
  else YAN_value := 0;

  if n4_sz<>0 then
    ZAN_value := (ARCCOS(ABS(n4_sx)/sqrt(n4_sx*n4_sx+n4_sz*n4_sz))
    *(n4_sz/ABS(n4_sz)))*RAD_TO_DEGREE_FACTOR
  else ZAN_value := 0;

  TTA_value := ARCCOS(n4_sx)*RAD_TO_DEGREE_FACTOR;

  result.YAN_value := YAN_value;
  result.ZAN_value := ZAN_value;
  result.TTA_value := TTA_value;
  result.n4_sx_value := n4_sx;
end;

function HorizonToEquatorial(AZI_value,ALT_value,LAT_value : double):TEquatorialCoordinates; STDCALL;
var
  x,y : double;
  DEC_value,HRA_value : double;
begin

  DEC_value := RAD_TO_DEGREE_FACTOR*arcsin(sin(ALT_value*DEGREE_TO_RAD_FACTOR)*sin(LAT_value*DEGREE_TO_RAD_FACTOR)+cos(ALT_value*DEGREE_TO_RAD_FACTOR)*cos(LAT_value*DEGREE_TO_RAD_FACTOR)*cos(AZI_value*DEGREE_TO_RAD_FACTOR));

  y := -cos(ALT_value*DEGREE_TO_RAD_FACTOR)*cos(LAT_value*DEGREE_TO_RAD_FACTOR)*sin(AZI_value*DEGREE_TO_RAD_FACTOR);
  x := sin(ALT_value*DEGREE_TO_RAD_FACTOR) - sin(LAT_value*DEGREE_TO_RAD_FACTOR)*sin(DEC_value*DEGREE_TO_RAD_FACTOR);

  HRA_value := arctan (y/x) * RAD_TO_DEGREE_FACTOR;

  if x<0 then
    HRA_value := HRA_value - 180;
  if (x<0) and (y>0) then
    HRA_value := HRA_value + 360;

  result.DEC_value := DEC_value;
  result.HRA_value := HRA_value;
end;


// Gueymard 1985
// "Une paramétrisation de la luminance énergétique du ciel clair en fonction
// de la turbidité"
function distributionFactorClearSkyRcs(beta, sunsAltitudeInDegree, AOI_value, zenithAngleOfSkyElementInDegree, AOI_lower_limit :double):double; STDCALL;
var
  A, B, F_N, B_C1, N_C1, t, D_c1,airMass,mtr : double;
  A_0, A_1, A_2, C_C1 : double;
  //

  D_C2, b0, b1, b2, b3, b4 : double;

  D_h,a0,a1,a2,a3,Z,D0,Kn,D_C : double;
  sunsZenithAngleInDegree : double;
begin
  sunsZenithAngleInDegree := 90-sunsAltitudeInDegree;

  
  
  // update 25/02/2019: Excluding the solar disk approx. 0.53° (= 2 x 0.26)!
  // limit to 5°
  if AOI_value < 0.5 then AOI_value := 0.5;
  

  AOI_lower_limit := min(5, max(0.26, AOI_lower_limit));
  if (AOI_value>AOI_lower_limit) and (AOI_value<=3) and (sunsAltitudeInDegree>0.01) then
  begin
    // Zone C1 Circumsolar (AOI <= 3°)
    // airMass := abs(1 / sin(sunsAltitudeInDegree*DEGREE_TO_RAD_FACTOR));      // Air mass
    airMass := abs(1 / (cos(sunsZenithAngleInDegree*DEGREE_TO_RAD_FACTOR) + 0.48353 * (Power(sunsZenithAngleInDegree,0.095846))/Power(96.741-sunsZenithAngleInDegree,1.754)));


    if (airMass*beta)<1 then
      t := 6.556 * sqrt(airMass*beta) - 3.346*(airMass*beta)
    else // (((m*beta)>=1) and ((m*beta)<3))
      t := 3.210 + 0.200 * sqrt((airMass*beta)-1);

    N_C1 := Power(AOI_value,-t);

    mtr := 0.11008 + (0.07260 + 0.04077*ln(airMass))*ln(airMass);

    if airMass <= 2.4 then
    begin
      A := 8.03  + (0.931 - 0.179*airMass)*airMass;
      B := 0.243 - 0.023*airMass;
    end else
    begin
      A := 8.093 + Power(airMass,0.15);
      B := 0.14 + exp(-2.211 - 0.342*airMass);
    end;

    F_N := exp(A * Power(beta,B));

    // Update 25/02/2019: Interpreation of ambiguous equation 3/8pi = 3/(8*pi)
    B_C1 := 3/(8 * Pi) * mtr * F_N;

    try
      A_0 := 3.17 - 1.98*airMass + 0.059*airMass*airMass + 11.36*sqrt(airMass);
      A_1 := 3.93 - 3.29*airMass + 0.061*airMass*airMass - 8.62*sqrt(airMass);
      A_2 := exp(2.780 + 0.121*airMass - 0.233*sqrt(airMass));

      C_C1 := (A_0 + A_1*beta)/(1 + A_2*beta);

      D_C1 := N_C1 * B_C1 * C_C1;

    except
      D_C1 := 0;
    end;

  end else
    D_C1 := 0;

 // D_C1 := 0;

  // Zone C2 (3° < AOI <= 20°)

  if (AOI_value>AOI_lower_limit) and (AOI_value > 3) and (AOI_value <= 20) then
  begin
    b0 := 0.109 + 0.029 * AOI_value + 0.005 * exp( -0.015 * sunsZenithAngleInDegree + 1.07E-5 * sunsZenithAngleInDegree * sunsZenithAngleInDegree * sunsZenithAngleInDegree );
    b1 := 0.02 - 6E-4 * AOI_value;
    b2 := 6.8E-5;
    b3 := 0.24;
    b4 := -0.054;

    D_C2 := (b0 + b1 * sunsZenithAngleInDegree + b2 * sunsZenithAngleInDegree * sunsZenithAngleInDegree) * exp(b3 + b4 * AOI_value);
  end else
    D_C2 := 0;
  // hemispherical

  D_h := (1 + 0.01*sunsZenithAngleInDegree) * (0.275 - 0.395 * cos(zenithAngleOfSkyElementInDegree*DEGREE_TO_RAD_FACTOR) + 0.170 * cos(zenithAngleOfSkyElementInDegree*DEGREE_TO_RAD_FACTOR)*cos(zenithAngleOfSkyElementInDegree*DEGREE_TO_RAD_FACTOR));

  // final model
  Z := 0.01 * sunsZenithAngleInDegree;
  a0 := 1.0156 + 0.0907*beta - 0.8644*beta*beta;
  a1 := -0.1966 + 1.5843*beta - 3.8185*beta*beta;
  a2 := 0.3651 - 4.8270*beta + 7.9650*beta*beta;
  a3 := -0.113 + 2.010*beta - 2.950*beta*beta;

  KN := a0 + a1 * Z + a2 * Z*Z + a3 * Z*Z*Z;

  D_c := max(D_c1, D_c2);

  D0 := KN* ( D_c + D_h );

  // function output

  result := D0;
end;
// [17] S. Fritz
// Illuminance and luminance under overcast skies
function distributionFactorOvercastSkyRov(zenithAngleOfSkyElementInDegree,cloudOpacity :double):double; STDCALL;
var
  b : double;
begin
  b := 0.5 + cloudOpacity;

  // Update 25/02/2019: Original equation for D1
  result := ( 1 + b*cos(zenithAngleOfSkyElementInDegree*DEGREE_TO_RAD_FACTOR))/ ( 1 + b);

  // Equation without denominator (1 + b) -> increasing calculation speed, due to normalisation
  // result := ( 1 + b*cos(zenithAngleOfSkyElementInDegree*DEGREE_TO_RAD_FACTOR));
end;

// angle between sun (AZI,ALT) and sky element (theta, zeta)
function AngleBetweenObjects(AZI_value,ALT_value,theta,zeta:double):double; STDCALL;
begin
  result := (RAD_TO_DEGREE_FACTOR)*arccos(sin(ALT_value*DEGREE_TO_RAD_FACTOR)*sin(zeta*DEGREE_TO_RAD_FACTOR)+cos(ALT_value*DEGREE_TO_RAD_FACTOR)
    *cos(zeta*DEGREE_TO_RAD_FACTOR)*cos((theta-AZI_value)*DEGREE_TO_RAD_FACTOR));
end;


function SlopeFactor(distributionFactorClearSky,distributionFactorOvercastSky,cloudOpacity : double):double; STDCALL;
begin
  result := (1 - cloudOpacity)*distributionFactorClearSky + cloudOpacity*distributionFactorOvercastSky;
end;


// This function calculates the cloud opacity Npt from measured weather data
// Update 15.12.2014 // Proposed by Matthias Winter
// Ib - beam irradiance
// Id - diffuse irradiance
// Ib_0 - Clear sky beam irradiance
// Id_0 - Clear sky diffuse irradiance
function CloudOpacityWinter(Ib, Id, Ib_0, Id_0 : double):double; STDCALL;
var
  rho, rho_0 : double;
begin
  try
    rho := Id / (Ib + Id);
    rho_0 := Id_0 / (Ib_0 + Id_0);
    result := (rho - rho_0)/(1 - rho_0);
  except
    result := 0;
  end;

  if result<0 then result := 0;
  if result>1 then result := 1;
end;

function SunAboveHorizont(ALT_value:double):boolean; STDCALL;
begin
  if (ALT_value>=0.1) then result:=true else result:=false
end;

function SunOnPanel(TTA_condition:boolean;n4_sx_value:double):boolean; STDCALL;
begin
 if TTA_condition=false then
  if (n4_sx_value>=0) and (n4_sx_value<=1) then result:=true else result:=false
 else result:=true;
end;

function Bifacial_Panels(TTA_condition:boolean;bifacialEfficiency:double;
  n4_sx_value:double):real; STDCALL;
begin
 if TTA_condition=false then result:=1 else
  if (n4_sx_value>=0) and (n4_sx_value<=1) then result:=1 else
    result:=bifacialEfficiency;
end;

function BoolToInt(a : boolean):integer; STDCALL;
begin
  if a then result:=1 else result:=0;
end;

function calculateCosinelawCorrectionFactor(isInputHorizontal, wantOutputNormal : boolean; incidenceAngleHorizontal, incidenceAnglePanel : double):double; STDCALL;
begin
  if isInputHorizontal then
    if wantOutputNormal then
      result := 1 / abs(cos(incidenceAngleHorizontal*DEGREE_TO_RAD_FACTOR))
    else
      if incidenceAnglePanel=incidenceAngleHorizontal then
        result := 1
      else
        result := abs(cos(incidenceAnglePanel*DEGREE_TO_RAD_FACTOR)) / abs(cos(incidenceAngleHorizontal*DEGREE_TO_RAD_FACTOR))
  else
    if wantOutputNormal then
      result := 1
    else
      result := abs(cos(incidenceAnglePanel*DEGREE_TO_RAD_FACTOR));
end;

EXPORTS
JulianDay,
EquationOfTime,
DayOfYear,
CoordinatesOfTheSun,
UniversalTime,
FaimanMills,
HorizonToEquatorial,
CloudOpacityWinter,
SlopeFactor,
AngleBetweenObjects,
distributionFactorClearSkyRcs,
distributionFactorOvercastSkyRov,
SunAboveHorizont,
SunOnPanel,
Bifacial_Panels,
BoolToInt,
calculateCosinelawCorrectionFactor;


begin
end.