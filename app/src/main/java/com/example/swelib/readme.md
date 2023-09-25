* This is a simple example demonstrating the use of the Swiss Ephemeris library in a real Android application.
* No library needed, just JNI.
* You can find more information and details about Swiss Ephemeris at https://www.astro.com.
* The code is based on krymlov's JNI implementation, available at https://github.com/krymlov, and includes a class wrapped to simplify usage.
* Using JNI directly offers greater flexibility, making it easy to add or remove functionality in C code.
* This demo attempts to mimic the functionality of pyswisseph (the Swiss Ephemeris library in Python) as closely as possible, making it user-friendly.
* Please note that this demo does not provide 100% API coverage; the PlaCalc function, which is now a legacy feature and not used in Swiss Ephemeris, is not included.



# Time
### Time input
1. Before 1972-1-1 (or in the system of java.time, before 1972-11-3/4), time is **UT1**
2. After that instant above, time is **UTC**

### Time's Name
* ET: Ephemeris time;
* TT: Terrestrial Time;
* ET = TT, used in calculation related to the solar system.

* UT1: Universal Time, mean solar time, used in calculation of ASC, MC, houses;
* UTC: COORDINATED UNIVERSAL TIME, the time read from a clock;
* TAI: International Atomic Time.
### Time Transform
1. DeltaT = **ET** - **UT1**, (use _calc_deltat_ or _swe_deltat_ex_ in swephlib.c; changes with time)
2. NLEAP_INIT = 10/11 (initial leap second is 10 at 1972-1-1, or 11 at 1972-11-4 (java.time))
3. nleap = **TAI** - **UTC** (obtained from the list of leap seconds)
4. DUT1 = **UT1** - **UTC** (measured and broadcast weekly by NIST, with an absolute value < 0.9 sec, adjusted via leap seconds)
5. ET - TAI = 32.184 sec (a const val)

|             |          ET           |               UT               | UTC  | 
|:-----------:|:---------------------:|:------------------------------:|:----:|  
| before 1972 |      tdj+DeltaT       |              tdj               |      |  
| after 1972  | tjd + 32.184s + nleap | tjd + 31.184s + nleap - DeltaT | tdj  |  

(tdj here is julian day number, get from SwissEphKt.utc_to_jdKT)
(pay attention to the unit, second/day)

### Other Time
* Local Mean Time/ Mean Solar Time: UT1, or UT1 + ZoneOffset(use java.time)
* Local Apparent Time/ True Solar Time: Local Mean Time + a parameter (use SwissEphKt.swe_lmt_to_lat)

### Other Notes
* TT sounds more modern, but ET was used in Swiss Ephemeris.
* in 寿星万年历, TT/ET is referred as 力学时
* in swisseph, Gregorian is actually Proleptic Gregorian

# Calendar
|Proleptic Gregorian Range (java.time)|Julian|
|---|---|---|
| 4-2-27 to 100-2-27| -2 |
| 100-2-28 to 200-2-28| -1 |
| 200-3-1 to 300-3-1| 0 |
| 300-3-2 to 500-3-2| 1 |
| 500-3-3 to 600-3-3| 2 |
| 600-3-4 to 700-3-4| 3 |
| 700-3-5 to 900-3-5| 4 |
| 900-3-6 to 1000-3-6| 5 |
| 1000-3-7 to 1100-3-7| 6 |
| 1100-3-8 to 1300-3-8| 7 |
| 1300-3-9 to 1400-3-9| 8 |
| 1400-3-10 to 1500-3-10| 9 |
| 1500-3-11 to 1582-10-14| 10 |
| 1582-10-15 to 1700-2-28|10|
Convert: use fun convertISO86ToGre and fun convertGreToISO86