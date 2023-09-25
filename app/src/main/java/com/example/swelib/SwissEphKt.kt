package com.example.swelib

import android.util.Log
import java.io.Closeable

import java.time.LocalDate
import java.time.ZoneOffset
import java.time.ZonedDateTime
import java.time.temporal.ChronoUnit
import kotlin.math.floor

private const val TAG = "serr"
private const val ERR = -1

/**
 * For a complete explanation, refer to pyswisseph and swiss ephemeris' docs.
 *
 * I try to make everything following pyswisseph.
 * Except:
 * - **PlaCalc functions**: PLACALC, the predecessor of SWISSEPH, had included several functions that we do not need for SWISSEPH anymore.
 * @see <a href="https://www.astro.com/swisseph/swephprg.htm">Swiss Ephemeris</a>
 * @see <a href="https://astrorigin.com/pyswisseph/pydoc/index.html">pyswisseph</a>
 */
object SwissEphKt: Closeable {
    internal val swe = SwissEph()
    private val julianDayUT = ZonedDateTime.of(-4713, 11, 24, 12, 0, 0, 0, ZoneOffset.UTC)
    /**
     * from utc to jd_utc; ut to jd_ut...
     */
    fun julday(zonedDateTime: ZonedDateTime): Double {
        return ChronoUnit.MILLIS.between(julianDayUT, zonedDateTime) / 86400000.0
    }
    /**
     * Calculate a Julian day number.
     * @param hour: the time of day, decimal with fraction
     * @param cal: either GREG_CAL (gregorian) or JUL_CAL (julian)
     * @return jd
     */
    fun julday(year:Int, month:Int = 1, day:Int = 1,
               hour:Double = 0.0, cal:Int = 1): Double {
        return swe.julday(year, month, day,
            hour, cal)
    }

    /**
     * from jd_utc to utc; jd_utc to utc...
     */
    fun revjul(tjut: Double): ZonedDateTime? {
        return julianDayUT.plus((tjut * 86400000).toLong(), ChronoUnit.MILLIS )
    }

    /**
     * Calculate year, month, day, hour from Julian day number.
     * @param jd: Julian day number
     * @param cal: either GREG_CAL (gregorian) or JUL_CAL (julian)
     * @return year: Int, month: Int, day: Int, hour: Double
     */
    fun revjul(jd: Double, cal:Int = 1):List<Any> {
        val ymd = IntArray(3)
        val hour = DoubleArray(1)
        swe.revjul(jd, cal, ymd,hour)
        return listOf(ymd[0], ymd[1], ymd[2], hour[0])
    }
    /**
     * from java.time to Pair(ET, UT1) in julian day number
     */
    fun utc_to_jd(zonedDateTime: ZonedDateTime): Pair<Double, Double> {
        val leap_seconds = listOf(
            19720630,
            19721231,
            19731231,
            19741231,
            19751231,
            19761231,
            19771231,
            19781231,
            19791231,
            19810630,
            19820630,
            19830630,
            19850630,
            19871231,
            19891231,
            19901231,
            19920630,
            19930630,
            19940630,
            19951231,
            19970630,
            19981231,
            20051231,
            20081231,
            20120630,
            20150630,
            20161231,
            0  /* keep this 0 as end mark */
        )
        //TAI-UTC in sec
        val NLEAP_INIT = 10
        //tt-TAI in sec
        val ttTAI = 32.184
        val year = zonedDateTime.year
        val serr = StringBuilder()
        val tjd = julday(zonedDateTime)
        //UT1-UTC in day
        val d = swe.deltat(tjd)
        if (serr.isNotEmpty()) {
            Log.d(TAG, "utc_to_jd: $serr")
        }
        val ndat = year * 10000 + zonedDateTime.monthValue * 100 + zonedDateTime.dayOfMonth
        //Log.d(TAG, "utc_to_jdKT: tdj  $tdj,  nleap ,  d:  $d")
        // according to the instant rule of java.time, before 1972-11-3 was ut1
        if (ndat <= 19721103 || year >= 2035) {
            return Pair(tjd + d, tjd)
        }

        var nleap = NLEAP_INIT  /* initial difference between UTC and TAI in 1972 */
        for (i in leap_seconds) {
            if (ndat <= i)
                break
            nleap++
        }

        /*
         * For input dates > today:
         * If leap seconds table is not up to date, we'd better interpret the
         * input time as UT1, not as UTC. How do we find out?
         * Check, if delta_t - nleap - 32.184 > 0.9
         */
        if (year >= 2023 && d * 86400.0 - nleap - ttTAI >= 0.9) {
            return Pair(tjd + d, tjd)
        }
        val tjd_et = (ttTAI + nleap) / 86400.0 + tjd
        //Well, a better way is to recalculate the d with
        val d1 = swe.deltat(tjd_et - d)
        if (d1!=d){
            val d2 = swe.deltat(tjd_et - d1)
            return Pair(tjd_et, tjd_et - d2)
        }else{
            return Pair(tjd_et, tjd_et -d)
        }



        //return Pair(tjd_et, tjd_et - d)
    }
    /**
     * Convert UTC to julian day.
     * @param year: Year
     * @param cal: either GREG_CAL or JUL_CAL
     * @return jdet: Julian day in ET (TT); jdut: Julian day in UT (UT1)
     */
    fun utc_to_jd(year: Int, month: Int =1,day: Int=1,hour:Int=0,minutes:Int=0,seconds: Double=0.0,
                  cal: Int =1): Pair<Double, Double>{
        val dret = DoubleArray(2)
        val serr = StringBuilder()
        val valueError = swe.utc_to_jd(year, month , day, hour,minutes, seconds,cal,dret, serr)
        if (valueError == -1){
            Log.d(TAG, "utc_to_jd: the cal is not GREG_CAL or JUL_CAL")
        }else if (serr.isNotEmpty()){
            Log.d(TAG, "utc_to_jd: $serr")
        }
        return Pair(dret[0],dret[1])
    }

    /**
     * Calculate planetary position with UT
     * return Pair(result[6], and return flag)
     */
    fun cal_ut(tdj_ut:Double,planet:Int, iflag: Int):Pair<List<Double>,Int>{
        val serr = StringBuilder()
        val xx = DoubleArray(6)
        val retflag = swe.calc_ut(tdj_ut, planet, iflag, xx, serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "cal_ut: $serr")
        }
        if (retflag<0){
            Log.d(TAG, "cal_ut: fatal error code $retflag")
        }
        return Pair(xx.toList(), retflag)
    }
    /**
     * convert ISO 8601 to Gregorian+Julian mix, non-existence date does not exist
     * @param localDate Date in ISO8601
     * @return before 1582-10-14 convert to Julian, else to Gregorian. Pair(year: Int, month: Int, day: Int) (Jan start at 1)
     * basically the same as: GregorianCalendar.from(zd).time, but 20* faster than GregorianCalendar -> reduce the useless timezone,
     * taking care of the timezone is the reason of most Calendar chaos
     */
    fun convertISO86ToGre(localDate: LocalDate):Triple<Int, Int, Int>{
        if (localDate.isAfter(LocalDate.of(1582,10,14))){
            return Triple(localDate.year,localDate.monthValue, localDate.dayOfMonth)
        }else{
            val epochDayRevised = localDate.toEpochDay() + 171595
            val wholeCycle = Math.floorDiv(epochDayRevised, 146100)
            return when (epochDayRevised - 146100 * wholeCycle){
                146099L-> Triple(localDate.year, 2,29)
                in 109575 .. 146098-> {val revisedDate = localDate.minusDays((12+3*wholeCycle))
                    Triple(revisedDate.year, revisedDate.monthValue, revisedDate.dayOfMonth)
                }
                109574L -> Triple(localDate.year, 2,29)
                in 73050 ..109573 -> {val revisedDate = localDate.minusDays((11+3*wholeCycle))
                    Triple(revisedDate.year, revisedDate.monthValue, revisedDate.dayOfMonth)
                }
                73049L -> Triple(localDate.year, 2,29)
                in 0 .. 73048-> {val revisedDate = localDate.minusDays((10+3*wholeCycle))
                    Triple(revisedDate.year, revisedDate.monthValue, revisedDate.dayOfMonth)
                }
                else -> Triple(Int.MAX_VALUE,Int.MAX_VALUE,Int.MAX_VALUE)
            }
        }
    }
    /**
     * convert Gregorian+Julian mix to ISO 8601, non-existence date will move forward
     * @param month make sure it is in 1 .. 12
     * @param day make sure it is in 1 .. 31
     * @return LocalDate in ISO8601
     * basically the same as GregorianCalendar(y,m-1,d).toZonedDateTime(), but 4 times faster,
     * except GregorianCalendar has different rules to deal with non-existence date.
     */
    fun convertGreToISO86(year: Int, month:Int, day:Int): LocalDate{
        fun toLocalDate(y: Int, m: Int, d: Int): LocalDate {
            try {
                return LocalDate.of(y, m, d)
            } catch(e: Exception){
                return LocalDate.of(y, m+1, 1)
            }
        }

        val timeNumber = year*10000+month*100+day
        if (timeNumber>=15821015){
            return toLocalDate(year,month,day)
        } else if(timeNumber in 15821005..15821014){
            return LocalDate.of(1582,10,15)
        } else {
            val timeNumberModified = timeNumber-229+1_000_000
            val cycle = Math.floorDiv(timeNumber-229+1_000_000, 4_000_000)
            val century =  Math.floorDiv(timeNumberModified,10000) - cycle * 400
            when (century) {
                in 0..99 -> {
                    print(cycle)
                    print("\t")
                    print(century)
                    return if (month ==2 && day == 29){
                        LocalDate.of(year,3,1).plusDays(3* cycle -3L)
                    } else {
                        toLocalDate(year,month,day).plusDays(3* cycle -2L)
                    }
                }
                in 100..199-> {
                    print(cycle)
                    print("\t")
                    print(century)
                    return toLocalDate(year,month,day).plusDays(3* cycle -2L)
                }
                in 200..299-> {
                    print(cycle)
                    print("\t")
                    print(century)
                    return toLocalDate(year,month,day).plusDays(3* cycle -1L)
                }
                in 300..399-> {
                    print(cycle)
                    print("\t")
                    print(century)
                    return toLocalDate(year,month,day).plusDays(3* cycle.toLong())
                }
                else-> {
                    print(cycle)
                    print("\t")
                    print(century)
                    return LocalDate.MAX
                }
            }
        }
    }
    private fun ganzhiSequence(gan: Int, zhi: Int): Int {
        return (6 * gan - 5 * zhi) % 60
    }
    /**
     * Convert Ganzhi to possible Julian Day Number in Universal Time.
     * @param yearGan 甲 is 1, 1-12
     * @param yearZhi 寅 is 1, 1-12
     * @param dayGan 甲 is 1, 1-12
     * @param dayZhi 寅 is 1, 1-12
     */
    fun convertBaziToDatesKT(yearGan: Int, yearZhi: Int, month: Int,
                             dayGan: Int, dayZhi:Int, hour: Int,
                             longitude: Double,
                             rangeMin: Int = 0,
                             rangeMax: Int = 2100,
                             ): List<Double> {
        val yearrange = (rangeMax - rangeMin)/60 + 1
        val result = DoubleArray(yearrange)
        val find = swe.convertBaziToDates(ganzhiSequence(yearGan,yearZhi), month, ganzhiSequence(dayGan,dayZhi), hour,
            longitude, tjd_ut = result)
        return result.take(find).toList()
    }
    /**
     * computes the horizontal coordinates (azimuth and altitude)
     *
     *  The **apparent altitude** of a body depends on the atmospheric pressure and temperature.
     *  If only the true altitude is required, these parameters can be neglected.
     *
     *  If [atpress] is given the value 0, the function estimates the pressure from the geographical altitude given in [xin] and [attemp]. If [xin] is 0, [atpress] will be estimated for sea level.
     * @param tjdut Julian day number, Universal Time
     * @param flag either from ecliptical coord (ECL2HOR) or equatorial (EQU2HOR)
     * @param geopos DoubleArray[3] with:
     *
     *      - 0: geographic longitude, in degrees (eastern positive)
     *      - 1: geographic latitude, in degrees (northern positive)
     *      - 2: geographic altitude, in meters above sea level
     * @param atpress atmospheric pressure in mbar (hPa), if set to 0.0 pressure will estimated by geopos
     * @param attemp atmospheric temperature in degrees Celsius, 15 degree here as mean temp
     * @param xin a DoubleArray(3) with:
     *
     *     - ECL2HOR: ecl. longitude, ecl. latitude, distance (0)
     *     - EQU2HOR: right ascension, declination, distance (1)
     * @return DoubleArray[3] with: azimuth, true_altitude, apparent_altitude:
     *  - azimuth: position degree, measured from south point to west
     *  - true_altitude: true altitude above horizon in degrees
     *  - apparent_altitude: apparent (refracted) altitude above horizon in   degrees
     *
     *
   */
    fun azalt(tjdut: Double,
              flag: Int = 1,
              geopos: DoubleArray,
              atpress: Double = 0.0,
              attemp: Double = 15.0,
              xin: DoubleArray?): DoubleArray {
        return swe.azalt(tjdut,flag,geopos,atpress,attemp,xin)
    }
    /**
     * Calculate either ecliptical or equatorial coordinates from azimuth and true altitude.
     *
     * This function is not precisely the reverse of [azalt]. It computes either ecliptical or equatorial coordinates from azimuth and true altitude. If only an apparent altitude is given, the true altitude has to be computed first with the function [refrac].
     * @param tjdut input time, Julian day number, Universal Time
     * @param flag either HOR2ECL (to ecliptical) or HOR2EQU (to equatorial)
     * @param geopos a DoubleArray(3) with:
     *
     *      - 0: geographic longitude, in degrees (eastern positive)
     *      - 1: geographic latitude, in degrees (northern positive)
     *      - 2: geographic altitude, in meters above sea level
     * @param azimuth position degree, measured from south point to west
     * @param attemp atmospheric temperature in degrees Celsius, 15 degree here as mean temp
     * @param true_altitude true altitude above horizon in degrees
     * @return ecliptical or equatorial coordinates, depending on flag
     *
     *
     */
    fun azalt_rev(
        tjdut: Double,
        flag: Int,
        geopos: DoubleArray,
        azimuth: Double,
        true_altitude: Double
    ): List<Double> {
        val xin = doubleArrayOf(azimuth, true_altitude)
        val xout = DoubleArray(2)
        swe.azalt_rev(tjdut,flag, geopos, xin, xout)
        return xout.toList()
    }

    /**
     * Calculate planetary positions (ET).
     * @param tjdet Julian day, Ephemeris Time, where tjdet == tjdut + deltat(tjdut)
     * @param planet body number
     * @param flags bit flags indicating what kind of computation is wanted
     * @return (xx), int retflags:
     * - xx: List of 6 Double for results: longitude, latitude, distance, speed in long., speed in lat., and speed in dist.
     * - retflags: bit flags indicating what kind of computation was done
     */
    fun calc(
        tjdet: Double,
        planet: Int,
        flags: Int
    ): Pair<List<Double>, Int> {
        val xx = DoubleArray(6)
        val serr = StringBuilder()
        val valueError = swe.calc(tjdet, planet, flags, xx, serr)
        if (valueError == ERR){
            Log.d(TAG, "calc: Fatal Exception")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "calc: $serr")
        }
        return Pair(xx.toList(), flags)
    }
    /**
     * Calculate planetocentric positions of planets (ET).
     *
     * @param tjdet Julian day, Ephemeris Time, where tjdet == tjdut + deltat(tjdut)
     * @param planet body number
     * @param center body number of center object
     * @param flags bit flags indicating what kind of computation is wanted
     * @return (xx), int retflags:
     * - xx: list of 6 double for results: longitude, latitude, distance, speed in long., speed in lat., and speed in dist.
     * - retflags: bit flags indicating what kind of computation was done
     */
    fun calc_pctr(
        tjdet: Double,
        planet: Int,
        center: Int,
        flags: Int
    ): Pair<List<Double>, Int> {
        val xx = DoubleArray(6)
        val serr = StringBuilder()
        val valueError = swe.calc_pctr(tjdet, planet, center, flags, xx, serr)
        if (valueError == -1){
            Log.d(TAG, "calc_pctr: Fatal Error")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "calc_pctr: $serr")
        }
        return Pair(xx.toList(), flags)
    }
    /**
     * Calculate planetary positions (UT).
     * @param tjdut Julian day, Ephemeris Time, where tjdet == tjdut + deltat(tjdut)
     * @param planet body number
     * @param flags bit flags indicating what kind of computation is wanted
     * @return (xx), int retflags:
     * - xx: List of 6 Double for results: longitude, latitude, distance, speed in long., speed in lat., and speed in dist.
     * - retflags: bit flags indicating what kind of computation was done
     */
    fun calc_ut(
        tjdut: Double,
        planet: Int,
        flags: Int = 258
    ): Pair<List<Double>, Int> {
        val xx = DoubleArray(6)
        val serr = StringBuilder()
        val valueError = swe.calc_ut(tjdut, planet, flags, xx, serr)
        if (valueError == -1){
            Log.d(TAG, "calc_ut: Fatal Error")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "calc_ut: $serr")
        }
        return Pair(xx.toList(), flags)
    }

    /**
     * Close Swiss Ephemeris.
     *
     * At the end of your computations you can release all resources (open files and allocated memory) used by the swisseph module.
     *
     *
     * After [close], no swisseph functions should be used unless you call [set_ephe_path] again and,
     * if required, [set_jpl_file].
     */
    override fun close(){
        swe.close()
    }

    /**
     * Coordinate transformation from ecliptic to equator or vice-versa.
     *
     * From equatorial to ecliptical, obliquity must be positive.
     *
     * From ecliptical to equatorial, obliquity must be negative.
     *
     * Longitude, latitude and obliquity are in positive degrees.
     *
     * @param coord DoubleArray(3) for coordinates:
     *
     *     - 0: longitude
     *     - 1: latitude
     *     - 2: distance (unchanged, can be set to 1)
     * @param eps obliquity of ecliptic, in degrees
     * @return retlon, retlat, retdist
     *  - retlon: converted longitude
     *  - retlat: converted latitude
     *  - retdist: converted distance
     */
    fun cotrans(coord: DoubleArray, eps: Double): List<Double> {
        val xpn = DoubleArray(3)
        swe.cotrans(coord, xpn, eps, )
        return xpn.toList()
    }

    /**
     * Coordinate transformation of position and speed, from ecliptic to equator or vice-versa.
     *
     * From equatorial to ecliptical, obliquity must be positive.
     *
     * From ecliptical to equatorial, obliquity must be negative.
     *
     * Longitude, latitude, their speeds and obliquity are in positive degrees.
     * @param coord DoubleArray(6) for coordinates:
     *
     *     - 0: longitude
     *     - 1: latitude
     *     - 2: distance
     *     - 3: longitude speed
     *     - 4: latitude speed
     *     - 5: distance speed
     * @param eps obliquity of ecliptic, in degrees
     * @return retlon, retlat, retdist
     *  - retlon: converted longitude
     *  - retlat: converted latitude
     *  - retdist: converted distance
     */
    fun cotrans_sp(coord: DoubleArray, eps: Double): List<Double> {
        val xpn = DoubleArray(6)
        swe.cotrans_sp(coord, xpn, eps, )
        return xpn.toList()
    }
    /**
     * Calculate Julian day number with check wether input date is correct.
     * @param year input year
     * @param month input month
     * @param day input day
     * @param hour input time, decimal with fraction
     * @param cal calendar type, gregorian ('g') or julian ('j')
     * @return isvalid: Boolean, jd: Double
     *  - isvalid: True if the input date and time are legal
     *  - jd: returned Julian day number
     */
    fun date_conversion(
        year:Int,
        month:Int,
        day:Int,
        hour:Double,
        cal:Char
    ): Pair<Boolean, Double> {
        val tjd = DoubleArray(1)
        val isValid = swe.date_conversion(year,month,day,hour,cal,tjd)!=ERR
        return Pair(isValid, tjd.first())
    }
    /**
     * Calculate day of week number [0;6] from Julian day number (monday is 0).
     */
    fun day_of_week(jd:Double): Int {
        return (Math.floorMod(floor(jd+0.5).toInt(),7))
    }
    /**
     * Calculate midpoint (in degrees).
     */
    fun deg_midp(x1:Double,x2:Double): Double {
        return swe.deg_midp(x1,x2)
    }

    /**
     * Normalization of any degree number to the range [0:360]
     */
    fun degnorm(x:Double): Double {
        return swe.degnorm(x)
    }

    /**
     * Calculate value of delta T from Julian day number.
     *
     * tjdet == tjdut + deltat(tjdut)
     *
     * This function is safe only if your application consistently uses the same ephemeris flags,
     * if your application consistently uses the same ephemeris files,
     * if you first call [set_ephe_path] (with flag ``FLG_SWIEPH``) or [set_jpl_file] (with flag ``FLG_JPLEPH``).
     *
     * Also, it is safe if you first call [set_tid_acc] with the tidal acceleration you want.
     * However, do not use that function unless you know what you are doing.
     *
     * For best control of the values returned, use function [deltat_ex] instead.
     *
     * The calculation of ephemerides in UT depends on Delta T,
     * which depends on the ephemeris-inherent value of the tidal acceleration of the Moon.
     * In default mode, the function [deltat] automatically tries to find the required values.
     *
     * Two warnings must be made, though:
     *  - It is not recommended to use a mix of old and new ephemeris files, because the old files were based on JPL Ephemeris DE406, whereas the new ones are based on DE431, and both ephemerides have a different inherent tidal acceleration of the Moon. A mixture of old and new ephemeris files may lead to inconsistent ephemeris output. Using old asteroid files ``se99999.se1`` together with new ones, can be tolerated, though.
     *  - The function [deltat] uses a default value of tidal acceleration (that of DE431).
     *  However, after calling some older ephemeris, like Moshier ephemeris, DE200, or DE406, [deltat] might provide slightly different values.
     *
     * In case of troubles related to these two points, it is recommended to either use function [deltat_ex], or control the value of the tidal acceleration using the functions [set_tid_acc] and [get_tid_acc].
     *
     * @param tjdut input time, Julian day number, Universal Time
     * @return delta T value
     */
    fun deltat(tjdut: Double): Double {
        return swe.deltat(tjdut)
    }

    /**
     * Calculate value of Delta T from Julian day number (extended).
     *
     * Call this function after a previous call of [set_ephe_path] or [set_jpl_file].
     *
     * The calculation of ephemerides in UT depends on the ephemeris-inherent value of the tidal acceleration of the Moon.
     * The function [deltat_ex] can provide ephemeris-dependent values of Delta T and is therefore better than the old function [deltat],
     * which has to make un uncertain guess of what ephemeris is being used. One warning must be made, though:
     *
     * It is not recommended to use a mix of old and new ephemeris files, because the old files were based on JPL Ephemeris DE406,
     * whereas the new ones are based on DE431, and both ephemerides have a different inherent tidal acceleration of the Moon.
     * A mixture of old and new ephemeris files may lead to inconsistent ephemeris output.
     * +Using old asteroid files ``se99999.se1`` together with new ones, can be tolerated, though.
     *
     * @param tjdut input time, Julian day number, Universal Time
     * @param flag ephemeris flag, **FLG_SWIEPH**, **FLG_JPLEPH**, **FLG_MOSEPH**
     * If iflag = -1, then the default tidal acceleration is ussed (i.e. that of DE431).
     * @return delta T value
     */
    fun deltat_ex(tjdut:Double, flag:Int = -1): Double {
        val serr = StringBuilder()
        val deltat =  swe.deltat_ex(tjdut,flag,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "deltat_ex: $serr")
        }
        return deltat
    }

    /**
     * Calculate fixed star positions (ET).
     * @param star name of fixed star to search for
     * @param tjdet input time, Julian day number,  Ephemeris Time
     * @param flags bit flags indicating what kind of computation is wanted
     * @return xx stnam, int retflags
     *  - xx: 6 doubles for results
     *  - stnam: returned star name
     *  - retflags: bit flags indicating what kind of computation was done
     */
    fun fixstar(
        star: String,
        tjd_et: Double,
        iflag: Int = 258,
    ): Triple<List<Double>, String, Int> {
        val starName = StringBuilder(star)
        val serr = StringBuilder()
        val xx = DoubleArray(6)
        val returncode = swe.fixstar(starName, tjd_et, iflag, xx,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "fixstar: $serr")
        }
        if (returncode == ERR){
            Log.d(TAG, "fixstar: Fatal Error")
        }
        return Triple(xx.toList(), starName.toString(),returncode)
    }
    /**
     * Calculate fixed star positions (faster version if many stars are calculated) (ET).
     * @param star name of fixed star to search for
     * @param tjdet input time, Julian day number,  Ephemeris Time
     * @param flags bit flags indicating what kind of computation is wanted
     * @return xx stnam, int retflags
     *  - xx: 6 doubles for results
     *  - stnam: returned star name
     *  - retflags: bit flags indicating what kind of computation was done
     */
    fun fixstar2(
        star: String,
        tjd_et: Double,
        iflag: Int = 258,
    ): Triple<List<Double>, String, Int> {
        val starName = StringBuilder(star)
        val serr = StringBuilder()
        val xx = DoubleArray(6)
        val returncode = swe.fixstar2(starName, tjd_et, iflag, xx,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "fixstar2: $serr")
        }
        if (returncode == ERR){
            Log.d(TAG, "fixstar2: Fatal Error")
        }
        return Triple(xx.toList(), starName.toString(),returncode)
    }

    /**
     * Get fixed star magnitude (faster version if many stars are calculated).
     * @param star name of fixed star
     * @return mag, stnam
     *  - mag: returned magnitude
     *  - stnam: returned star name
     */
    fun fixstar2_mag(
        star: String,
    ): Pair<Double, String> {
        val starName = StringBuilder(star)
        val serr = StringBuilder()
        val mag = DoubleArray(1)
        val returncode = swe.fixstar2_mag(starName, mag,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "fixstar2_mag: $serr")
        }
        if (returncode == ERR){
            Log.d(TAG, "fixstar2_mag: Fatal Error")
        }
        return Pair(mag.first(),starName.toString())
    }
    /**
     * Calculate fixed star positions (faster version if many stars are calculated) (UT).
     * @param star name of fixed star to search for
     * @param tjdut input time, Julian day number,  Ephemeris Time
     * @param flags bit flags indicating what kind of computation is wanted
     * @return xx stnam, int retflags
     *  - xx: 6 doubles for results
     *  - stnam: returned star name
     *  - retflags: bit flags indicating what kind of computation was done
     */
    fun fixstar2_ut(
        star: String,
        tjd_ut: Double,
        iflag: Int = 258,
    ): Triple<List<Double>, String, Int> {
        val starName = StringBuilder(star)
        val serr = StringBuilder()
        val xx = DoubleArray(6)
        val returncode = swe.fixstar2_ut(starName, tjd_ut, iflag, xx,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "fixstar2_ut: $serr")
        }
        if (returncode == ERR){
            Log.d(TAG, "fixstar2_ut: Fatal Error")
        }
        return Triple(xx.toList(), starName.toString(),returncode)
    }

    /**
     * Get fixed star magnitude.
     * @param star name of fixed star
     * @return mag, stnam
     *  - mag: returned magnitude
     *  - stnam: returned star name
     */
    fun fixstar_mag(
        star: String,
    ): Pair<Double, String> {
        val starName = StringBuilder(star)
        val serr = StringBuilder()
        val mag = DoubleArray(1)
        val returncode = swe.fixstar_mag(starName, mag,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "fixstar_mag: $serr")
        }
        if (returncode == ERR){
            Log.d(TAG, "fixstar_mag: Fatal Error")
        }
        return Pair(mag.first(),starName.toString())
    }
    /**
     * Calculate fixed star positions (UT).
     * @param star name of fixed star to search for
     * @param tjdut input time, Julian day number,  Ephemeris Time
     * @param flags bit flags indicating what kind of computation is wanted
     * @return xx stnam, int retflags
     *  - xx: 6 doubles for results
     *  - stnam: returned star name
     *  - retflags: bit flags indicating what kind of computation was done
     */
    fun fixstar_ut(
        star: String,
        tjd_ut: Double,
        iflag: Int = 258,
    ): Triple<List<Double>, String, Int> {
        val starName = StringBuilder(star)
        val serr = StringBuilder()
        val xx = DoubleArray(6)
        val returnCode = swe.fixstar_ut(starName, tjd_ut, iflag, xx,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "fixstar_ut: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "fixstar_ut: Fatal Error")
        }
        return Triple(xx.toList(), starName.toString(),returnCode)
    }

    /**
     * Calculate Gauquelin sector position of a body (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param body planet number (number in String) or fixed star name (String)
     * @param method number indicating which computation method is wanted:
     *
     *     - 0 with latitude
     *     - 1 without latitude
     *     - 2 from rising and setting times of the disc center of planet
     *     - 3 from rising and setting times of disc center, incl. refraction
     *     - 4 from rising and setting times of the disk edge of planet
     *     - 5 from rising and setting times of disk edge, incl. refraction
     * @param geopos a DoubleArray(3) containing
     *
     *     - 0: geographic longitude, in degrees (eastern positive)
     *     - 1: geographic latitude, in degrees (northern positive)
     *     - 2: geographic altitude, in meters above sea level
     * @param atpress atmospheric pressure (if 0, the default 1013.25 mbar is used)
     * @param attemp atmospheric temperature in degrees Celsius
     * @param flags bit flags for ephemeris and FLG_TOPOCTR, etc
     * @return sector: Double
     * - sector: Range(1.0,37.0). Gauquelin sectors are numbered in clockwise direction.
     */
    fun gauquelin_sector(
        tjdut:Double,
        body:String,
        method: Int,
        geopos: DoubleArray,
        atpress: Double = 0.0,
        attemp: Double,
        flags: Int
    ): Double {
        val ipl = body.toIntOrNull()
        val serr = StringBuilder()
        val dgsect = DoubleArray(1)
        val returnCode = if(ipl!= null){
           swe.gauquelin_sector(
                tjdut,ipl, StringBuilder(), flags,method,geopos,atpress,attemp,dgsect,serr)
        }else{
            swe.gauquelin_sector(
                tjdut,0, StringBuilder(body), flags,method,geopos,atpress,attemp,dgsect,serr)
        }
        if (serr.isNotEmpty()){
            Log.d(TAG, "gauquelin_sector: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "gauquelin_sector: Fatal Error")
        }
        return dgsect.first()
    }

    /**
     * Calculate ayanamsa (ET).
     * @param tjdet input time, Julian day number, Ephemeris Time
     * @return aya: ayanamsa value, without nutation
     */
    fun get_ayanamsa(tjdet:Double): Double {
        return swe.get_ayanamsa(tjdet)
    }

    /**
     * Calculate ayanamsa, extended version (ET).
     * @param tjdet input time, Julian day number, Ephemeris Time
     * @param flags ephemeris flag, etc
     * @return Pair(retflags: Int, aya: Double)
     *  - retflags: returned bit flags
     *  - aya: ayanamsa value
     */
    fun get_ayanamsa_ex(
        tjdet:Double,
        flags:Int,
    ): Pair<Int, Double> {
        val serr = StringBuilder()
        val daya = DoubleArray(1)
        val returnCode = swe.get_ayanamsa_ex(tjdet,flags,daya,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "get_ayanamsa_ex: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "get_ayanamsa_ex: Fatal Error")
        }
        return Pair(returnCode,daya.first())
    }
    /**
     * Calculate ayanamsa, extended version (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param flags ephemeris flag, etc
     * @return Pair(retflags: Int, aya: Double)
     *  - retflags: returned bit flags
     *  - aya: ayanamsa value
     */
    fun get_ayanamsa_ex_ut(
        tjdut:Double,
        flags:Int,
    ): Pair<Int, Double> {
        val serr = StringBuilder()
        val daya = DoubleArray(1)
        val returnCode = swe.get_ayanamsa_ex_ut(tjdut,flags,daya,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "get_ayanamsa_ex_ut: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "get_ayanamsa_ex_ut: Fatal Error")
        }
        return Pair(returnCode,daya.first())
    }
    /**
     * Get ayanamsa name from sidereal mode constant.
     *
     * If sidmode is not found (incorrect), returned string is empty.
     * @param sidmode
     * @return name
     */
    fun get_ayanamsa_name(sidmode:Int): String? {
        return swe.get_ayanamsa_name(sidmode)
    }
    /**
     * Calculate ayanamsa (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @return aya: ayanamsa value, without nutation
     */
    fun get_ayanamsa_ut(tjdut:Double): Double {
        return swe.get_ayanamsa_ut(tjdut)
    }

    /**
     * Find start and end date of an se1 ephemeris file after a function call.
     *
     * This can be used to find out the start and end date of an ``se1`` ephemeris file after a call of [calc].
     *
     * The function returns data from internal file structures ``sweph.fidat`` used in the last call to [calc] or [fixstar].
     * Data returned are (currently) 0 with JPL files and fixed star files.
     * Thus, the function is only useful for ephemerides of planets or asteroids that are based on ``se1`` files.
     * @param fno an integer indicating what type of file is searched:
     *
     *     - 0: planet file sepl_xxx, used for Sun etc, or jpl file
     *     - 1: moon file semo_xxx
     *     - 2: main asteroid file seas_xxx, if such an object was computed
     *     - 3: other asteroid or planetary moon file, if such object was computed
     *     - 4: star file
     * @return Triple(path: String, Pair(start: Double, end: Double), denum: Int)
     *  - path: full file path, or empty string if no data
     *  - start: start date of file
     *  - end: end date of file
     *  - denum: jpl ephemeris number 406 or 431 from which file was derived
     */
    fun get_current_file_data(fno: Int): Triple<String?, Pair<Double, Double>, Int> {
        val tfstart = DoubleArray(1)
        val tfend = DoubleArray(1)
        val denum = IntArray(1)
        val path = swe.get_current_file_data(fno, tfstart,tfend,denum)
        return Triple(path, Pair(tfstart.first(),tfend.first()), denum.first())
    }

    /**
     * Find the path of the executable or swisseph library (dll) actually in use.
     *
     * This function may fail depending the systems.
     * @return path: String
     */
    fun get_library_path(): String? {
        return swe.get_library_path()
    }

    /**
     * Calculate osculating elements (Kepler elements) and orbital periods.
     * @param tjdet input time, Julian day number, Ephemeris Time (TT)
     * @param planet identifier of planet or object
     * @param flags bit flags indicating what computation is wanted:
     *
     *     - ephemeris flag: FLG_JPLEPH, FLG_SWIEPH, FLG_MOSEPH, etc
     *     - center:
     *        - Sun: FLG_HELCTR (assumed as default) or
     *        - SS Barycentre: FLG_BARYCTR (rel. to solar system barycentre)
     *          Only possible for planets beyond Jupiter.
     *          For elements of the Moon, the calculation is geocentric.
     *     - sum all masses inside the orbit to be computed (method of
     *       Astronomical Almanac): FLG_ORBEL_AA
     *     - reference ecliptic: FLG_J2000
     * @return elements
     *  - elements: a List of 50 Double, of which:
     *     - 0: semimajor axis (a)
     *     - 1: eccentricity (e)
     *     - 2: inclination (in)
     *     - 3: longitude of ascending node (upper-case omega OM)
     *     - 4: argument of periapsis (lower-case omega om)
     *     - 5: longitude of periapsis (peri)
     *     - 6: mean anomaly at epoch (M0)
     *     - 7: true anomaly at epoch (N0)
     *     - 8: eccentric anomaly at epoch (E0)
     *     - 9: mean longitude at epoch (LM)
     *     - 10: sidereal orbital period in tropical years
     *     - 11: mean daily motion
     *     - 12: tropical period in years
     *     - 13: synodic period in days, negative for inner planets or Moon
     *     - 14: time of perihelion passage
     *     - 15: perihelion distance
     *     - 16: aphelion distance
     */
    fun get_orbital_elements(
        tjdet: Double,
        planet: Int,
        flags: Int
    ): List<Double> {
        val serr = StringBuilder()
        val dret = DoubleArray(50)
        val  returnCode = swe.get_orbital_elements(tjdet,planet,flags,dret,serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "get_orbital_elements: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "get_orbital_elements: Fatal Error")
        }
        return dret.toList()
    }

    /**
     * Get a planet or asteroid name.
     * @param planet identifier of planet or object
     * @return name found or empty string
     */
    fun get_planet_name(planet: Int): String? {
        return swe.get_planet_name(planet)
    }

    /**
     * Get current value of the tidal acceleration.
     * @return tidacc: Double
     */
    fun get_tid_acc(): Double {
        return swe.get_tid_acc()
    }

    /**
     * Provides data that are relevant for the calculation of heliacal risings and settings.
     * @param tjdut input time, Julian day number, Universal Time
     * @param geopos a DoubleArray with:
     *
     *     - 0: geographic longitude (eastern positive)
     *     - 1: geographic latitude (northern positive)
     *     - 2: altitude above sea level, in meters
     * @param atmo a DoubleArray with:
     *
     *     - 0: atmospheric pressure in mbar (hPa)
     *     - 1: atmospheric temperature in degrees Celsius
     *     - 2: relative humidity in %
     *     - 3: if >= 1, Meteorological Range (km).
     *       Between 1 and 0, total atmospheric coefficient (ktot).
     *       If = 0, the other atmospheric parameters determine the total
     *       atmospheric coefficient (ktot)
     * @param observer a DoubleArray with:
     *
     *     - 0: age of observer in years (default = 36)
     *     - 1: snellen ratio of observers eyes (default = 1 = normal)
     *     - The following parameters are only relevant if HELFLAG_OPTICAL_PARAMS
     *       is set:
     *     - 2: (0) = monocular, (1) = binocular (boolean)
     *     - 3: telescope magnification, (0) = default to naked eye (binocular),
     *       (1) = naked eye
     *     - 4: optical aperture (telescope diameter) in mm
     *     - 5: optical transmission
     * @param objname name of planet or fixed star
     * @param eventtype either:
     *
     *     - HELIACAL_RISING: morning first, for all visible planets and stars
     *     - HELIACAL_SETTING: evening last, for all visible planets and stars
     *     - EVENING_FIRST: evening first, for Mercury, Venus, Moon
     *     - MORNING_LAST: morning last, for Mercury, Venus, Moon
     * @param flags bit flags for ephemeris, and also:
     *
     *     - HELFLAG_OPTICAL_PARAMS: for optical instruments
     *     - HELFLAG_NO_DETAILS: provide date, without details
     *     - HELFLAG_VISLIM_DARK: behave as if Sun is at nadir
     *     - HELFLAG_VISLIM_NOMOON: behave as if Moon is at nadir, i.e. the Moon as
     *       a factor disturbing the observation is excluded, useful if one is not
     *       interested in the heliacal date of that particular year, but in the
     *       heliacal date of that epoch
     * @return  - dret: List of 50 Double, of which:
     *     - 0: AltO [deg] topocentric altitude of object (unrefracted)
     *     - 1: AppAltO [deg] apparent altitude of object (refracted)
     *     - 2: GeoAltO [deg] geocentric altitude of object
     *     - 3: AziO [deg] azimuth of object
     *     - 4: AltS [deg] topocentric altitude of Sun
     *     - 5: AziS [deg] azimuth of Sun
     *     - 6: TAVact [deg] actual topocentric arcus visionis
     *     - 7: ARCVact [deg] actual (geocentric) arcus visionis
     *     - 8: DAZact [deg] actual difference between object's and sun's azimuth
     *     - 9: ARCLact [deg] actual longitude difference between object and sun
     *     - 10: kact [-] extinction coefficient
     *     - 11: minTAV [deg] smallest topocentric arcus visionis
     *     - 12: TfistVR [JDN] first time object is visible, according to VR
     *     - 13: TbVR [JDN] optimum time the object is visible, according to VR
     *     - 14: TlastVR [JDN] last time object is visible, according to VR
     *     - 15: TbYallop [JDN] best time the object is visible, according to Yallop
     *     - 16: WMoon [deg] crescent width of Moon
     *     - 17: qYal [-] q-test value of Yallop
     *     - 18: qCrit [-] q-test criterion of Yallop
     *     - 19: ParO [deg] parallax of object
     *     - 20: Magn [-] magnitude of object
     *     - 21: RiseO [JDN] rise/set time of object
     *     - 22: RiseS [JDN] rise/set time of Sun
     *     - 23: Lag [JDN] rise/set time of object minus rise/set time of Sun
     *     - 24: TvisVR [JDN] visibility duration
     *     - 25: LMoon [deg] crescent length of Moon
     *     - 26: CVAact [deg]
     *     - 27: Illum [%] new
     *     - 28: CVAact [deg] new
     *     - 29: MSk [-]
     */
    fun heliacal_pheno_ut(
        tjdut:Double,
        geopos: DoubleArray,
        atmo: DoubleArray,
        observer:DoubleArray,
        objname:String,
        eventtype:Int,
        flags: Int
    ): List<Double> {
        val serr = StringBuilder()
        val darr = DoubleArray(50)
        val returnCode = swe.heliacal_pheno_ut(tjdut, geopos, atmo, observer,objname,eventtype,flags,darr, serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "heliacal_pheno_ut: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "heliacal_pheno_ut: Fatal Error")
        }
        return darr.toList()
    }

    /**
     * Find the Julian day of the next heliacal phenomenon.
     * @param tjdut input time, Julian day number, Universal Time
     * @param geopos a DoubleArray with:
     *
     *     - 0: geographic longitude (eastern positive)
     *     - 1: geographic latitude (northern positive)
     *     - 2: altitude above sea level, in meters
     * @param atmo a DoubleArray with:
     *
     *     - 0: atmospheric pressure in mbar (hPa)
     *     - 1: atmospheric temperature in degrees Celsius
     *     - 2: relative humidity in %
     *     - 3: if >= 1, Meteorological Range (km).
     *       Between 1 and 0, total atmospheric coefficient (ktot).
     *       If = 0, the other atmospheric parameters determine the total
     *       atmospheric coefficient (ktot)
     * @param observer a DoubleArray with:
     *
     *     - 0: age of observer in years (default = 36)
     *     - 1: snellen ratio of observers eyes (default = 1 = normal)
     *     - The following parameters are only relevant if HELFLAG_OPTICAL_PARAMS
     *       is set:
     *     - 2: (0) = monocular, (1) = binocular (boolean)
     *     - 3: telescope magnification, (0) = default to naked eye (binocular),
     *       (1) = naked eye
     *     - 4: optical aperture (telescope diameter) in mm
     *     - 5: optical transmission
     * @param objname name of planet or fixed star
     * @param eventtype either:
     *
     *     - HELIACAL_RISING: morning first, for all visible planets and stars
     *     - HELIACAL_SETTING: evening last, for all visible planets and stars
     *     - EVENING_FIRST: evening first, for Mercury, Venus, Moon
     *     - MORNING_LAST: morning last, for Mercury, Venus, Moon
     * @param flags bit flags for ephemeris, and also:
     *
     *     - HELFLAG_OPTICAL_PARAMS: for optical instruments
     *     - HELFLAG_NO_DETAILS: provide date, without details
     *     - HELFLAG_VISLIM_DARK: behave as if Sun is at nadir
     *     - HELFLAG_VISLIM_NOMOON: behave as if Moon is at nadir, i.e. the Moon as
     *       a factor disturbing the observation is excluded, useful if one is not
     *       interested in the heliacal date of that particular year, but in the
     *       heliacal date of that epoch
     * @return  - dret: List<Double> of 3 Julian Days, of which:
     *     - 0: start visibility
     *     - 1: optimum visibility, 0 if flags >= HELFLAG_AV
     *     - 2: end of visibility, 0 if flags >= HELFLAG_AV
     *
     */
    fun heliacal_ut(
        tjdut: Double,
        geopos: DoubleArray,
        atmo: DoubleArray,
        observer: DoubleArray,
        objname: String,
        eventtype:Int,
        flags: Int
    ): List<Double> {
        val serr = StringBuilder()
        val darr = DoubleArray(3)
        val returnCode = swe.heliacal_ut(tjdut, geopos, atmo, observer,objname,eventtype,flags,darr, serr)
        if (serr.isNotEmpty()){
            Log.d(TAG, "heliacal_ut: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "heliacal_ut: Fatal Error")
        }
        return darr.toList()
    }

    /**
     * Compute a planet heliocentric crossing over some longitude (ET).
     * @param planet planet number
     * @param x2cross longitude to search
     * @param tjdet start time of search, as Julian day number, Ephemeris Time
     * @param flags bit flags indicating what computation is wanted
     * @param backwards a boolean indicating if we search back in time
     * @return jdcross: Julian day found
     */
    fun helio_cross(
        planet: Int,
        x2cross: Double,
        tjdet: Double,
        flags:Int,
        backwards: Boolean
        ): Double {

        val jx = DoubleArray(1)
        val serr = StringBuilder()
        val dir = if(backwards) -1 else 1
        val returnCode =  swe.helio_cross(planet,x2cross,tjdet,flags,dir,jx,serr )
        if (serr.isNotEmpty()){
            Log.d(TAG, "helio_cross: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "helio_cross: Fatal Error")
        }
        return jx.first()
    }
    /**
     * Compute a planet heliocentric crossing over some longitude (UT).
     * @param planet planet number
     * @param x2cross longitude to search
     * @param tjdut start time of search, as Julian day number, Universal Time
     * @param flags bit flags indicating what computation is wanted
     * @param backwards a boolean indicating if we search back in time
     * @return jdcross: Julian day found
     */
    fun helio_cross_ut(
        planet: Int,
        x2cross: Double,
        tjdut: Double,
        flags:Int,
        backwards: Boolean
    ): Double {

        val jx = DoubleArray(1)
        val serr = StringBuilder()
        val dir = if(backwards) -1 else 1
        val returnCode =  swe.helio_cross_ut(planet,x2cross,tjdut,flags,dir,jx,serr )
        if (serr.isNotEmpty()){
            Log.d(TAG, "helio_cross_ut: $serr")
        }
        if (returnCode == ERR){
            Log.d(TAG, "helio_cross_ut: Fatal Error")
        }
        return jx.first()
    }

    /**
     * Get the name of a house method.
     * @param hsys house system identifier
     * @return hsysname: house system name, empty string if not found
     */
    fun house_name(hsys: Int): String? {
        return swe.house_name(hsys)
    }

    /**
     * Calculate house position of a body.
     * @param armc ARMC
     * @param geolat geographic latitude, in degrees (northern positive)
     * @param eps obliquity, in degrees
     * @param objcoord DoubleArray(2) for ecl. longitude and latitude of the planet, in degrees
     * @param hsys house method identifier (1 byte)
     * @return hpos: value in Range(1:13) (Gauquelin: Range(1,37)) indicating the house position
     */
    fun house_pos(
        armc: Double,
        geolat: Double,
        eps: Double,
        objcood: DoubleArray,
        hsys: Int = 'P'.code
    ): Double {
        val serr = StringBuilder()
        val position = swe.house_pos(armc, geolat, eps, hsys, objcood, serr)
        if(serr.isNotEmpty()){
            Log.d(TAG, "house_pos: $serr")
        }
        return position
    }

    /** Calculate houses cusps (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param lat geographic latitude, in degrees (northern positive)
     * @param lon geographic longitude, in degrees (eastern positive)
     * @param hsys house method: ie. 'G'.toInt(), 'P'.toInt()
     * @return (cusps), (ascmc)
     *  - cusps: List of 12 Double for cusps (except Gauquelin: 36 Double)
     *  - ascmc: List of 8 Double for additional points (@see astro.com/swisseph/swephprg.htm)
     */
    fun houses(
        tjdut: Double,
        lat: Double,
        lon: Double,
        hsys:Int = 'P'.code
    ): Pair<List<Double>, List<Double>> {
        val size =  if (hsys == 'G'.code ||hsys == 'g'.code) 37 else 13
        val cusps = DoubleArray(size)
        val ascmc = DoubleArray(10)
        val returnCode = swe.houses(tjdut, lat, lon, hsys, cusps, ascmc)
        if(returnCode==ERR){
            Log.d(TAG, "houses: ERR")
        }
        return Pair(cusps.drop(1).toList(), ascmc.take(8).toList())
    }
    /** Calculate houses cusps with ARMC.
     * @param armc ARMC
     * @param lat geographic latitude, in degrees (northern positive)
     * @param eps obliquity, in degrees
     * @param hsys house method: ie. 'G'.toInt(), 'P'.toInt()
     * @param ascmc9 optional parameter for Sunshine house system
     * @return (cusps), (ascmc)
     *  - cusps: List of 12 Double for cusps (except Gauquelin: 36 Double)
     *  - ascmc: List of 8 Double for additional points (@see astro.com/swisseph/swephprg.htm)
     */
    fun houses_armc(
        armc: Double,
        lat: Double,
        eps: Double,
        hsys:Int = 'P'.code,
        ascmc9: Double = 0.0
    ): Pair<List<Double>, List<Double>> {
        val size =  if (hsys == 'G'.code ||hsys == 'g'.code) 37 else 13
        val cusps = DoubleArray(size)
        val ascmc = DoubleArray(10)
        ascmc[9] = ascmc9
        val returnCode = swe.houses_armc(armc, lat, eps, hsys, cusps, ascmc)
        if(returnCode==ERR){
            Log.d(TAG, "houses_armc: ERR")
        }
        return Pair(cusps.drop(1).toList(), ascmc.take(8).toList())
    }
    /** Calculate houses cusps and their speeds with ARMC.
     * @param armc ARMC
     * @param lat geographic latitude, in degrees (northern positive)
     * @param eps obliquity, in degrees
     * @param hsys house method: ie. 'G'.toInt(), 'P'.toInt()
     * @param ascmc9 optional parameter for Sunshine house system
     * @return (cusps), (ascmc), (cuspsspeed), (ascmcspeed)
     *  - cusps: List of 12 Double for cusps (except Gauquelin: 36 Double)
     *  - ascmc: List of 8 Double for additional points (@see astro.com/swisseph/swephprg.htm)
     *  - cuspsspeed: List of 12 Double for cusps speeds
     *  - ascmcspeed: List of 8 Double for speeds of additional points
     */
    fun houses_armc_ex2(
        armc: Double,
        lat: Double,
        eps: Double,
        hsys:Int = 'P'.code,
        ascmc9: Double = 0.0
    ): Pair<List<Double>, List<Double>> {
        val size =  if (hsys == 'G'.code ||hsys == 'g'.code) 37 else 13
        val cusps = DoubleArray(size)
        val ascmc = DoubleArray(10)
        val cusps_speed = DoubleArray(size)
        val ascmc_speed = DoubleArray(10)
        ascmc[9] = ascmc9
        val serr = StringBuilder()
        val returnCode = swe.houses_armc_ex2(armc, lat, eps, hsys, cusps, ascmc, cusps_speed,ascmc_speed,serr)
        if(returnCode==ERR){
            Log.d(TAG, "houses_armc_ex2: ERR")
        }
        if (serr.isNotEmpty()){
            Log.d(TAG, "houses_armc_ex2: $serr")
        }
        return Pair(cusps.drop(1).toList(), ascmc.take(8).toList())
    }

    /** Calculate houses cusps (extended) (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param lat geographic latitude, in degrees (northern positive)
     * @param lon geographic longitude, in degrees (eastern positive)
     * @param hsys house method: ie. 'G'.toInt(), 'P'.toInt()
     * @param flags ephemeris flag, etc
     * @return (cusps), (ascmc)
     *  - cusps: List of 12 Double for cusps (except Gauquelin: 36 Double)
     *  - ascmc: List of 8 Double for additional points (@see astro.com/swisseph/swephprg.htm)
     */
    fun houses_ex(
        tjdut: Double,
        lat: Double,
        lon: Double,
        hsys:Int = 'P'.code,
        flags: Int
    ): Pair<List<Double>, List<Double>> {
        val size =  if (hsys == 'G'.code ||hsys == 'g'.code) 37 else 13
        val cusps = DoubleArray(size)
        val ascmc = DoubleArray(10)
        val returnCode = swe.houses_ex(tjdut,flags, lat, lon, hsys, cusps, ascmc)
        if(returnCode==ERR){
            Log.d(TAG, "houses_ex: ERR")
        }
        return Pair(cusps.drop(1).toList(), ascmc.take(8).toList())
    }
    /** Calculate houses cusps and cusps speeds (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param lat geographic latitude, in degrees (northern positive)
     * @param lon geographic longitude, in degrees (eastern positive)
     * @param hsys house method: ie. 'G'.toInt(), 'P'.toInt()
     * @param flags ephemeris flag, etc
     * @return (cusps), (ascmc)
     *  - cusps: List of 12 Double for cusps (except Gauquelin: 36 Double)
     *  - ascmc: List of 8 Double for additional points (@see astro.com/swisseph/swephprg.htm)
     *  - cuspsspeed: List of 12 Double for cusps speeds
     *  - ascmcspeed: List of 8 Double for speeds of additional points
     */
    fun houses_ex2(
        tjdut: Double,
        lat: Double,
        lon: Double,
        hsys:Int = 'P'.code,
        flags: Int
    ): Pair<List<Double>, List<Double>> {
        val size =  if (hsys == 'G'.code ||hsys == 'g'.code) 37 else 13
        val cusps = DoubleArray(size)
        val ascmc = DoubleArray(10)
        val cusps_speed = DoubleArray(size)
        val ascmc_speed = DoubleArray(10)
        val serr = StringBuilder()
        val returnCode = swe.houses_ex2(tjdut,flags, lat, lon, hsys, cusps, ascmc, cusps_speed, ascmc_speed, serr)
        if(returnCode==ERR){
            Log.d(TAG, "houses_ex2: ERR")
        }
        if (serr.isNotEmpty()){
            Log.d(TAG, "houses_ex2: $serr")
        }
        return Pair(cusps.drop(1).toList(), ascmc.take(8).toList())
    }

    /**
     * Convert ET Julian day number to UTC.
     * @param tjdet Julian day number in ET (TT)
     * @param cal calendar flag, either GREG_CAL or JUL_CAL
     * @return Pair(listOf(year: Int, month: Int, day: Int, hour: Int, mins: Int), second: Double)
     */
    fun jdet_to_utc(tjdet: Double, cal: Int): Pair<List<Int>, Double> {
        val date = IntArray(5)
        val dsec = DoubleArray(1)
        swe.jdet_to_utc(tjdet, cal, date, dsec)
        return Pair(date.toList(), dsec.first())
    }
    /**
     * Convert UT1 Julian day number to UTC.
     * @param tjdut Julian day number, in UT (UT1)
     * @param cal calendar flag, either GREG_CAL or JUL_CAL
     * @return Pair(listOf(year: Int, month: Int, day: Int, hour: Int, mins: Int), second: Double)
     */
    fun jdut1_to_utc(tjdut: Double, cal: Int): Pair<List<Int>, Double> {
        val date = IntArray(5)
        val dsec = DoubleArray(1)
        swe.jdut1_to_utc(tjdut, cal, date, dsec)
        return Pair(date.toList(), dsec.first())
    }

    /**
     * Translate local apparent time (LAT) to local mean time (LMT).
     * @param tjdlat Julian day number, local apparent time
     * @param geolon geographic longitude, in degrees (eastern positive)
     * @return tjdlmt: returned Julian day number, local mean time
     */
    fun lat_to_lmt(tjdlat: Double, geolon: Double): Double {
        val tjdlmtPointer = DoubleArray(1)
        val serr = StringBuilder()
        swe.lat_to_lmt(tjdlat,geolon,tjdlmtPointer,serr)
        if(serr.isNotEmpty()){
            Log.d(TAG, "lat_to_lmt: $serr")
        }
        return tjdlmtPointer.first()
    }
    /**
     * Translate local mean time (LMT) to local apparent time (LAT).
     * @param tjdlmt Julian day number, local mean time
     * @param geolon geographic longitude, in degrees (eastern positive)
     * @return tjdlat: returned Julian day number, local apparent time
     */
    fun lmt_to_lat(tjdlmt: Double, geolon: Double): Double {
        val tjdlatPointer = DoubleArray(1)
        val serr = StringBuilder()
        swe.lmt_to_lat(tjdlmt,geolon,tjdlatPointer,serr)
        if(serr.isNotEmpty()){
            Log.d(TAG, "lmt_to_lat: $serr")
        }
        return tjdlatPointer.first()
    }

    /**
     * Calculate attributes of a lunar eclipse (UTC).
     * @param tjdut input time, Julian day number, Universal Time
     * @param geopos a sequence with:
     *
     *     - geographic longitude, in degrees (eastern positive)
     *     - geographic latitude, in degrees (northern positive)
     *     - geographic altitude above sea level, in meters
     * @param flags ephemeris flag, etc
     * @return retflag, (attr)
     *  - retflag: returned bit flags:
     *     - 0 if there is no eclipse
     *     - SE_ECL_TOTAL or ECL_PENUMBRAL or ECL_PARTIAL
     *  - attr: list of 20 Double, of which:
     *     - 0: umbral magnitude at tjd
     *     - 1: penumbral magnitude
     *     - 2: ?
     *     - 3: ?
     *     - 4: azimuth of moon at tjd
     *     - 5: true altitude of moon above horizon at tjd
     *     - 6: apparent altitude of moon above horizon at tjd
     *     - 7: distance of moon from opposition in degrees
     *     - 8: eclipse magnitude (equals attr[0])
     *     - 9: saros series number (if available, otherwise -99999999)
     *     - 10: saros series member number (if available, otherwise -99999999)
     */
    fun lun_eclipse_how(
        tjdut: Double,
        geopos: DoubleArray,
        flags: Int = 2
    ): Pair<Int, List<Double>> {
        val serr = StringBuilder()
        val attr = DoubleArray(20)
        val retflag = swe.lun_eclipse_how(tjdut, flags, geopos, attr, serr)
        if (retflag == ERR){
            Log.d(TAG, "lun_eclipse_how: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "lun_eclipse_how: $serr")
        }
        return Pair(retflag, attr.toList())
    }

    /**
     * Find the next lunar eclipse globally (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param flags ephemeris flag
     * @param ecltype Int for eclipse type wanted:
     *
     *     - ECL_TOTAL ECL_PARTIAL ECL_PENUMBRAL
     *     - ECL_ALLTYPES_LUNAR or 0 for any type
     * @param backwards boolean, set to True to search back in time
     * @return retflag, (tret)
     *  - retflag: returned bit flag:
     *     - ECL_TOTAL ECL_PARTIAL ECL_PENUMBRAL
     *  - tret: List of 10 Double, of which:
     *     - 0: time of maximum eclipse
     *     - 1: ?
     *     - 2: time of partial phase begin (indices consistent with solar eclipses)
     *     - 3: time of partial phase end
     *     - 4: time of totality begin
     *     - 5: time of totality end
     *     - 6: time of penumbral phase begin
     *     - 7: time of penumbral phase end
     */
    fun lun_eclipse_when(
        tjdut: Double,
        flags: Int = 2,
        ecltype: Int = 0,
        backwards: Boolean = false
    ): Pair<Int, List<Double>> {
        val serr = StringBuilder()
        val tret = DoubleArray(10)
        val backwardsInt = if(backwards) 1 else 0
        val retflag = swe.lun_eclipse_when(tjdut,flags,ecltype,tret, backwardsInt,serr)
        if(retflag ==ERR){
            Log.d(TAG, "lun_eclipse_when: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "lun_eclipse_when: $serr")
        }
        return Pair(retflag, tret.toList())
    }
    /**
     * Find the next lunar eclipse observable from a given geographic position (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param geopos a DoubleArray with:
     *
     *     - geographic longitude, in degrees (eastern positive)
     *     - geographic latitude, in degrees (northern positive)
     *     - geographic altitude, in meters above sea level
     * @param flags ephemeris flag
     * @param backwards boolean, set to True to search back in time
     * @return retflag, (tret), (attr)
     *  - retflag: returned bit flag:
     *     - ECL_TOTAL ECL_PARTIAL ECL_PENUMBRAL
     *  - tret: List of 10 Double, of which:
     *     - 0: time of maximum eclipse
     *     - 1: ?
     *     - 2: time of partial phase begin (indices consistent with solar eclipses)
     *     - 3: time of partial phase end
     *     - 4: time of totality begin
     *     - 5: time of totality end
     *     - 6: time of penumbral phase begin
     *     - 7: time of penumbral phase end
     *     - 8: time of moonrise, if it occurs during the eclipse
     *     - 9: time of moonset, if it occurs during the eclipse
     *  - attr: list of 20 Double, of which:
     *     - 0: umbral magnitude at tjd
     *     - 1: penumbral magnitude
     *     - 2: ?
     *     - 3: ?
     *     - 4: azimuth of moon at tjd
     *     - 5: true altitude of moon above horizon at tjd
     *     - 6: apparent altitude of moon above horizon at tjd
     *     - 7: distance of moon from opposition in degrees (separation angle)
     *     - 8: umbral magnitude at tjd (equals attr[0])
     *     - 9: saros series number (if available; otherwise -99999999)
     *     - 10: saros series member number (if available; otherwise -99999999)
     */
    fun lun_eclipse_when_loc(
        tjd_start: Double,
        geopos: DoubleArray,
        flags: Int =2,
        backwards: Boolean = false
    ): Triple<Int, List<Double>, List<Double>> {
        val serr = StringBuilder()
        val tret = DoubleArray(10)
        val attr = DoubleArray(20)
        val backwardsInt = if(backwards) 1 else 0
        val retflag = swe.lun_eclipse_when_loc(tjd_start,flags,geopos,tret,attr, backwardsInt,serr)
        if(retflag ==ERR){
            Log.d(TAG, "lun_eclipse_when_loc: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "lun_eclipse_when_loc: $serr")
        }
        return Triple(retflag, tret.toList(), attr.toList())
    }

    /**
     * Find the next occultation of a planet or star by the moon globally (UT)
     * @param tjdut input time, Julian day number, Universal Time
     * @param body planet identifier or star name (str)
     * @param flags ephemeris flag, eventually ECL_ONE_TRY, etc
     * @param ecltype bit flags for eclipse type wanted:
     *
     *     - ECL_CENTRAL ECL_NONCENTRAL ECL_TOTAL ECL_ANNULAR ECL_PARTIAL
     *     - ECL_ANNULAR_TOTAL (equals ECL_HYBRID)
     *     - 0 for any type
     * @param backwards boolean, set to True to search back in time
     * @return retflag, (tret)
     *  - retflag: returned bit flag:
     *      - 0 if no occultation or eclipse found
     *     - ECL_TOTAL or ECL_ANNULAR or ECL_PARTIAL or ECL_ANNULAR_TOTAL    - ECL_CENTRAL
     *     - ECL_NONCENTRAL
     *  - tret: list of 10 Double, of which:
     *      - 0: time of maximum occultation/eclipse
     *      - 1: time when occultation takes place at local apparent noon
     *      - 2: time of occultation begin
     *      - 3: time of occultation end
     *      - 4: time of of totality begin
     *      - 5: time of totality end
     *      - 6: time of center line begin
     *      - 7: time of center line end
     *      - 8: time when annular-total occultation becomes total
     *      - 9: time when annular-total occultation becomes annular again
     */
    fun lun_occult_when_glob(
        tjdut: Double,
        body: String,
        flags: Int = 2,
        ecltype: Int = 0,
        backwards: Boolean = false
    ): Pair<Int, List<Double>> {
        val serr = StringBuilder()
        val tret = DoubleArray(10)
        val backwardsInt = if(backwards) 1 else 0
        val bodyNumber = body.toIntOrNull()

        val retflag = if (bodyNumber!=null){
            swe.lun_occult_when_glob(tjdut,bodyNumber,StringBuilder(),flags,ecltype,tret, backwardsInt,serr)
        }else{
            swe.lun_occult_when_glob(tjdut,0,StringBuilder(body),flags,ecltype,tret, backwardsInt,serr)
        }

        if(retflag ==ERR){
            Log.d(TAG, "lun_occult_when_glob: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "lun_occult_when_glob: $serr")
        }
        return Pair(retflag, tret.toList())
    }

    /**
     * Find next occultation of a planet or star by the moon for a given geographic position (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param body planet identifier or star name (str)
     * @param geopos a sequence with:
     *
     *     - geographic longitude, in degrees (eastern positive)
     *     - geographic latitude, in degrees (northern positive)
     *     - geographic altitude above sea level, in meters
     * @param flags ephemeris flag, eventually ECL_ONE_TRY, etc
     * @param backwards boolean, set to True to search back in time
     * @return retflags, (tret), (attr)
     *  - retflags: returned bit flags:
     *     - 0 if no occultation or eclipse found
     *     - ECL_TOTAL or ECL_ANNULAR or ECL_PARTIAL,
     *     - ECL_VISIBLE, ECL_MAX_VISIBLE, ECL_1ST_VISIBLE, ECL_2ND_VISIBLE,
     *       ECL_3RD_VISIBLE, ECL_4TH_VISIBLE
     *  - tret: List of 10 Double, of which:
     *     - 0: time of maximum occultation
     *     - 1: time of first contact
     *     - 2: time of second contact
     *     - 3: time of third contact
     *     - 4: time of fourth contact
     *  - attr: List of 20 Double, of which:
     *     - 0: fraction of planet diameter covered by moon (magnitude)
     *     - 1: ratio of lunar diameter to planet one
     *     - 2: fraction of planet disc covered by moon (obscuration)
     *     - 3: diameter of core shadow in km
     *     - 4: azimuth of planet at tjd
     *     - 5: true altitude of planet above horizon at tjd
     *     - 6: apparent altitude of planet above horizon at tjd
     *     - 7: elongation of moon in degrees (separation angle)
     */
    fun lun_occult_when_loc(
        tjdut: Double,
        body: String,
        geopos: DoubleArray,
        flags: Int = 2,
        backwards: Boolean = false
    ): Triple<Int, List<Double>, List<Double>> {
        val serr = StringBuilder()
        val tret = DoubleArray(10)
        val attr = DoubleArray(20)

        val backwardsInt = if(backwards) 1 else 0
        val bodyNumber = body.toIntOrNull()

        val retflag = if (bodyNumber!=null){
            swe.lun_occult_when_loc(tjdut,bodyNumber,StringBuilder(),flags,geopos,tret,attr, backwardsInt,serr)
        }else{
            swe.lun_occult_when_loc(tjdut,0,StringBuilder(body),flags,geopos,tret,attr, backwardsInt,serr)
        }
        if(retflag ==ERR){
            Log.d(TAG, "lun_occult_when_loc: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "lun_occult_when_loc: $serr")
        }
        return Triple(retflag, tret.toList(),attr.toList())
    }

    /**
     * Find where a lunar occultation is central or maximal (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param body planet identifier (int) or star name (str)
     * @param flags ephemeris flag
     * @return retflags, (geopos), (attr)
     *  - retflags: returned bit flags:
     *     - 0 if there is no occultation at tjd
     *     - ECL_TOTAL
     *     - ECL_ANNULAR
     *     - ECL_TOTAL | ECL_CENTRAL
     *     - ECL_TOTAL | ECL_NONCENTRAL
     *     - ECL_ANNULAR | ECL_CENTRAL
     *     - ECL_ANNULAR | ECL_NONCENTRAL
     *     - ECL_PARTIAL
     *  - geopos: List of 10 Double, of which:
     *     - 0: geographic longitude of central line
     *     - 1: geographic latitude of central line
     *     - 2: geographic longitude of northern limit of umbra
     *     - 3: geographic latitude of northern limit of umbra
     *     - 4: geographic longitude of southern limit of umbra
     *     - 5: geographic latitude of southern limit of umbra
     *     - 6: geographic longitude of northern limit of penumbra
     *     - 7: geographic latitude of northern limit of penumbra
     *     - 8: geographic longitude of southern limit of penumbra
     *     - 9: geographic latitude of southern limit of penumbra
     *  - attr: List of 20 Double, of which:
     *     - 0: fraction of object's diameter covered by moon (magnitude)
     *     - 1: ratio of lunar diameter to object's diameter
     *     - 2: fraction of object's disc covered by moon (obscuration)
     *     - 3: diameter of core shadow in km
     *     - 4: azimuth of object at tjd
     *     - 5: true altitude of object above horizon at tjd
     *     - 6: apparent altitude of object above horizon at tjd
     *     - 7: angular distance of moon from object in degrees
     */
    fun lun_occult_where(
        tjdut: Double,
        body: String,
        flags: Int = 2,
    ): Triple<Int, List<Double>, List<Double>> {
        val attr = DoubleArray(20)
        val serr = StringBuilder()
        val bodyNumber = body.toIntOrNull()
        val geopos = DoubleArray(2)
        val retflag = if (bodyNumber!=null){
            swe.lun_occult_where(tjdut,bodyNumber,StringBuilder(),flags,geopos,attr,serr)
        }else{
            swe.lun_occult_where(tjdut,0,StringBuilder(body),flags,geopos,attr,serr)
        }
        if(retflag ==ERR){
            Log.d(TAG, "lun_occult_where: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "lun_occult_where: $serr")
        }
        return Triple(retflag, geopos.toList(),attr.toList())
    }

    /**
     * Compute Moon crossing over some longitude (ET).
     * @param x2cross longitude to search
     * @param tjdet start time of search, Julian day number, Ephemeris Time
     * @param flags bit flags indicating what computation is wanted
     * @return jd_cross: Julian day number found
     */
    fun mooncross(
        x2cross: Double,
        tjdet:Double,
        flags:Int = 2
    ): Double {
        val serr = StringBuilder()
        val jx = swe.mooncross(x2cross,tjdet,flags,serr)
        if(jx < tjdet){
            Log.d(TAG, "mooncross: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "mooncross: $serr")
        }
        return jx
    }

    /**
     * Compute next Moon crossing over node, by finding zero latitude crossing (ET).
     * @param tjdet start time of search, Julian day number, Ephemeris Time
     * @param flags bit flags indicating what computation is wanted
     * @return jd_cross: Double, xlon: Double, xlat: Double
     *  - jd_cross: Julian day number found
     *  - xlon: Moon longitude found
     *  - xlat: Moon latitude found
     *
     */
    fun mooncross_node(tjdet: Double,flags:Int = 2): Triple<Double, Double, Double> {
        val serr = StringBuilder()
        val xlon = DoubleArray(1)
        val xlat = DoubleArray(1)
        val jd_cross = swe.mooncross_node(tjdet,flags,xlon,xlat,serr)
        if(jd_cross<tjdet){
            Log.d(TAG, "mooncross_node: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "mooncross_node: $serr")
        }
        return Triple(jd_cross,xlon.first(),xlat.first())
    }
    /**
     * Compute next Moon crossing over node, by finding zero latitude crossing (UT).
     * @param tjdut start time of search, Julian day number, Universal Time
     * @param flags bit flags indicating what computation is wanted
     * @return jd_cross: Double, xlon: Double, xlat: Double
     *  - jd_cross: Julian day number found
     *  - xlon: Moon longitude found
     *  - xlat: Moon latitude found
     *
     */
    fun mooncross_node_ut(tjdut: Double,flags:Int = 2): Triple<Double, Double, Double> {
        val serr = StringBuilder()
        val xlon = DoubleArray(1)
        val xlat = DoubleArray(1)
        val jd_cross = swe.mooncross_node_ut(tjdut,flags,xlon,xlat,serr)
        if(jd_cross<tjdut){
            Log.d(TAG, "mooncross_node_ut: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "mooncross_node_ut: $serr")
        }
        return Triple(jd_cross,xlon.first(),xlat.first())
    }
    /**
     * Compute Moon crossing over some longitude (UT).
     * @param x2cross longitude to search
     * @param tjdet start time of search, Julian day number, Universal Time
     * @param flags bit flags indicating what computation is wanted
     * @return jd_cross: Julian day number found
     */
    fun mooncross_ut(
        x2cross: Double,
        tjdut:Double,
        flags:Int = 2
    ): Double {
        val serr = StringBuilder()
        val jx = swe.mooncross_ut(x2cross,tjdut,flags,serr)
        if(jx < tjdut){
            Log.d(TAG, "mooncross_ut: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "mooncross_ut: $serr")
        }
        return jx
    }

    /**
     * Calculate planetary nodes and apsides (ET).
     * @param tjdet input time, Julian day number, Ephemeris Time
     * @param planet identifer of planet or object
     * @param method bit flags NODBIT_MEAN, NODBIT_OSCU, NODBIT_OSCU_BAR, NODBIT_FOPOINT
     * @param flags bit flags indicating what type of computation is wanted
     * @return (xnasc)(xndsc)(xperi)(xaphe)
     *  - xnasc: List of 6 Double for ascending node
     *  - xndsc: List of 6 Double for descending node
     *  - xperi: List of 6 Double for perihelion
     *  - xaphe: List of 6 Double for aphelion
     */
    fun nod_aps(
        tjdet: Double,
        planet: Int,
        method: Int,
        flags: Int = 258
    ): List<List<Double>> {
        val serr = StringBuilder()
        val xnasc = DoubleArray(6)
        val xndsc = DoubleArray(6)
        val xperi = DoubleArray(6)
        val xaphe = DoubleArray(6)
        val returnCode = swe.nod_aps(tjdet,planet,flags,method,
            xnasc, xndsc, xperi, xaphe, serr)
        if(returnCode==ERR){
            Log.d(TAG, "nod_aps: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "nod_aps: $serr")
        }
        return listOf(xnasc.toList(),xndsc.toList(),xperi.toList(),xaphe.toList())
    }
    /**
     * Calculate planetary nodes and apsides (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param planet identifer of planet or object
     * @param method bit flags NODBIT_MEAN, NODBIT_OSCU, NODBIT_OSCU_BAR, NODBIT_FOPOINT
     * @param flags bit flags indicating what type of computation is wanted
     * @return (xnasc)(xndsc)(xperi)(xaphe)
     *  - xnasc: List of 6 Double for ascending node
     *  - xndsc: List of 6 Double for descending node
     *  - xperi: List of 6 Double for perihelion
     *  - xaphe: List of 6 Double for aphelion
     */
    fun nod_aps_ut(
        tjdut: Double,
        planet: Int,
        method: Int,
        flags: Int = 258
    ): List<List<Double>> {
        val serr = StringBuilder()
        val xnasc = DoubleArray(6)
        val xndsc = DoubleArray(6)
        val xperi = DoubleArray(6)
        val xaphe = DoubleArray(6)
        val returnCode = swe.nod_aps_ut(tjdut,planet,flags,method,
            xnasc, xndsc, xperi, xaphe, serr)
        if(returnCode==ERR){
            Log.d(TAG, "nod_aps_ut: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "nod_aps_ut: $serr")
        }
        return listOf(xnasc.toList(),xndsc.toList(),xperi.toList(),xaphe.toList())
    }

    /**
     * Calculate the maximum possible distance, the minimum possible distance, and the current true distance of planet, the EMB, or an asteroid.
     * @param tjdet input time, Julian day number, Ephemeris Time
     * @param planet identifier of planet or object
     * @param flags bit flags indicating what computation is wanted:
     *
     *     - ephemeris flag
     *     - optional heliocentric flag FLG_HELIOCTR
     * @return max_dist: Double, min_dist: Double, true_dist: Double
     *  - max_dist: maximum distance
     *  - min_dist: minimum_distance
     *  - true_dist: true distance
     */
    fun orbit_max_min_true_distance(
        tjdet: Double,
        planet: Int,
        flags: Int
    ): Triple<Double, Double, Double> {
        val serr = StringBuilder()
        val dmax = DoubleArray(1)
        val dmin = DoubleArray(1)
        val dtrue = DoubleArray(1)
        val returnCode = swe.orbit_max_min_true_distance(tjdet,planet,flags,dmax,dmin,dtrue,serr)
        if(returnCode==ERR){
            Log.d(TAG, "orbit_max_min_true_distance: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "orbit_max_min_true_distance: $serr")
        }
        return Triple(dmax.first(), dmin.first(),dtrue.first())
    }

    /**
     * Calculate planetary phenomena (ET).
     * @param tjdet input time, Julian day number, Ephemeris Time
     * @param planet object identifier
     * @param flags ephemeris flag, etc
     * @return  - attr: List of 20 Double, of which:
     *     - 0: phase angle (earth-planet-sun)
     *     - 1: phase (illuminated fraction of disc)
     *     - 2: elongation of planet
     *     - 3: apparent diameter of disc
     *     - 4: apparent magnitude
     *     - 5: geocentric horizontal parallax (Moon)
     */
    fun pheno(
        tjdet: Double,
        planet: Int,
        flags: Int = 2
    ): List<Double> {
        val serr = StringBuilder()
        val attr = DoubleArray(20)
        val returnCode = swe.pheno(tjdet, planet, flags, attr, serr)
        if(returnCode==ERR){
            Log.d(TAG, "pheno: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "pheno: $serr")
        }
        return attr.toList()
    }
    /**
     * Calculate planetary phenomena UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param planet object identifier
     * @param flags ephemeris flag, etc
     * @return  - attr: List of 20 Double, of which:
     *     - 0: phase angle (earth-planet-sun)
     *     - 1: phase (illuminated fraction of disc)
     *     - 2: elongation of planet
     *     - 3: apparent diameter of disc
     *     - 4: apparent magnitude
     *     - 5: geocentric horizontal parallax (Moon)
     */
    fun pheno_ut(
        tjdut: Double,
        planet: Int,
        flags: Int = 2
    ): List<Double> {
        val serr = StringBuilder()
        val attr = DoubleArray(20)
        val returnCode = swe.pheno_ut(tjdut, planet, flags, attr, serr)
        if(returnCode==ERR){
            Log.d(TAG, "pheno_ut: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "pheno_ut: $serr")
        }
        return attr.toList()
    }

    /**
     * Calculate midpoint (in radians).
     */
    fun rad_midp(x:Double, y:Double): Double {
        return swe.rad_midp(x,y)
    }

    /**
     * Normalization of any radian number to the range [0;2*pi].
     */
    fun radnorm(x:Double): Double {
        return swe.radnorm(x)
    }

    /**
     * Calculate true altitude from apparent altitude, or vice-versa.
     * @param alt altitude of object above geometric horizon in degrees,
     * where geometric horizon = plane perpendicular to gravity
     * @param atpress atmospheric pressure in mbar (hPa)
     * @param attemp atmospheric temperature in degrees Celsius
     * @param flag either TRUE_TO_APP or APP_TO_TRUE
     * @return retalt: converted altitude
     */
    fun refrac(
        alt: Double,
        atpress: Double,
        attemp: Double,
        flag: Int
    ): Double {
        return swe.refrac(alt,atpress,attemp,flag)
    }
    /**
     * Calculate true altitude from apparent altitude, or vice-versa (extended).
     * @param alt altitude of object above geometric horizon in degrees,
     * where geometric horizon = plane perpendicular to gravity
     * @param geoalt altitude of observer above sea level, in meters,
     * @param atpress atmospheric pressure in mbar (hPa)
     * @param attemp atmospheric temperature in degrees Celsius
     * @param lapserate dattemp/dgeoalt [deg K/m]
     * @param flag either TRUE_TO_APP or APP_TO_TRUE
     * @return
     * - retalt: converted altitude
     * - dret: list of 4 Double:
     *     - 0: true altitude if possible, otherwise input value
     *     - 1: apparent altitude if possible, otherwise input value
     *     - 2: refraction
     *     - 3: dip of the horizon
     */
    fun refrac_extended(
        alt: Double,
        geoalt: Double,
        atpress: Double,
        attemp: Double,
        lapserate: Double,
        flag: Int
    ): Pair<Double, List<Double>> {
        val dret = DoubleArray(4)
        val retalt = swe.refrac_extended(alt,geoalt, atpress,attemp,lapserate,flag,dret)
        return Pair(retalt, dret.toList())
    }

    /**
     * Calculate times of rising, setting and meridian transits.
     * @param tjdut input time, Julian day number, Universal Time
     * @param body planet identifier (int) or fixed star name (str)
     * @param rsmi bit flag for rise, set, or one of the two meridian transits, etc
     * @param geopos a sequence for:
     *
     *   - 0: geographic longitude, in degrees (eastern positive)
     *   - 1: geographic latitude, in degrees (northern positive)
     *   - 2: geographic altitude, in meters above sea level
     * @param atpress atmospheric pressure in mbar/hPa
     * @param attemp atmospheric temperature in degrees Celsius
     * @param flags ephemeris flags etc
     * @return
     *  - res: integer indicating:
     *     - (0) event found
     *     - (-2) event not found because the object is circumpolar
     *  - tret: list of 10 Double, of which:
     *     - 0: tjd of event
     */
    fun rise_trans(
        tjdut:Double,
        body: String,
        rsmi: Int,
        geopos: DoubleArray,
        atpress : Double = 0.0,
        attemp: Double = 0.0,
        flags: Int = 2
    ): Pair<Int, List<Double>> {
        val serr = StringBuilder()
        val bodyNumber = body.toIntOrNull()
        val tret  = DoubleArray(10)
        val retflag = if (bodyNumber!=null){
            swe.rise_trans(tjdut,bodyNumber,StringBuilder(),flags,rsmi,geopos, atpress,attemp,tret,serr)
        }else{
            swe.rise_trans(tjdut,0,StringBuilder(body),flags,rsmi,geopos, atpress,attemp,tret,serr)
        }
        if (retflag == ERR){
            Log.d(TAG, "rise_trans: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "rise_trans: $serr")
        }
        return Pair(retflag, tret.toList())
    }
    /**
     * Calculate times of rising, setting and meridian transits (with altitude).
     * @param tjdut input time, Julian day number, Universal Time
     * @param body planet identifier (int) or fixed star name (str)
     * @param rsmi bit flag for rise, set, or one of the two meridian transits, etc
     * @param geopos a sequence for:
     *
     *   - 0: geographic longitude, in degrees (eastern positive)
     *   - 1: geographic latitude, in degrees (northern positive)
     *   - 2: geographic altitude, in meters above sea level
     * @param atpress atmospheric pressure in mbar/hPa
     * @param attemp atmospheric temperature in degrees Celsius
     * @param horhgt height of local horizon in degrees (where body rises or sets)
     * @param flags ephemeris flags etc
     * @return
     *  - res: integer indicating:
     *     - (0) event found
     *     - (-2) event not found because the object is circumpolar
     *  - tret: list of 10 Double, of which:
     *     - 0: tjd of event
     */
    fun rise_trans_true_hor(
        tjdut:Double,
        body: String,
        rsmi: Int,
        geopos: DoubleArray,
        atpress : Double = 0.0,
        attemp: Double = 0.0,
        horhgt: Double = 0.0,
        flags: Int = 2
    ): Pair<Int, List<Double>> {
        val serr = StringBuilder()
        val bodyNumber = body.toIntOrNull()
        val tret  = DoubleArray(10)
        val retflag = if (bodyNumber!=null){
            swe.rise_trans_true_hor(tjdut,bodyNumber,StringBuilder(),flags,rsmi,geopos, atpress,attemp,horhgt,tret,serr)
        }else{
            swe.rise_trans_true_hor(tjdut,0,StringBuilder(body),flags,rsmi,geopos, atpress,attemp,horhgt,tret,serr)
        }
        if (retflag == ERR){
            Log.d(TAG, "rise_trans_true_hor: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "rise_trans_true_hor: $serr")
        }
        return Pair(retflag, tret.toList())
    }

    /**
     * Set a fixed Deltat T value.
     *
     * This function allows the user to set a fixed Delta T value that will be returned by [deltat] or [deltat_ex]. The same Delta T value will then be used by [calc_ut], eclipse functions, heliacal functions, and all functions that require UT as input time.
     *
     * In order to return to automatic Delta T, call this function with the following value: set_delta_t_userdef(DELTAT_AUTOMATIC)
     * @param acc
     *
     */
    fun set_delta_t_userdef(acc: Double){
        swe.set_delta_t_userdef(acc)
    }

    /**
     * Set ephemeris files path.
     * @param path It is possible to pass None as path, which is equivalent to an empty string.
     */
    fun set_ephe_path(path:String = "/data/user/0/com.example.swelib/files/ephe"){
        swe.set_ephe_path(path)
    }

    /**
     * Set name of JPL ephemeris file.
     *
     * This file must reside in the ephemeris path you are using for all your ephemeris files.
     */
    fun set_jpl_file(name: String){
        swe.set_jpl_file(name)
    }

    /**
     * Set lapse rate.
     */
    fun set_lapse_rate(lrate: Double){
        swe.set_lapse_rate(lrate)
    }

    /**
     * Set sidereal mode.
     *
     * If a predefined mode is wanted, the variable sid_mode has to be set, while t0 and ayan_t0 are not considered, i.e. can be 0.
     */
    fun set_sid_mode(sidmode: Int, t0:Double=0.0, ayan_t0: Double=0.0){
        swe.set_sid_mode(sidmode,t0,ayan_t0)
    }

    /**
     * Set value of the tidal acceleration.
     * @param acc the values possible are:
     *
     *     - TIDAL_DE200
     *     - TIDAL_DE403
     *     - TIDAL_DE404
     *     - TIDAL_DE405
     *     - TIDAL_DE406
     *     - TIDAL_DE421
     *     - TIDAL_DE422
     *     - TIDAL_DE430
     *     - TIDAL_DE431
     *     - TIDAL_DE441
     *     - TIDAL_26
     *     - TIDAL_STEPHENSON_2016
     *     - TIDAL_DEFAULT (equals TIDAL_DE431)
     *     - TIDAL_MOSEPH (equals TIDAL_DE404)
     *     - TIDAL_SWIEPH (equals TIDAL_DEFAULT)
     *     - TIDAL_JPLEPH (equals TIDAL_DEFAULT)
     */
    fun set_tid_acc(acc: Double){
        swe.set_tid_acc(acc)
    }

    /**
     * Set topocentric parameters.
     * @param lon geographic longitude, in degrees (eastern positive)
     * @param lat geographic latitude, in degrees (northern positive)
     * @param alt geographic altitude in meters above sea level
     */
    fun set_topo(lon:Double, lat:Double,alt:Double = 0.0){
        swe.set_topo(lon,lat,alt)
    }

    /**
     * Calculate sidereal time (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @return sidtime: Double
     */
    fun sidtime(tjdut: Double): Double {
        return swe.sidtime(tjdut)
    }
    /**
     * Calculate sidereal time, given obliquity and nutation (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param eps obliquity, in degrees
     * @param nutation nutation in longitude, in degrees
     * @return sidtime: Double
     */
    fun sidtime(tjdut: Double, eps:Double, nutation: Double): Double {
        return swe.sidtime0(tjdut,eps, nutation)
    }

    /**
     * Calculate attributes of a solar eclipse for a given geographic position and time.
     * @param tjdut input time, Julian day number, Universal Time
     * @param geopos a sequence with:
     *
     *   - geographic longitude, in degrees (eastern positive)
     *   - geographic latitude, in degrees (northern positive)
     *   - geographic altitude above sea level, in meters
     * @param flags ephemeris flag, etc
     * @return
     *  - retflags: returned bit flags:
     *     - 0 if no eclipse is visible at position
     *     - ECL_TOTAL ECL_ANNULAR ECL_PARTIAL
     *  - attr: list of 20 Double, of which:
     *     - 0: fraction of solar diameter covered by moon
     *     - 1: ratio of lunar diameter to solar one
     *     - 2: fraction of solar disc covered by moon (obscuration)
     *     - 3: diameter of core shadow in km
     *     - 4: azimuth of sun at tjd
     *     - 5: true altitude of sun above horizon at tjd
     *     - 6: apparent altitude of sun above horizon at tjd
     *     - 7: elongation of moon in degrees (separation angle)
     *     - 8: magnitude acc. to NASA (equals attr[0] for partial and attr[1] for      annular and total eclipses)
     *     - 9: saros series number (if available, otherwise -99999999)
     *     - 10: saros series member number (if available, otherwise -99999999)
     */
    fun sol_eclipse_how(
        tjdut: Double,
        geopos: DoubleArray,
        flags: Int = 2
    ): Pair<Int, List<Double>> {
        val attr = DoubleArray(20)
        val serr = StringBuilder()
        val returnCode = swe.sol_eclipse_how(tjdut,flags,geopos,attr, serr)
        if(returnCode==ERR){
            Log.d(TAG, "sol_eclipse_how: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "sol_eclipse_how: $serr")
        }
        return Pair(returnCode,attr.toList())
    }

    /**
     * Find the next solar eclipse globally (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param flags ephemeris flag
     * @param ecltype bit flags for eclipse type wanted:
     *
     *    - ECL_CENTRAL ECL_NONCENTRAL ECL_TOTAL ECL_ANNULAR ECL_PARTIAL
     *    - ECL_ANNULAR_TOTAL (equals ECL_HYBRID)
     *    - ECL_ALLTYPES_SOLAR or 0 for any type
     * @param backwards boolean, set to True to search back in time
     * @return
     *  - res: returned bit flags:
     *     - ECL_TOTAL ECL_ANNULAR ECL_PARTIAL ECL_ANNULAR_TOTAL
     *     - ECL_CENTRAL
     *     - ECL_NONCENTRAL
     *
     *  - tret: list of 10 Double, of which:
     *     - 0: time of maximum eclipse
     *     - 1: time when eclipse takes place at local apparent noon
     *     - 2: time of eclipse begin
     *     - 3: time of eclipse end
     *     - 4: time of totality begin
     *     - 5: time of totality end
     *     - 6: time of center line begin
     *     - 7: time of center line end
     *     - 8: time when annular-total eclipse becomes total
     *     - 9: time when annular-total eclipse becomes annular again
     */
    fun sol_eclipse_when_glob(
        tjdut: Double,
        flags: Int,
        ecltype: Int,
        backwards: Boolean
    ): Pair<Int, List<Double>> {
        val serr = StringBuilder()
        val tret = DoubleArray(10)
        val returnCode = swe.sol_eclipse_when_glob(tjdut,flags,ecltype,tret,if (backwards) 1 else 0,serr)
        if(returnCode==ERR){
            Log.d(TAG, "sol_eclipse_when_glob: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "sol_eclipse_when_glob: $serr")
        }
        return Pair(returnCode,tret.toList())
    }
    /**
     * Find the next solar eclipse for a given geographic position (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param geopos a sequence with:
     *
     *     - geographic longitude, in degrees (eastern positive)
     *     - geographic latitude, in degrees (northern positive)
     *     - geographic altitude above sea level, in meters
     * @param flags ephemeris flag
     * @param backwards boolean, set to True to search back in time
     * @return
     *  - retflags: returned bit flags:
     *     - ECL_TOTAL or ECL_ANNULAR or ECL_PARTIAL,
     *     - ECL_VISIBLE, ECL_MAX_VISIBLE, ECL_1ST_VISIBLE, ECL_2ND_VISIBLE,
     *       ECL_3RD_VISIBLE, ECL_4TH_VISIBLE
     *  - tret: list of 10 Double, of which:
     *     - 0: time of maximum eclipse
     *     - 1: time of first contact
     *     - 2: time of second contact
     *     - 3: time of third contact
     *     - 4: time of fourth contact
     *     - 5: time of sunrise between first and fourth contact
     *     - 6: time of sunset between first and fourth contact
     *  - attr: list of 20 Double, of which:
     *     - 0: fraction of solar diameter covered by moon; with total/annular
     *       eclipse, it results in magnitude acc. to IMCCE.
     *     - 1: ratio of lunar diameter to solar one
     *     - 2: fraction of solar disc covered by moon (obscuration)
     *     - 3: diameter of core shadow in km
     *     - 4: azimuth of sun at tjd
     *     - 5: true altitude of sun above horizon at tjd
     *     - 6: apparent altitude of sun above horizon at tjd
     *     - 7: elongation of moon in degrees (separation angle)
     *     - 8: magnitude acc. to NASA (equals attr[0] for partial and attr[1] for      annular and total eclipses)
     *     - 9: saros series number (if available, otherwise -99999999)
     *     - 10: saros series member number (if available, otherwise -99999999)
     */
    fun sol_eclipse_when_loc(
        tjdut: Double,
        geopos: DoubleArray,
        flags: Int,
        backwards: Boolean
    ): Triple<Int, List<Double>, List<Double>> {
        val serr = StringBuilder()
        val tret = DoubleArray(10)
        val attr = DoubleArray(20)
        val returnCode = swe.sol_eclipse_when_loc(tjdut,flags,geopos,tret,attr,if (backwards) 1 else 0,serr)
        if(returnCode==ERR){
            Log.d(TAG, "sol_eclipse_when_loc: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "sol_eclipse_when_loc: $serr")
        }
        return Triple(returnCode,tret.toList(),attr.toList())
    }

    /**
     * Find where a solar eclipse is central or maximal (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @param flags ephemeris flag
     * @return
     *  - retflags: returned bit flags:
     *     - ECL_TOTAL
     *     - ECL_ANNULAR
     *     - ECL_TOTAL | ECL_CENTRAL
     *     - ECL_TOTAL | ECL_NONCENTRAL
     *     - ECL_ANNULAR | ECL_CENTRAL
     *     - ECL_ANNULAR | ECL_NONCENTRAL
     *     - ECL_PARTIAL
     *  - geopos: list of 10 Double for longitude/latitude of:
     *     - 0: geographic longitude of central line
     *     - 1: geographic latitude of central line
     *     - 2: geographic longitude of northern limit of umbra
     *     - 3: geographic latitude of northern limit of umbra
     *     - 4: geographic longitude of southern limit of umbra
     *     - 5: geographic latitude of southern limit of umbra
     *     - 6: geographic longitude of northern limit of penumbra
     *     - 7: geographic latitude of northern limit of penumbra
     *     - 8: geographic longitude of southern limit of penumbra
     *     - 9: geographic latitude of southern limit of penumbra
     *  - attr: list of 20 Double, of which:
     *     - 0: fraction of solar diameter covered by moon; with total/annular
     *       eclipse, it results in magnitude acc. to IMCCE.
     *     - 1: ratio of lunar diameter to solar one
     *     - 2: fraction of solar disc covered by moon
     *     - 3: diameter of core shadow in km
     *     - 4: azimuth of sun at tjd
     *     - 5: true altitude of sun above horizon at tjd
     *     - 6: apparent altitude of sun above horizon at tjd
     *     - 7: elongation of moon in degrees (separation angle)
     *     - 8: magnitude acc. to NASA (equals attr[0] for partial and attr[1] for      annular and total eclipses)
     *     - 9: saros series number (if available, otherwise -99999999)
     *     - 10: saros series member number (if available, otherwise -99999999)
     */
    fun sol_eclipse_where(
        tjdut: Double,
        flags: Int = 2
    ): Triple<Int, List<Double>, List<Double>> {
        val geopos = DoubleArray(2)
        val attr = DoubleArray(20)
        val serr = StringBuilder()
        val returnCode = swe.sol_eclipse_where(tjdut,flags,geopos,attr, serr)
        if(returnCode==ERR){
            Log.d(TAG, "sol_eclipse_where: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "sol_eclipse_where: $serr")
        }
        return Triple(returnCode, geopos.toList(),attr.toList())
    }

    /**
     * Compute next Sun crossing over some longitude (ET).
     * @param x2cross longitude to search
     * @param tjdet start time of search, Julian day number, Ephemeris Time
     * @param flags bit flags indicating what computation is wanted
     * @return jdcross: Julian day number found
     */
    fun solcross(
        x2cross: Double,
        tjdet: Double,
        flags: Int
    ): Double {
        val serr = StringBuilder()
        val jx = swe.solcross(x2cross,tjdet,flags,serr)
        if(jx<tjdet){
            Log.d(TAG, "solcross: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "solcross: $serr")
        }
        return jx
    }
    /**
     * Compute next Sun crossing over some longitude (ET).
     * @param x2cross longitude to search
     * @param tjdut start time of search, Julian day number, Ephemeris Time
     * @param flags bit flags indicating what computation is wanted
     * @return jdcross: Julian day number found
     */
    fun solcross_ut(
        x2cross: Double,
        tjdut: Double,
        flags: Int
    ): Double {
        val serr = StringBuilder()
        val jx = swe.solcross_ut(x2cross,tjdut,flags,serr)
        if(jx<tjdut){
            Log.d(TAG, "solcross_ut: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "solcross_ut: $serr")
        }
        return jx
    }

    /**
     * Provide sign or nakshatra, degree, minutes, seconds and fraction of second from decimal degree. Can also round to seconds, minutes, degrees.
     * @param degree position in decimal degrees
     * @param roundflag bit flags combination indicating how to round:
     *
     *     - 0: dont round
     *     - SPLIT_DEG_ROUND_SEC: round to seconds
     *     - SPLIT_DEG_ROUND_MIN: round to minutes
     *     - SPLIT_DEG_ROUND_DEG: round to degrees
     *     - SPLIT_DEG_ZODIACAL: with zodiac sign number
     *     - SPLIT_DEG_NAKSHATRA: with nakshatra number
     *     - SPLIT_DEG_KEEP_SIGN: dont round to next zodiac sign/nakshatra
     *     - SPLIT_DEG_KEEP_DEG: dont round to next degree
     * @return Triple(listOf(degree, minute, second), degree of second, sign)
     *  - deg: returned degree
     *  - min: returned minutes
     *  - sec: returned seconds
     *  - secfr: returned fraction of second
     *  - sign: returned sign/nakshatra number
     */
    fun split_deg(degree: Double, roundflag:Int): Triple<List<Int>, Double, Int> {
        val dms = IntArray(3)
        val dsecfr = DoubleArray(1)
        val isign = IntArray(1)
        swe.split_deg(degree, roundflag,dms,dsecfr,isign)
        return Triple(dms.toList(),dsecfr.first(),isign.first())
    }

    /**
     * Calculate equation of time (UT).
     * @param tjdut input time, Julian day number, Universal Time
     * @return e: difference between local apparent time and local mean time in days
     */
    fun time_equa(tjdut:Double): Double {
        val e = DoubleArray(1)
        val serr = StringBuilder()
        swe.time_equ(tjdut,e,serr)
        if(serr.isNotEmpty()){
            Log.d(TAG, "time_equa: $serr")
        }
        return e.first()
    }

    /**
     * Transform local time to UTC or UTC to local time.
     * @param offset timezone offset in hours (east of greenwich positive)
     * @return Pair(List(year, month, day, hour, second), second: Double)
     */
    fun utc_time_zone(year: Int, month: Int,day: Int,hour: Int,minutes: Int,seconds: Double,offset:Double): Pair<List<Int>, Double> {
        val ymdhm_out = IntArray(5)
        val dsec_out = DoubleArray(1)
        swe.utc_time_zone(year,month,day,hour,minutes,seconds,offset,ymdhm_out,dsec_out)
        return Pair(ymdhm_out.toList(),dsec_out.first())
    }

    /**
     * Find the limiting visual magnitude in dark skies.
     * @param tjdut input time, Julian day number, Universal Time
     * @param geopos a sequence with:
     *
     *     - 0: geographic longitude (eastern positive)
     *     - 1: geographic latitude (northern positive)
     *     - 2: altitude above sea level, in meters
     * @param atmo a sequence with:
     *     - 0: atmospheric pressure in mbar (hPa)
     *     - 1: atmospheric temperature in degrees Celsius
     *     - 2: relative humidity in %
     *     - 3: if >= 1, Meteorological Range (km).
     *       Between 1 and 0, total atmospheric coefficient (ktot).
     *       If = 0, the other atmospheric parameters determine the total
     *       atmospheric coefficient (ktot)
     * @param observer a sequence with:
     *
     *     - 0: age of observer in years (default = 36)
     *     - 1: snellen ratio of observers eyes (default = 1 = normal)
     *     - The following parameters are only relevant if HELFLAG_OPTICAL_PARAMS
     *       is set:
     *     - 2: (0) = monocular, (1) = binocular (boolean)
     *     - 3: telescope magnification, (0) = default to naked eye (binocular),
     *       (1) = naked eye
     *     - 4: optical aperture (telescope diameter) in mm
     *     - 5: optical transmission
     * @param objname name of planet or fixed star
     * @param flags bit flags for ephemeris, and also:
     *
     *     - HELFLAG_OPTICAL_PARAMS: for optical instruments
     *     - HELFLAG_NO_DETAILS: provide date, without details
     *     - HELFLAG_VISLIM_DARK: behave as if Sun is at nadir
     *     - HELFLAG_VISLIM_NOMOON: behave as if Moon is at nadir, i.e. the Moon as
     *       a factor disturbing the observation is excluded, useful if one is not
     *       interested in the heliacal date of that particular year, but in the
     *       heliacal date of that epoch
     * @return 
     *  - res: result:
     *     - (-2): object is below horizon
     *     - (0): OK, photopic vision
     *     - (1): OK, scotopic vision
     *     - (2): OK, near limit photopic/scotopic vision
     *  - dret: a list of 10 Double, of which:
     *     - 0: limiting visual magnitude (if > magnitude of object, then the
     *       object is visible)
     *     - 1: altitude of object
     *     - 2: azimuth of object
     *     - 3: altitude of sun
     *     - 4: azimuth of sun
     *     - 5: altitude of moon
     *     - 6: azimuth of moon
     *     - 7: magnitude of object
     */
    fun vis_limit_mag(
        tjdut: Double,
        geopos: DoubleArray,
        atmo: DoubleArray,
        observer: DoubleArray,
        objname: String,
        flags: Int
    ): Pair<Int, List<Double>> {
        val serr = StringBuilder()
        val dret = DoubleArray(10)
        val returnCode = swe.vis_limit_mag(tjdut, geopos,atmo,observer,StringBuilder(objname),flags,dret, serr)
        if(returnCode==ERR){
            Log.d(TAG, "vis_limit_mag: ERR")
        }
        if(serr.isNotEmpty()){
            Log.d(TAG, "vis_limit_mag: $serr")
        }
        return Pair(returnCode, dret.toList())
    }






}

