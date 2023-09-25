package com.example.swelib


class SwissEph(ephePath: String = "/data/user/0/com.example.swelib/files/ephe") {
    init {
        System.loadLibrary("swelib")
        set_ephe_path(ephePath)
    }


    /**
     * computes the horizontal coordinates (azimuth and altitude)
     *
     * @param tjdut Julian day number, Universal Time
     * @param flag either from ecliptical coord (ECL2HOR) or equatorial (EQU2HOR)
     * @param geopos DoubleArray[3] with:
     *      - 0: geographic longitude, in degrees (eastern positive)
     *      - 1: geographic latitude, in degrees (northern positive)
     *      - 2: geographic altitude, in meters above sea level
     * @param atpress atmospheric pressure in mbar (hPa), if set to 0.0 pressure will estimated by geopos
     * @param attemp atmospheric temperature in degrees Celsius, 15 degree here as mean temp
     * @param xin: a sequence with:
     *     - ECL2HOR: ecl. longitude, ecl. latitude, distance (0)
     *     - EQU2HOR: right ascension, declination, distance (1)
     * @return DoubleArray[3] with: float azimuth, true_altitude, apparent_altitude
     *  - azimuth: position degree, measured from south point to west
     *  - true_altitude: true altitude above horizon in degrees
     *  - apparent_altitude: apparent (refracted) altitude above horizon in   degrees
     */
    external fun azalt(
        tjdut: Double,
        flag: Int,
        geopos: DoubleArray?,
        atpress: Double,
        attemp: Double,
        xin: DoubleArray?
    ): DoubleArray

    external fun azalt_rev(
        tjd_ut: Double,
        calc_flag: Int,
        geopos: DoubleArray?,
        xin: DoubleArray?,
        xout: DoubleArray?
    )

    external fun heliacal_ut(
        tjdstart: Double,
        dgeo: DoubleArray?,
        datm: DoubleArray?,
        dobs: DoubleArray?,
        objectname: String?,
        event_type: Int,
        helflag: Int,
        dret: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun heliacal_pheno_ut(
        tjd_ut: Double,
        dgeo: DoubleArray?,
        datm: DoubleArray?,
        dobs: DoubleArray?,
        objectname: String?,
        event_type: Int,
        helflag: Int,
        darr: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun vis_limit_mag(
        tjdut: Double,
        dgeo : DoubleArray?,
        datm: DoubleArray?,
        dobs: DoubleArray?,
        objectname: StringBuilder?,
        helflag: Int,
        dret: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun heliacal_angle(
        var0: Double,
        var2: DoubleArray?,
        var3: DoubleArray?,
        var4: DoubleArray?,
        var5: Int,
        var6: Double,
        var8: Double,
        var10: Double,
        var12: Double,
        var14: Double,
        var16: DoubleArray?,
        var17: StringBuilder?
    ): Int

    external fun topo_arcus_visionis(
        var0: Double,
        var2: DoubleArray?,
        var3: DoubleArray?,
        var4: DoubleArray?,
        var5: Int,
        var6: Double,
        var8: Double,
        var10: Double,
        var12: Double,
        var14: Double,
        var16: Double,
        var18: DoubleArray?,
        var19: StringBuilder?
    ): Int

    external fun set_astro_models(var0: StringBuilder?, var1: Int)

    external fun get_astro_models(var0: StringBuilder?, var1: StringBuilder?, var2: Int)

    external fun version(): String?

    external fun get_library_path(): String?

    external fun calc(
        tjd_et: Double,
        ipl: Int,
        flag: Int,
        result: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun calc_ut(
        tjd_ut: Double,
        ipl: Int,
        flag: Int,
        result: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun calc_pctr(
        tjd: Double,
        ipl: Int,
        iplctr: Int,
        iflag: Int,
        xxret: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun fixstar(
        star: StringBuilder?,
        tjd_et: Double,
        iflag: Int,
        xx: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun fixstar_ut(
        star: StringBuilder?,
        tjd_ut: Double,
        iflag: Int,
        xx: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun fixstar_mag(
        star: StringBuilder?,
        mag: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun fixstar2(
        star: StringBuilder?,
        tjd_et: Double,
        iflag: Int,
        xx: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun fixstar2_ut(
        star: StringBuilder?,
        tjd_ut: Double,
        iflag: Int,
        xx: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun fixstar2_mag(
        star: StringBuilder?,
        mag: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun close()


    external fun set_jpl_file(fname: String?)

    external fun get_planet_name(ipl: Int): String?

    external fun set_topo(geolon: Double, geolat: Double, altitude: Double)

    external fun set_sid_mode(sid_mode: Int, t0: Double, ayan_t0: Double)

    external fun get_ayanamsa_ex(
        tjd_ut: Double,
        ephe_flag: Int,
        daya: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun get_ayanamsa_ex_ut(
        tjd_ut: Double,
        ephe_flag: Int,
        daya: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun get_ayanamsa(tjd_et: Double): Double

    external fun get_ayanamsa_ut(tjd_ut: Double): Double

    external fun get_ayanamsa_name(isidmode: Int): String?

    external fun get_current_file_data(
        ifno: Int,
        tfstart: DoubleArray?,
        tfend: DoubleArray?,
        denum: IntArray?
    ): String?

    external fun date_conversion(
        year: Int,
        month: Int,
        day: Int,
        hour: Double,
        calendar: Char,                 /* calendar ‘g’[regorian]|’j’[ulian] */
        tjdPointer: DoubleArray?
    ): Int

    external fun julday(
        year: Int,
        month: Int,
        day: Int,
        hour: Double,
        gregFlag: Int
    ): Double

    external fun revjul(
        tjd: Double,
        gregflag: Int,              /* Gregorian calendar: 1, Julian calendar: 0 */
        yearmonthdayPointer: IntArray?,     /* target addresses for year, etc. */
        hourPointer: DoubleArray?
    )
    /**
     * Convert UTC to julian day.
     *
     * int year, int month, int day, int hour, int minutes, float seconds, int cal=GREG_CAL
     * @param cal GREG_CAL:1 or JUL_CAL:0
     * @return DoubleArray[3] with: float azimuth, true_altitude, apparent_altitude
     *  - azimuth: position degree, measured from south point to west
     *  - true_altitude: true altitude above horizon in degrees
     *  - apparent_altitude: apparent (refracted) altitude above horizon in   degrees
     */
    external fun utc_to_jd(
        iyear: Int,
        imonth: Int,
        iday: Int,
        ihour: Int,
        imin: Int,
        dsec : Double,
        cal: Int,              /* Gregorian calendar: 1, Julian calendar: 0 */
        dret: DoubleArray?,         /* return array, two doubles:
                                    * dret[0] = Julian day in ET (TT)
                                    * dret[1] = Julian day in UT (UT1) */
        serr: StringBuilder?        /* error string */
    ): Int

    external fun jdet_to_utc(tjd_et: Double, gregflag: Int, date: IntArray?, dsec: DoubleArray?)

    external fun jdut1_to_utc(tjd_ut: Double, gregflag: Int, date: IntArray?, dsec: DoubleArray?)

    external fun utc_time_zone(
        iyear: Int,
        imonth: Int,
        iday: Int,
        ihour : Int,
        imin: Int,
        dsec: Double,
        d_timezone: Double,
        ymdhm_out: IntArray?,
        dsec_out: DoubleArray?
    )

    external fun houses(
        tjd_ut: Double,
        geolat: Double,
        geolon: Double,
        hsys: Int,
        cusps: DoubleArray?,
        ascmc: DoubleArray?
    ): Int

    external fun houses_ex(
        tjd_ut: Double,
        iflag: Int,
        geolat: Double,
        geolon: Double,
        hsys: Int,
        cusps: DoubleArray?,
        ascmc: DoubleArray?
    ): Int

    external fun houses_ex2(
        var0: Double,
        iflag: Int,
        geolat: Double,
        geolon: Double,
        hsys: Int,
        cusps: DoubleArray?,
        ascmc: DoubleArray?,
        cusp_speed: DoubleArray?,
        ascmc_speed: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun houses_armc(
        armc: Double,
        geolat: Double,
        eps: Double,
        hsys: Int,
        cusps: DoubleArray?,
        ascmc: DoubleArray?
    ): Int

    external fun houses_armc_ex2(
        armc: Double,
        geolat: Double,
        eps: Double,
        hsys: Int,
        cusps: DoubleArray?,
        ascmc: DoubleArray?,
        cusp_speed: DoubleArray?,
        ascmc_speed: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun house_pos(
        armc: Double,
        geolat: Double,
        eps: Double,
        hsys: Int,
        xpin: DoubleArray?,
        serr: StringBuilder?
    ): Double
    external fun solcross_ut(
        x2cross: Double, jd_ut:Double, flag:Int, serr:StringBuilder?
    ):Double
    external fun house_name(hsys: Int): String?

    external fun gauquelin_sector(
        tjd_ut: Double,
        ipl: Int,
        starname: StringBuilder?,
        iflag: Int,
        imeth: Int,
        geopos: DoubleArray?,
        atpress: Double,
        attemp: Double,
        dgsect: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun sol_eclipse_where(
        tjd_ut: Double,
        ifl: Int,
        geopos: DoubleArray?,
        attr: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun lun_occult_where(
        tjd_ut: Double,
        ipl: Int,
        starname: StringBuilder?,
        ifl: Int,
        geopos: DoubleArray?,
        attr: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun sol_eclipse_how(
        tjd_ut: Double,
        ifl: Int,
        geopos: DoubleArray?,
        attr: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun sol_eclipse_when_loc(
        tjd_start: Double,
        ifl: Int,
        geopos: DoubleArray?,
        tret: DoubleArray?,
        attr: DoubleArray?,
        backward: Int,
        serr: StringBuilder?
    ): Int

    external fun lun_occult_when_loc(
        tjd_start: Double,
        ipl: Int,
        starname: StringBuilder?,
        ifl: Int,
        geopos: DoubleArray?,
        tret: DoubleArray?,
        attr: DoubleArray?,
        backward: Int,
        serr: StringBuilder?
    ): Int

    external fun sol_eclipse_when_glob(
        tjd_start: Double,
        ifl: Int,
        ifltype: Int,
        tret: DoubleArray?,
        backward: Int,
        serr: StringBuilder?
    ): Int

    external fun lun_occult_when_glob(
        tjd_start: Double,
        ipl: Int,
        starname: StringBuilder?,
        ifl: Int,
        ifltype: Int,
        tret: DoubleArray?,
        backward: Int,
        serr: StringBuilder?
    ): Int

    external fun lun_eclipse_how(
        tjd_ut: Double,
        ifl: Int,
        geopos: DoubleArray?,
        attr: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun lun_eclipse_when(
        tjd_start: Double,
        ifl: Int,
        ifltype: Int,
        tret: DoubleArray?,
        backward: Int,
        serr: StringBuilder?
    ): Int

    external fun lun_eclipse_when_loc(
        tjd_start: Double,
        ifl: Int,
        geopos: DoubleArray?,
        tret: DoubleArray?,
        attr: DoubleArray?,
        backward: Int,
        serr: StringBuilder?
    ): Int

    external fun pheno(
        tjd_et: Double,
        ipl: Int,
        iflag: Int,
        attr: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun pheno_ut(
        tjd_ut: Double,
        ipl: Int,
        iflag: Int,
        attr: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun refrac(inalt: Double, atpress: Double, attemp: Double, calc_flag: Int): Double

    external fun refrac_extended(
        inalt: Double,
        atpress: Double,
        geoalt: Double,
        attemp: Double,
        lapse_rate: Double,
        calc_flag: Int,
        dret: DoubleArray?
    ): Double

    external fun set_lapse_rate(lrate: Double)



    external fun rise_trans_true_hor(
        tjd_ut: Double,
        ipl: Int,
        starname: StringBuilder?,
        epheflag: Int,
        rsmi: Int,
        geopos: DoubleArray?,
        atpress: Double,
        attemp: Double,
        horhgt: Double,
        tret: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun rise_trans(
        tjd_ut: Double,
        ipl: Int,
        starname: StringBuilder?,
        epheflag: Int,
        rsmi: Int,
        geopos: DoubleArray?,
        atpress: Double,
        attemp: Double,
        tret: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun nod_aps(
        tjd_et: Double,
        ipl: Int,
        iflag: Int,
        method: Int,
        xnasc: DoubleArray?,
        xndsc: DoubleArray?,
        xperi: DoubleArray?,
        xaphe: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun nod_aps_ut(
        tjd_ut: Double,
        ipl: Int,
        iflag: Int,
        method: Int,
        xnasc: DoubleArray?,
        xndsc: DoubleArray?,
        xperi: DoubleArray?,
        xaphe: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun get_orbital_elements(
        tjd_et: Double,
        ipl: Int,
        iflag: Int,
        dret: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun orbit_max_min_true_distance(
        tjd_et: Double,
        ipl: Int,
        iflag: Int,
        dmax: DoubleArray?,
        dmin: DoubleArray?,
        dtrue: DoubleArray?,
        serr: StringBuilder?
    ): Int

    external fun deltat(tjd: Double): Double
    /**
     * DeltaT (ET - UT) in days
     *
     */
    external fun deltat_ex(tdj: Double, epheFlag: Int, serr: StringBuilder?): Double

    //9.4. Mean solar time versus True solar time:
    /** Equation of Time
     * The function returns the difference between local apparent and local mean time in days.
     * E = LAT - LMT
     * Input variable tjd is UT.
     */
    external fun time_equ(
        tjd: Double,
        ePointer: DoubleArray?,
        serrPointer: StringBuilder?
    ): Int  //return OK:0; ERR:-1
    /* converts Local Mean Time (LMT) to Local Apparent Time (LAT) */
    /* tjd_lmt and tjd_lat are a Julian day number
    * geolon is geographic longitude, where eastern longitudes are positive,
    * western ones negative */
    external fun lmt_to_lat(
        tjdLmt: Double,
        geoLon: Double,
        tjdLatPointer: DoubleArray?,
        serrPointer: StringBuilder?
    ): Int

    external fun lat_to_lmt(
        tjdLat: Double,
        geoLon: Double,
        tjdLmtPointer: DoubleArray?,
        serrPointer: StringBuilder?
    ): Int
    
    external fun sidtime0(tjd_ut: Double, eps: Double, nut: Double): Double

    external fun sidtime(tjd_ut: Double): Double

    external fun set_ephe_path(ephePath:String)
    
    external fun set_interpolate_nut(var0: Int)

    external fun cotrans(xpo: DoubleArray?, xpn: DoubleArray?, eps: Double)

    external fun cotrans_sp(xpo: DoubleArray?, xpn: DoubleArray?, eps: Double)

    external fun get_tid_acc(): Double

    external fun set_tid_acc(t_acc: Double)

    external fun set_delta_t_userdef(t_acc: Double)

    external fun degnorm(x: Double): Double

    external fun radnorm(x: Double): Double

    external fun rad_midp(x1: Double, x2: Double): Double

    external fun deg_midp(x1: Double, x2: Double): Double

    external fun split_deg(
        ddeg: Double,
        roundflag: Int,
        dms: IntArray?,
        dsecfr : DoubleArray?,
        int32 : IntArray?
    )

    external fun solcross_ut_jieqi_list_bazi(
        x12: DoubleArray,
        year: Int,
        serr: StringBuilder?
    )

    /**
     * 反推八字
     * 甲寅 is number 1, 1-60
     * for month and time, 甲/寅 is 1, 1-12
     */
    external fun convertBaziToDates(
        yearGanZhi: Int,
        month: Int,
        dayGanZhi: Int,
        hour: Int,
        longitude: Double,
        rangeMin: Int = 0,
        rangeMax: Int = 2100,
        tjd_ut: DoubleArray?
    ):Int
    external fun helio_cross(
        ipl: Int,
        x2cross: Double,
        tjd_et:Double,
        iflag:Int,
        dir:Int,
        jx:DoubleArray?,
        serr:StringBuilder?
    ): Int
    external fun helio_cross_ut(
        ipl: Int,
        x2cross: Double,
        tjd_ut:Double,
        iflag:Int,
        dir:Int,
        jx:DoubleArray?,
        serr:StringBuilder?
    ): Int
    external fun solcross(
        x2cross:Double,
        tjd_et:Double,
        iflag:Int,
        serr: StringBuilder?
    ):Double
    external fun mooncross(
        x2cross:Double,
        tjd_et:Double,
        iflag:Int,
        serr: StringBuilder?
    ):Double
    external fun mooncross_ut(
        x2cross:Double,
        tjd_ut:Double,
        iflag:Int,
        serr: StringBuilder?
    ):Double
    external fun mooncross_node(
        tjd_et:Double,
        iflag:Int,
        xlon: DoubleArray?,
        xlat: DoubleArray?,
        serr: StringBuilder?
    ):Double
    external fun mooncross_node_ut(
        tjd_ut:Double,
        iflag:Int,
        xlon: DoubleArray?,
        xlat: DoubleArray?,
        serr: StringBuilder?
    ):Double
}

