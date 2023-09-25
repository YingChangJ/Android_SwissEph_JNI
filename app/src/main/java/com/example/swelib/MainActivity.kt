package com.example.swelib

import android.os.Bundle
import androidx.activity.ComponentActivity
import androidx.activity.compose.setContent
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.Surface
import androidx.compose.material3.Text
import androidx.compose.ui.Modifier
import androidx.compose.ui.platform.LocalContext
import com.example.swelib.ui.theme.MyApplicationTheme
import java.time.ZonedDateTime

class MainActivity : ComponentActivity() {
    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)

        setContent {
            MyApplicationTheme {
                val context = LocalContext.current
                val filesDir = context.filesDir
                Surface(
                    modifier = Modifier.fillMaxSize(),
                    color = MaterialTheme.colorScheme.background
                ) {
                    Column {
                        val swe = SwissEphKt
                        //With some modification in the SwissEph(change the initial filedir there,
                        //no need to set this ephe_path each time the app run in Android
                        swe.set_ephe_path(filesDir.toString())
                        //time now in ZonedDateTime
                        val time = ZonedDateTime.now()
                        //julianDay_ut of this instant
                        val julianDayUT = swe.julday(time)
                        //calculate the sun at this time, 258 = FLG_SPEED|FLG_SWIEPH
                        val (xxSun,retflagsSun) = swe.calc_ut(julianDayUT, SUN,258)
                        //calculate a fixstar, say Spica, at this time
                        val (xxSpica, stnamSpica, retflagsSpica) = swe.fixstar_ut("Spica",julianDayUT)
                        //calculate the houses' cusps
                        val(cusps,ascmc) = swe.houses(julianDayUT,39.9,116.9)
                        //great, with all the information, it is easily to create a beautiful chart with canvas
                        swe.close()
                        //Sun's longitude, latitude, distance, speed in long., speed in lat., and speed in dist.
                        Text(text = "Sun:\n$xxSun")
                        //Did the phone find my ephemeris? If it is 258 then, yes, it could find my ephemeris and use them
                        Text(text = "Flags:\n $retflagsSun")
                        //Spica's longitude, latitude, distance, speed in long., speed in lat., and speed in dist.
                        Text(text = "Spica:\n $xxSpica")
                        //Houses Cusps
                        Text(text = "Houses (Placidus):\n $cusps")

                    //if you have no idea about where the ephe file should be in the android phone, check this out with these lines.
                    //Text(text = filesDir.toString())
                    }
                    
                }
            }
        }
    }
}


