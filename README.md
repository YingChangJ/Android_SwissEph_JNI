# Android_SwissEph_JNI
Really simple android app, using swiss ephemeris through jni.
Before everything, you need to add the ephe files to a place where the app could find them, check with the MainActivity.kt to get this location.

* This is a simple example demonstrating the use of the Swiss Ephemeris library in a real Android application.
* No library needed, just JNI.
* You can find more information and details about Swiss Ephemeris at https://www.astro.com.
* The code is based on krymlov's JNI implementation, available at https://github.com/krymlov, and includes a class wrapped to simplify usage.
* Using JNI directly offers greater flexibility, making it easy to add or remove functionality in C code.
* This demo attempts to mimic the functionality of pyswisseph (the Swiss Ephemeris library in Python) as closely as possible, making it user-friendly.
* Please note that this demo does not provide 100% API coverage; the PlaCalc function, which is now a legacy feature and not used in Swiss Ephemeris, is not included.

Last word: the interface of krymlov's original work should be better. What I did, is make the use of the swisseph functions easy to use, as pyswisseph.
