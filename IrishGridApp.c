
secs = secs + 1;
/* if GPS fix is not 100% show only the GPS status*/
if (SUUNTO_GPS_STATE < 100) 
{
  RESULT = SUUNTO_GPS_STATE;
    prefix="GPS";
    postfix="%";
}
/* show the Lattitude in DD.MM format, for lack of  jave.round fuction suunto.mod had to be used*/
else if (secs < 8) 
{
 /*
  RESULT = (Suunto.abs(SUUNTO_GPS_LATITUDE) - Suunto.mod(Suunto.abs(SUUNTO_GPS_LATITUDE),1) + Suunto.mod(Suunto.abs(SUUNTO_GPS_LATITUDE),1)*0.6);
 */
 RESULT = SUUNTO_GPS_LATITUDE*10000; 
 prefix="N";
 postfix="";
}
/* next show the Longitude in DD.MM format, for lack of  jave.round fuction suunto.mod had to be used*/
else if (secs < 16) 
{
  /*
 RESULT = (Suunto.abs(SUUNTO_GPS_LONGITUDE) - Suunto.mod(Suunto.abs(SUUNTO_GPS_LONGITUDE),1) + Suunto.mod(Suunto.abs(SUUNTO_GPS_LONGITUDE),1)*0.6);
 */
  RESULT = SUUNTO_GPS_LONGITUDE*10000;
  prefix="W";
 postfix="!";
}
/* go back to Lattitude in DD.MM */
else if (secs >= 16 ) {
  
  RESULT = Suunto.abs((SUUNTO_GPS_LATITUDE - Suunto.mod(SUUNTO_GPS_LATITUDE,1) + Suunto.mod(SUUNTO_GPS_LATITUDE,1)*0.6));
 prefix="Lat";
 postfix="d.m";
 secs = 0;
}