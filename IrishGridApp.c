/* While in sport mode do this once per second */
/* 
  Convert Lat/Long values to RD values
  "RijksDriehoeksmeting"-coordinates
  For the X and the Y coordinate there are separate apps
  This one is for the RD-Y coordinate
  INPUT: SUUNTO_GPS_LATITUDE, SUUNTO_GPS_LONGITUDE
  RESULT: RDY
  Lugnaquilla: N 52.96712 / -6.46464 => 303207 / 191773
*/
var digits = 8;

var pi = 3.14159265358979;
var deg2rad = pi/180;

/*
var lat=SUUNTO_GPS_LATITUDE;
var lon=SUUNTO_GPS_LONGITUDE;
*/
var lat = 52.96712;
var lon = -6.46464;

var latRad = lat * deg2rad; 
var lonRad = lon * deg2rad;

var a = 6378137; 
var b = 6356752.3142;
var tx = -482.53;
var ty = 130.596;
  
var tz = -564.557;
var rx = 1.042/3600 * deg2rad;  
var ry = 0.214/3600 * deg2rad;

var rz = 0.631/3600 * deg2rad;

var s1 = -8.15/1000000 + 1;    
  
var F0 = 1.000035;                 

var lat0 = (53.5) * deg2rad;		
var lon0 = (-8) * deg2rad;  		
var N0 = 250000;
var E0 = 200000;     

var F0 = 1.000035;        
var lat0 = (53.5) * deg2rad;		
var lon0 = (-8) * deg2rad;
var N0 = 250000;
var E0 = 200000; 

var sinPhi = Suunto.sin(latRad);
var cosPhi = Suunto.cos(latRad);
var sinLambda = Suunto.sin(lonRad);
var cosLambda = Suunto.cos(lonRad);
var H = 0;

var eSq = (a*a - b*b) / (a*a);
var nu = a / Math.sqrt(1 - eSq*sinPhi*sinPhi);

var x1 = (nu+H) * cosPhi * cosLambda;
var y1 = (nu+H) * cosPhi * sinLambda;
var z1 = ((1-eSq)*nu + H) * sinPhi;
 
var x2 = tx + x1*s1 - y1*rz + z1*ry;
var y2 = ty + x1*rz + y1*s1 - z1*rx;
var z2 = tz - x1*ry + y1*rx + z1*s1;

var a = 6377340.189;
var b = 6356034.447;
var precision = 4 / a; 
 
var eSq = (a*a - b*b) / (a*a);
var p = Suunto.sqrt(x2*x2 + y2*y2);

var phi = Suunto.atan2(z2, p*(1-eSq));
  
var phiP = 2*pi;

 while (Math.abs(phi-phiP) > precision) {
    nu = a / Math.sqrt(1 - eSq*Math.sin(phi)*Math.sin(phi));
    phiP = phi;
    phi = Math.atan2(z2 + eSq*nu*Math.sin(phi), p);
  }

  var lambda = Math.atan2(y2, x2);
  H = p/Math.cos(phi) - nu;
 
	lat = phi * 180 / pi;
	lon = lambda * 180 / pi;

 var n = (a-b)/(a+b);
var n2 = n*n;
var n3 = n*n*n;

  var cosLat = Math.cos(latRad);
var sinLat = Math.sin(latRad);
  nu = a*F0/Math.sqrt(1-eSq*sinLat*sinLat);              
  var rho = a*F0*(1-eSq)/Math.pow(1-eSq*sinLat*sinLat, 1.5);  
  var eta2 = nu/rho-1;


  var Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (latRad-lat0);
  var Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(latRad-lat0) * Math.cos(latRad+lat0);
  var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(latRad-lat0)) * Math.cos(2*(latRad+lat0));
  var Md = (35/24)*n3 * Math.sin(3*(latRad-lat0)) * Math.cos(3*(latRad+lat0));
  var M = b * F0 * (Ma - Mb + Mc - Md);              
/* n 0.001673220384152102 GOOD */
/* n2 0.000002799666453942108 GOOD */
/* latRad 0.9244456947119299 GOOD */
/*lat0 0.9337511498169663 GOOD */
/*Mb  0.000013262492645875477 XX*/
/*Mc  8.212845405154541e-8 XX*/
/*Md   -1.447834787272973e-10 XX*/


  var cos3lat = cosLat*cosLat*cosLat;
  var cos5lat = cos3lat*cosLat*cosLat;
  var tan2lat = Math.tan(latRad)*Math.tan(latRad);
  var tan4lat = tan2lat*tan2lat;


/* M -59330.814905801046 */

  var I = M + N0;

	var II = (nu/2)*sinLat*cosLat;
  var III = (nu/24)*sinLat*cos3lat*(5-tan2lat+9*eta2);
  var IIIA = (nu/720)*sinLat*cos5lat*(61-58*tan2lat+tan4lat);
  var IV = nu*cosLat;
  var V = (nu/6)*cos3lat*(nu/rho-tan2lat);
  var VI = (nu/120) * cos5lat * (5 - 18*tan2lat + tan4lat + 14*eta2 - 58*tan2lat*eta2);


/* I should be 190669.185094198 GOOD*/
  var dLon = lonRad-lon0;
  var dLon2 = dLon*dLon;
var dLon3 = dLon2*dLon;
var dLon4 = dLon3*dLon;
var dLon5 = dLon4*dLon;
var dLon6 = dLon5*dLon;

  var N = I + II*dLon2 + III*dLon4 + IIIA*dLon6;
  var E = E0 + IV*dLon + V*dLon3 + VI*dLon5;

  
/*RESULT = Suunto.mod(N/10,10000);*/
/* RESULT = Suunto.mod(E/10,10000); */

  E = Math.floor((Suunto.mod(E,10000000))/Math.pow(10,5-digits/2));
  N = Math.floor((Suunto.mod(N,10000000))/Math.pow(10,5-digits/2));




/* pick the 3 most significant digits */
E = Suunto.mod(E/1000,10000);
N = Suunto.mod(N/100,100);

var gridRef = E * 1000 ;
/*
if(E < 100){
	Esz = "0";
}
*/
RESULT=N;



/*
  var es = E.toString();
   

  for (var i=0; i<(digits/2)-es.length; i++) es = '0' + es;

  var ns = N.toString();
  for (var i=0; i<(digits/2)-ns.length; i++) ns = '0' + ns;

  var gridRef = es + ns;
*/
/*
secs = secs + 1;
if (SUUNTO_GPS_STATE < 100){
  RESULT = SUUNTO_GPS_STATE;
    prefix="GPS";
    postfix="%";
}else if (secs <= 4){
 RESULT = E; 
 prefix="E";
 postfix="";
}else if (secs > 4 && secs <= 8){
  RESULT = N;
  prefix="N";
 postfix="";
}else if (secs > 8 ) {
 
  RESULT = N;
 prefix="N";
 postfix="";
  
 secs = 0;
}
*/



/*John Costigan code End*/
/*
e = 0.081696831222;   
m = 0.003773953832;  
n = 1.00047585668;   
k = 0.9999079;      
r = 6382644.571;
pi= 3.141592653589793;

lambda0 = 5.387638889 * deg2rad;
y0 = 463000.0;

phi=SUUNTO_GPS_LATITUDE * deg2rad;
lambda=SUUNTO_GPS_LONGITUDE * deg2rad;

tsin = Suunto.sin(phi);
q = 0.5*Suunto.log((1+tsin)/(1-tsin)) - (e * 0.5*Suunto.log((1+(e*tsin))/(1-(e*tsin))));
w = n * q + m;
b = 2*Suunto.atan2(Suunto.exp(w),1) - pi/2;  
b0= 52.121097249 * deg2rad;

dl = n * (lambda - lambda0);
d = 1 + Suunto.sin(b) * Suunto.sin(b0) + Suunto.cos(b) * Suunto.cos(b0) * Suunto.cos(dl);
rdy = 2 * k * r * ((Suunto.sin(b) * Suunto.cos(b0) - Suunto.cos(b) * Suunto.sin(b0) * Suunto.cos(dl)) / d) + y0;

RESULT=rdy;
  */
