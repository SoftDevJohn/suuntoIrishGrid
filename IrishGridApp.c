/* 
  Convert Lat/Long values to Irish Grid Reference values
  Example:
  Lugnaquilla: N 52.96712 / -6.46464 => 303207 / 191773
  Lugnaquilla: N 52.96712 / -6.46464 => 32.17 (to the nearest 100 m2 using 4 digits, 
  by losing the most significant digit of the Easting and Northing)
*/
var digits = 6;


var pi = 3.141592653589793238462643383279;
var deg2rad = pi/180;

var lat=SUUNTO_GPS_LATITUDE;
var lon=SUUNTO_GPS_LONGITUDE;
/* Locationof OSO
lat =  53+21/60+50.5441/3600;
lon = -(6+20/60+52.9181/3600);
*/
/* spot on
lat =  53+21/60+50.5441/3600;
lon = -(6+20/60+52.9181/3600);
 */
 
/*
var lat = 52.9668366960;
var lon = -6.4637116624;

var lat = 54.75408;
var lon = -8.08118;
*/

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
  var Mb = (3*n + 3*n*n + (21/8)*n3) * Suunto.sin(latRad-lat0) * Suunto.cos(latRad+lat0);
  var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(latRad-lat0)) * Math.cos(2*(latRad+lat0));
  var Md = (35/24)*n3 * Math.sin(3*(latRad-lat0)) * Math.cos(3*(latRad+lat0));

var M = b * F0 * (Ma - Mb + Mc - Md);
  var cos3lat = cosLat*cosLat*cosLat;
  var cos5lat = cos3lat*cosLat*cosLat;
  var tan2lat = Math.tan(latRad)*Math.tan(latRad);
  var tan4lat = tan2lat*tan2lat;

  var I = M + N0;
	var II = (nu/2)*sinLat*cosLat;
  var III = (nu/24)*sinLat*cos3lat*(5-tan2lat+9*eta2);
  var IIIA = (nu/720)*sinLat*cos5lat*(61-58*tan2lat+tan4lat);
  var IV = nu*cosLat;
  var V = (nu/6)*cos3lat*(nu/rho-tan2lat);
  var VI = (nu/120) * cos5lat * (5 - 18*tan2lat + tan4lat + 14*eta2 - 58*tan2lat*eta2);

  var dLon = lonRad-lon0;
  var dLon2 = dLon*dLon;
var dLon3 = dLon2*dLon;
var dLon4 = dLon3*dLon;
var dLon5 = dLon4*dLon;
var dLon6 = dLon5*dLon;

  var N = I + II*dLon2 + III*dLon4 + IIIA*dLon6;
  var E = E0 + IV*dLon + V*dLon3 + VI*dLon5;
/*
  E = Math.floor((Suunto.mod(E,10000000))/Math.pow(10,5-digits/2));
  N = Math.floor((Suunto.mod(N,10000000))/Math.pow(10,5-digits/2));
*/
  E = Math.round((Suunto.mod(E,10000000))/Math.pow(10,5-digits/2));
  N = Math.round((Suunto.mod(N,10000000))/Math.pow(10,5-digits/2));


/*  GOOD
Give a 6-figure grid reference pin-pointing 100 m2 in and 100 km.
This gives an an area that is 1/10 of the grid square.
var gridRef = Suunto.mod(E,1000)*1000+Suunto.mod(N,1000);
*/

/* 
Give a 6-figure grid reference pin-pointing 100 m2 in and 10 km2 area.
Remove the most signifcant figure, leaving two digits which are the second 
digit of the margins and the place within 1/10 of the grid square. 
*/
var gridRef = Suunto.mod(E,100)+Suunto.mod(N,100)/100;

  
RESULT = gridRef;




