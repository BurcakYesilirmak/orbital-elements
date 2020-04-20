// data comes from https://ssd.jpl.nasa.gov/horizons.cgi#results
// Object: 02 Pallas 
// Date 2448504.500000000 JD. 1991-SEP-05  00:00:00.0000 TDB
// Obtaining Orbital Elements from Position and Velocity Vectors

#include<stdio.h>
#include<math.h>
#define PI 3.1415926535897932 
int main()
{
double k3=0.01720209895,k2,x,y,z,i,j,k,rc,x2,y2;
double dot1,dot2,norm_r,norm_v,norm_h,norm_N,cosE,sinE,tau;
double eccentricity,a,inclination,omega,w,E,M,n;
double r_vec[]= {-2.139572741328206E+00,-1.166561515936117E+00,9.803478519193857E-01};    //  r_vec[3]={X,Y,Z}  V_vec[3]={Vx,Vy,Vz} ;
double v_vec[]= {4.069295114303854E-03,-8.446241616177807E-03,5.506625366696030E-03};
double t=2448504.500000000,Tp=2448203.487344466150; 
double EC= 2.348199771398780E-01 ;
double A = 2.770202581754357E+00;
double IN= 3.481273005219021E+01;
double OM= 1.733078069342644E+02;
double W = 3.096925819937945E+02;
double MA= 6.434599701373317E+01;
k2=k3*k3;
   norm_r = sqrt(r_vec[0]* r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]);
   norm_v = sqrt(v_vec[0] * v_vec[0] + v_vec[1] * v_vec[1] + v_vec[2] * v_vec[2]);
   dot1 = (r_vec[0] * v_vec[0]) + (r_vec[1] * v_vec[1]) + (r_vec[2] * v_vec[2]);
   dot2 = (v_vec[0] * v_vec[0]) + (v_vec[1] * v_vec[1]) + (v_vec[2] * v_vec[2]);
   x = r_vec[0]*( dot2/k2 - (1.0/norm_r)) - v_vec[0]* (dot1/k2); 
   y = r_vec[1]*( dot2/k2 - (1.0/norm_r)) - v_vec[1]* (dot1/k2);
   z = r_vec[2]*( dot2/k2 - (1.0/norm_r)) - v_vec[2]* (dot1/k2);  
   eccentricity = sqrt (x*x +y*y +z*z);
   a = 1/((2.0/norm_r) - dot2/k2);
   i = r_vec[1] * v_vec[2] - r_vec[2] * v_vec[1] ;
   j = r_vec[2] * v_vec[0] - r_vec[0] * v_vec[2] ;
   k = r_vec[0] * v_vec[1] - r_vec[1] * v_vec[0];
   norm_h = sqrt(i*i+j*j+k*k);
   inclination = acos(k/norm_h)*180 /PI;
   norm_N = sqrt(j*j+i*i);
   omega = acos(-j/norm_N)*180/PI;
   w = 360 - acos((-j*x+i*y)/(eccentricity*norm_N))*180/PI;
   rc=norm_h*norm_h/k2; // rc = a*(1-eccentricity**2)
   x2 = (rc-norm_r)/eccentricity ;
   y2 = dot1/eccentricity * sqrt(rc/k2);
   cosE = (x2/a)+eccentricity ;
   sinE = y2/(a*sqrt(1-pow(eccentricity,2)));
   if (fabs(sinE) <= 0.707107)
   E = asin(fabs(sinE));
   if (fabs(cosE) <= 0.707107)
   E = acos(fabs(cosE));
   if (cosE >= 0 && sinE >=0) 
   E = acos(cosE);
   if (cosE<0 && sinE >= 0)
   E = PI-E ;
   if (cosE<0 && sinE < 0)
   E = PI+E ;
   if (cosE>=0 && sinE < 0) 
   E = 2*PI-E ;
   M = (E- eccentricity*sin(E))*180/PI ;
   n = 360/(pow(a,1.5)*365.25); // deg/day unit
   tau = t - (M*180/PI)/n ;

printf(" norm_r= %.16f norm_v= %.16f \n",norm_r,norm_v);
printf(" eccentricity = %.16f  error = %.16f \n ", eccentricity,fabs(eccentricity-EC)); 
printf(" semi major axis = %.16f error = %.16f \n",a,fabs(a-A));
printf(" inclination = %.16f error = %.16f \n",inclination,inclination-IN);
printf(" longitude of ascending node = %.16f error = %.16f \n",omega,OM-omega);
printf(" argument of perifocus = %.16f error = %.16f \n",w,W-w);
printf(" rc error control = %.16f \n", rc-a*(1-pow(eccentricity,2)));
printf(" x = %.16f AU \n ",x2);
printf(" y = %.16f AU \n ",y2);
printf(" r error = %.16f \n ",sqrt(x2*x2+y2*y2)-norm_r);
printf(" cosE = %.16f rad sinE= %.16f rad \n",cosE,sinE);
printf(" E = %.16f \n",E);
printf(" M = %.16f error = %.16f \n",M,MA-M);
printf(" Tp= %.16f error = %.12f \n",tau,Tp-tau);
    return (0);
}
