#include<stdio.h>
#include<math.h>
#define g 9.8
#define PI 3.1415926535
#define k 0.38
/* Note that height is z variable */
double Force(double z,double v)
{
	return -k*v;
}
double Fn1(double t,double x,double v)/*x',y'=*/
{
	return v;
}
double Fn2(double t,double x,double v,double z)/*x"=*/
{
	return Force(z,v);
}
double Fn3(double t,double y,double v,double z)/*y"=*/
{
	return Force(z,v);
}
double Fn4(double t,double z,double v)/*z"=*/
{
	return -g+Force(z,v);
}
/* Spherical coordinate*/
double Projectile3D(double Fn1(double,double,double),double Fn2(double,double,double,double),double Fn3(double,double,double,double)
,double Fn4(double,double,double),double t0,double x0,double y0,double z0,double v0,double dt,double angle,double phi)
{
	FILE *fptw;
	double ti=t0,xi=x0,yi=y0,zi=z0;
	double vxi=v0*cos(angle*PI/180)*sin(phi*PI/180),vyi=v0*sin(angle*PI/180)*sin(phi*PI/180),vzi=v0*cos(phi*PI/180);
	double k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,p1,p2,p3,p4,q1,q2,q3,q4;
	int i=0;
	
	fptw=fopen("Projectile.txt","wt");
	fprintf(fptw,"  ti\txi\tyi\n");

	printf("\n\n\t   Projectile Motion\n-----------RK4 (%4.2lf , %4.2lf , %4.2lf, %4.2lf , %4.2lf , %4.2lf)\n\n",x0,y0,z0,angle,phi);
	printf(" n\t    ti\t  	    xi\t  	    yi\t	   zi\n");
	printf(" %d\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\n",i++,ti,xi,yi,zi);
	do{
		fprintf(fptw,"%f\t%f\t%f\n",xi,yi,zi);
		
		k1=Fn1(ti,xi,vxi);
		l1=Fn2(ti,xi,vxi,zi);
	 	k2=Fn1(ti+dt/2,xi+k1*dt/2,vxi+l1*dt/2);
	 	l2=Fn2(ti+dt/2,xi+k1*dt/2,vxi+l1*dt/2,zi);
	  	k3=Fn1(ti+dt/2,xi+k2*dt/2,vxi+l2*dt/2);
	  	l3=Fn2(ti+dt/2,xi+k2*dt/2,vxi+l2*dt/2,zi);
		k4=Fn1(ti+dt,xi+k3*dt,vxi+l3*dt);
		l4=Fn2(ti+dt,xi+k3*dt,vxi+l3*dt,zi);
		
		m1=Fn1(ti,yi,vyi);
		n1=Fn3(ti,yi,vyi,zi);
		m2=Fn1(ti+dt/2,yi+m1*dt/2,vyi+n1*dt/2);
	 	n2=Fn3(ti+dt/2,yi+m1*dt/2,vyi+n1*dt/2,zi);
	 	m3=Fn1(ti+dt/2,yi+m2*dt/2,vyi+n2*dt/2);
	  	n3=Fn3(ti+dt/2,yi+m2*dt/2,vyi+n2*dt/2,zi);
	  	m4=Fn1(ti+dt,yi+m3*dt,vyi+n3*dt);
		n4=Fn3(ti+dt,yi+m3*dt,vyi+n3*dt,zi);
		
		p1=Fn1(ti,zi,vzi);
		q1=Fn4(ti,zi,vzi);
		p2=Fn1(ti+dt/2,zi+p1*dt/2,vzi+q1*dt/2);
	 	q2=Fn4(ti+dt/2,zi+p1*dt/2,vzi+q1*dt/2);
	 	p3=Fn1(ti+dt/2,zi+p2*dt/2,vzi+q2*dt/2);
	  	q3=Fn4(ti+dt/2,zi+p2*dt/2,vzi+q2*dt/2);
	  	p4=Fn1(ti+dt,zi+p3*dt,vzi+q3*dt);
		q4=Fn4(ti+dt,zi+p3*dt,vzi+q3*dt);
		 
		xi+=(k1+2.*k2+2.*k3+k4)*dt/6;
	 	yi+=(m1+2.*m2+2.*m3+m4)*dt/6;
	 	zi+=(p1+2.*p2+2.*p3+p4)*dt/6;
	 	vxi+=(l1+2.*l2+2.*l3+l4)*dt/6;
	 	vyi+=(n1+2.*n2+2.*n3+n4)*dt/6;
	 	vzi+=(q1+2.*q2+2.*q3+q4)*dt/6;
	 	ti+=dt;
	 	printf(" %d\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\n",i++,ti,xi,yi,zi);
	}while( zi>0);
	fclose(fptw);
	return sqrt(pow(xi,2.)+pow(yi,2.));
}	

int main(void)
{
	printf("--------- Distance max = %f ",Projectile3D(Fn1,Fn2,Fn3,Fn4,0,0,0,1.4,35.5,0.01,45,60));
	return 1;
}

