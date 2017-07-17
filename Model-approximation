#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.000001
#define MAXT7_1 10
#define MAXT7_2 25

#define BETA 8/3
#define RHO 28
#define SIGMA 10.0

#define FD(t,u1,u2) ((-1.5*pow(u2,2)+1/pow(u1,3)-1)/u1)
#define F(t,u1,u2) (u2)
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


void RK_fourthO(double, double);
void Lorenz(double);

double DX(double, double, double, double);
double DY(double, double, double, double);
double DZ(double, double, double, double);

int main() {
	double r = 0.2, dt = 0.01;
	RK_fourthO(r, dt);
	system("pause");
	Lorenz(dt);
	system("pause");
	return 0;
}

void RK_fourthO(double r, double dt) {
	int i = 0;
	double  t = 0.0, un = 1 + r, ustar, udstar, utstar,
			un1, TOL, u2star, un2 = 0, u2dstar, u2tstar, 
			un2_1, delta;

	FILE *file = fopen("7_1.csv", "w");
	fprintf(file, "time, R\n");

	while (t < MAXT7_1) {

		dt = MIN(dt, MAXT7_1 - t);

		ustar = un + (dt*0.5)*F(t, un, un2);
		u2star = un2 + (dt*0.5)*FD(t, un, un2);

		udstar = un + dt*(0.5 * F(t + 0.5*dt, ustar, u2star));
		u2dstar = un2 + dt*(0.5 * FD(t + 0.5*dt, ustar, u2star));

		utstar = un + dt*F(t, udstar, u2dstar);
		u2tstar = un2 + dt*FD(t, udstar, u2dstar);

		un1 = un + (dt / 6)*(F(t, un, un2) + 2 * F(t + 0.5*dt, ustar, u2star) + 2 * F(t + (1 / 2)*dt, udstar, u2dstar) + F(t, utstar, u2tstar));
		un2_1 = un2 + (dt / 6)*(FD(t, un, un2) + 2 * FD(t + 0.5*dt, ustar, u2star) + 2 * FD(t + (1 / 2)*dt, udstar, u2dstar) + FD(t, utstar, u2tstar));

		TOL = fabs(un1 - un)*dt;
		delta = 0.84*(pow(EPSILON / TOL, 0.25));

		if (TOL <= EPSILON) {
			t += dt;
			un = un1;
			un2 = un2_1;
			dt *= delta;
			i++;
			printf("step = %d\tR = %.5lf\tV = %.5lf\ttime = %lf\n ", i, un1, un2_1, t);
			fprintf(file, "%lf,%lf\n", t, un1);
		}
		else {
			dt *= delta;
		}
	}
	fclose(file);
}
void Lorenz(double dt) {
	int i = 1;

	double  t = 0.0, usx, ussx, usssx, ux1, TOL;
	double  usy, ussy, usssy, uy1;
	double  usz, ussz, usssz, uz1;
	double x = 1.0, y = 1.0, z = 1.0;

	FILE *file = fopen("data2.csv", "w");
	fprintf(file, "t,x,y,z\n");

	while (t < MAXT7_2) {
		dt = MIN(dt, MAXT7_2 - t);

		usx = x + (dt*0.5)*DX(t, x, y, z);
		usy = y + (dt*0.5)*DY(t, x, y, z);
		usz = z + (dt*0.5)*DZ(t, x, y, z);

		ussx = x + dt*(0.5 * DX(t + 0.5*dt, usx, usy, usz));
		ussy = y + dt*(0.5 * DY(t + 0.5*dt, usx, usy, usz));
		ussz = z + dt*(0.5 * DZ(t + 0.5*dt, usx, usy, usz));

		usssx = x + dt*DX(t, ussx, ussy, ussz);
		usssy = y + dt*DY(t, ussx, ussy, ussz);
		usssz = z + dt*DZ(t, ussx, ussy, ussz);

		ux1 = x + (dt / 6)*(DX(t, x, y, z) + 2 * DX(t + 0.5*dt, usx, usy, usz) + 2 * DX(t + (1 / 2)*dt, ussx, ussy, ussz) + DX(t, usssx, usssy, usssz));
		uy1 = y + (dt / 6)*(DY(t, x, y, z) + 2 * DY(t + 0.5*dt, usx, usy, usz) + 2 * DY(t + (1 / 2)*dt, ussx, ussy, ussz) + DY(t, usssx, usssy, usssz));
		uz1 = z + (dt / 6)*(DZ(t, x, y, z) + 2 * DZ(t + 0.5*dt, usx, usy, usz) + 2 * DZ(t + (1 / 2)*dt, ussx, ussy, ussz) + DZ(t, usssx, usssy, usssz));


		TOL = fabs(ux1 - x)*dt;
		double	delta = 0.5*(pow(EPSILON / TOL, 0.25));

		if (TOL <= EPSILON) {
			t += dt;
			x = ux1;
			y = uy1;
			z = uz1;
			dt *= delta;
			printf("step = %d\tX = %.5lf\tY = %.5lf\tZ = %.5lf\tTime = %.5lf\t\n", i, ux1, uy1, uz1, t);
			i++;
			fprintf(file, "%lf,%lf,%lf,%lf\n", t, ux1, uy1, uz1);
		}
		else {
			dt *= delta;
		}
	}
	fclose(file);
}
double DX(double t, double x, double y, double z) {
	return (SIGMA*(y - x));
}
double DY(double t, double x, double y, double z) {
	return (x*(RHO - z)) - y;
}
double DZ(double t, double x, double y, double z) {
	return (x*y) - (BETA*z);
}

