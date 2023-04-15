#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.141592653589793
#define N_t 1000				//Number of trajectories
#define N 100					//Number of layers
#define lambda 5				//Mean free path (unspecified)
#define P_abs 0.1				//Absorption probability
#define eps 6.0				//Screening coefficient
#define OPTION 1

/**
 * #USAGE:
 *		(OPTION)
 *		(1) Conventional sampling;
 *		(2) "M" sampling;
 *		(3) Conventional sampling (corrected);
 * #NOTE: the following piece of code is intended for didactical purpose only, in this sense, the layer width D is taken to be 1.
 *		Further implementation is required in order to describe inelastic collisions, three-dimensional trajectories and much more.
 *
 * #GiulioLoBello
 */


int ABSORBED;
int INELASTIC;
int BACKWARD;
int BACKSCATTERED;
int NEWLAYER;

//step-length, interaction probability, horizontal displacement, depth, direction angle, cosine value of new angular displacement
double s, I, x, z, theta, cth_new;

//random number in range (0,1)
double rnd()
{
	return (double) rand()/RAND_MAX;
}

//angular conversion tools
double toRad(double x)
{
	return x*pi/180;
}
double toDeg(double x)
{
	return x*180/pi;
}

//initial conditions
void initialize()
{
	ABSORBED = 0;
	INELASTIC = 0;
	BACKSCATTERED = 0;
	NEWLAYER = 1;
	s = 0;
	theta = 0;
	x = 0;
	z = 0;
}

//direction switch
void invertDirection()
{
	if (!BACKWARD)
		BACKWARD = 1;
	else 
		BACKWARD = 0; 
}

//sets all elements to zero
void zeroArray(int *a, int lengthOfArray) 
{
	int i;
	for (i = 0; i < lengthOfArray; i++)
		a[i] = 0;
}

void interaction()
{
	double R, cth;
	R = rnd();
	if (R < P_abs)
		ABSORBED = 1;
	else {
		R = rnd();
		cth_new = 1 + (R - 1)/(eps*R + 0.5);
		if (rnd() > 0.5) {
			theta += acos(cth_new);
			if (theta > pi)
				theta -= 2*pi;
		} else {
			theta -= acos(cth_new);
			if (theta < -pi)
				theta += 2*pi;
		}
	}
}

//particle step
void singleStep() 
{
	double R, delta;
	
	if (OPTION == 1) {
		s = -lambda*log(1 - rnd());
	} else if (OPTION == 2) {
		s = 1;
	} else if (OPTION == 3) {
		s = -(1./(1 - exp(-1./lambda)))*log(1 - rnd());
	}
	I = 1 - exp(-s/lambda);
	
	x += s*sin(theta);
	z += s*cos(theta);
	
	if (z < 0) {
		BACKSCATTERED = 1;
		z = 0;
	}
	if (z > N - 1)
		z = N - 1;
		
	R = rnd();
	if (R < I) {
		interaction();
		if (INELASTIC) {
		}
	}
}

/******************************************************************************************************/

int main()
{
	int i, j, step;
	int absorbed[N], incident[N];
	FILE *out, *trk;
	srand(215255);
	
	if (OPTION == 1)
		out = fopen("C.dat", "w+");
	else if (OPTION == 2)
		out = fopen("M.dat", "w+");
	else if (OPTION == 3)
		out = fopen("C_corr.dat", "w+");
	else {
		printf("Invalid (OPTION)\n");
		return 1;
	}	
	trk = fopen("track.dat", "w+");
	
	zeroArray(absorbed, N);
	zeroArray(incident, N);
	
	for (i = 0; i < N_t; i++) {
		initialize();
		step = 0;
		fprintf(trk, "%d	%d\n", 0, 0);
		while (!ABSORBED && !BACKSCATTERED) {
			if ((int) z != j) {
				j = z;
				NEWLAYER = 1;
			}
			if (NEWLAYER) {
				incident[j]++;
				NEWLAYER = 0;
			}
			singleStep();
			step++;
			fprintf(trk, "%lg	%lg\n", x, z);	
			if (ABSORBED)
				absorbed[j]++;
			if (z == N - 1)
				break;			//END OF THE ROAD		
		}
		fprintf(trk, "\n");
	}
	
	for (j = 0; j < N; j++)
		fprintf(out,"%d	%lg\n", j + 1, (double) absorbed[j]/N_t);
	
	fclose(out);
	fclose(trk);
	return 0;
}



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	




