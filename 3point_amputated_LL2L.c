// Usage: ./"program_name" <fixed degree> <field> <random seed>
// N_min and N_max are the numbers that identify the first and last data respectively
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10000000

#define N_BP 1000000

#define N_min 300
#define N_max 401

#define MAXRAND (4294967296ULL)
#define FNORM (2.3283064365e-10)
#define RANDOM ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM) 
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)
#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)
#define sign(x) (x > 0 ? 1 : (x < 0 ? -1 : 0))

double u[N_BP], new_u[N_BP];



struct Juu_type 
{ 
	double J;
	double u0;
	double uL;
};
typedef struct Juu_type Juu_type;

struct ww_type
{
	double u0;
	double uL;
};
typedef struct ww_type ww_type;

struct params
{
	double probJ;	 
	double prob_con; 
	int n_con;	
	int n_dis;	
};
typedef struct params params;

int z;
double field; 

Juu_type *new_Juu; 
ww_type *new_ww_con, *new_ww_dis;

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

unsigned rand4init(void)
{
	unsigned long long y;

	y = myrand * 16807LL;
	myrand = (y & 0x7fffffff) + (y >> 31);
	if (myrand & 0x80000000)
		myrand = (myrand & 0x7fffffff) + 1;
	return myrand;
}

void Init_Random(void)
{
	int i;

	ip = 128;
	ip1 = ip - 24;
	ip2 = ip - 55;
	ip3 = ip - 61;

	for (i = ip3; i < ip; i++)
		ira[i] = rand4init();
}

double gaussRan(double sigma)
{
	double fac, rsq, v1, v2;

	do
	{
		v1 = 2.0 * FRANDOM - 1.0;
		v2 = 2.0 * FRANDOM - 1.0;
		rsq = v1 * v1 + v2 * v2;
	} while (rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0 * log(rsq) / rsq);
	return v1 * fac * sigma;
}

int theta(double a, double b)
{ 

	return (a - b > 0 ? 1 : 0);
}

void oneStep_BP(int degree, double field)
{ 

	int i, j, k, ran;
	double sum, J, sigmaJ;

	sigmaJ = 1. / sqrt((double)z - 1.);
	k = degree - 1;
	for (i = 0; i < N_BP; i++)
	{
		sum = field;
		for (j = 0; j < k; j++)
		{
			ran = (int)(FRANDOM * N_BP);
			sum += u[ran];
		}
		J = gaussRan(sigmaJ);
		new_u[i] = (fabs(sum) < fabs(J) ? sum * (double)sign(J) : (double)sign(sum) * J);
	}
	for (i = 0; i < N_BP; i++)
	{
		u[i] = new_u[i];
	}
}

void Init_triplet(Juu_type *Juu, params *p)
{

	int i;
	double sigmaJ;

	sigmaJ = 1. / sqrt((double)z - 1.);

	for (i = 0; i < N; i++)
	{
		//    Juu[i].J=pm1;
		Juu[i].J = gaussRan(sigmaJ);
		Juu[i].u0 = 0;
		Juu[i].uL = 0;
	}
	p->n_con = 0;
	p->n_dis = 0;
	p->probJ = 1;
	p->prob_con = 0;
}

void new_triplet_or_couple(double J1, double J2, double u0, double h, int *n_triplet, int *n_couple_con)
{
	double A, B, C, D;

	A = fabs(J1 + J2 + h);
	B = fabs(J1 + J2 - h);
	C = fabs(J1 - J2 + h);
	D = fabs(-J1 + J2 + h);
	if (fabs(h) > fabs(J1) + fabs(J2))
	{
		new_ww_con[*n_couple_con].u0 = u0 + (A + C - B - D) / 4.;
		new_ww_con[*n_couple_con].uL = (A + D - C - B) / 4.; 
		(*n_couple_con)++;
	}
	else
	{
		new_Juu[*n_triplet].u0 = u0 + (A + C - B - D) / 4.;
		new_Juu[*n_triplet].uL = (A + D - C - B) / 4.;
		new_Juu[*n_triplet].J = (A + B - C - D) / 4.;
		(*n_triplet)++;
	}
}

void just_triplet(double J1, double J2, double u0, double h, int *n_triplet)
{
	double A, B, C, D;

	if (fabs(h) < fabs(J1) + fabs(J2))
	{
		A = fabs(J1 + J2 + h);
		B = fabs(J1 + J2 - h);
		C = fabs(J1 - J2 + h);
		D = fabs(-J1 + J2 + h);
		new_Juu[*n_triplet].u0 = u0 + (A + C - B - D) / 4.;
		new_Juu[*n_triplet].uL = (A + D - C - B) / 4.;
		new_Juu[*n_triplet].J = (A + B - C - D) / 4.;
		(*n_triplet)++;
	}
}

void new_couple_con_or_dis(double J1, double J2, double u0, double h, int *i, int *j, int *prob_con)
{

	double A, B, C, D;

	A = fabs(J1 + J2 + h);
	B = fabs(J1 + J2 - h);
	C = fabs(J1 - J2 + h);
	D = fabs(-J1 + J2 + h);
	if (fabs(h) < fabs(J1) + fabs(J2))
	{
		new_ww_con[*i].u0 = u0 + (A + C - B - D) / 4.;
		new_ww_con[*i].uL = (A + D - C - B) / 4.;
		(*i)++;
		(*prob_con)++;
	}
	else
	{
		new_ww_dis[*j].u0 = u0 + (A + C - B - D) / 4.;
		new_ww_dis[*j].uL = (A + D - C - B) / 4.;
		(*j)++;
	}
}

void just_couple_con(double J1, double J2, double u0, double h, int *i, int *prob_con)
{

	double A, B, C, D;
	if (fabs(h) < fabs(J1) + fabs(J2))
	{
		A = fabs(J1 + J2 + h);
		B = fabs(J1 + J2 - h);
		C = fabs(J1 - J2 + h);
		D = fabs(-J1 + J2 + h);
		new_ww_con[*i].u0 = u0 + (A + C - B - D) / 4.;
		new_ww_con[*i].uL = (A + D - C - B) / 4.;
		(*i)++;
		(*prob_con)++;
	}
}

void new_couple_dis(double J1, double J2, double u0, double h, int *i)
{

	double A, B, C, D;

	A = fabs(J1 + J2 + h);
	B = fabs(J1 + J2 - h);
	C = fabs(J1 - J2 + h);
	D = fabs(-J1 + J2 + h);
	new_ww_dis[*i].u0 = u0 + (A + C - B - D) / 4.;
	new_ww_dis[*i].uL = (A + D - C - B) / 4.;
	(*i)++;
}

double Generate_field(double u2, double u3)
{

	double h;
	int i;

	h = u2 + u3 + field;
	//  h=u2+u3+field*pm1;
	for (i = 0; i < z - 2; i++)
	{
		h += u[(int)(FRANDOM * N_BP)];
	}

	return h;
}

int extract_random(Juu_type *Juu, ww_type *ww_con, ww_type *ww_dis, double *u1, double *u2, double *J1, double *J2, params *p)
{ 
	int chosen, check_triplet = 2;
	double ran;

	ran = FRANDOM;
	if (ran < p->probJ)
	{
		chosen = (int)(FRANDOM * N);

		

		*J1 = Juu[chosen].J;
		*u1 = Juu[chosen].u0;
		*u2 = Juu[chosen].uL;
	}
	else
	{
		if (ran < p->probJ + p->prob_con)
		{
			chosen = (int)(FRANDOM * p->n_con); 
			*J1 = 0;
			*u1 = ww_con[chosen].u0;
			*u2 = ww_con[chosen].uL;
			check_triplet--;
		}
		else
		{
			chosen = (int)(FRANDOM * p->n_dis); 
			*J1 = 0;
			*u1 = ww_dis[chosen].u0;
			*u2 = ww_dis[chosen].uL;
			check_triplet -= 2;
		}
	}

	*J2 = gaussRan(1. / sqrt((double)z - 1.));

	return check_triplet;
}

void extract_triplet(Juu_type *Juu, double *u1, double *u2, double *J1, double *J2)
{
	int chosen;

	chosen = (int)(FRANDOM * N);
	*J1 = Juu[chosen].J;
	*u1 = Juu[chosen].u0;
	*u2 = Juu[chosen].uL;
	*J2 = gaussRan(1. / sqrt((double)z - 1.));
}

void copy(Juu_type *Juu, ww_type *ww_con, ww_type *ww_dis)
{
	int i;

	for (i = 0; i < N; i++)
	{
		ww_con[i] = new_ww_con[i];
		ww_dis[i] = new_ww_dis[i];
		Juu[i] = new_Juu[i];
	}
}
int extract_triplet_or_con(Juu_type *Juu, ww_type *ww_con, double *u1, double *u2, double *J1, double *J2, params *p)
{
	int chosen, check_triplet = 2;
	double ran, prob_extr_J;

	prob_extr_J = p->probJ / (p->probJ + p->prob_con);

	ran = FRANDOM;
	if (ran < prob_extr_J)
	{
		chosen = (int)(FRANDOM * N);
		*J1 = Juu[chosen].J;
		*u1 = Juu[chosen].u0;
		*u2 = Juu[chosen].uL;
	}
	else
	{
		chosen = (int)(FRANDOM * p->n_con); 
		*J1 = 0;
		*u1 = ww_con[chosen].u0;
		*u2 = ww_con[chosen].uL;
		check_triplet--;
	}
	*J2 = gaussRan(1. / sqrt((double)z - 1.));

	return check_triplet;
}



void Lplus1(Juu_type *Juu, ww_type *ww_con, ww_type *ww_dis, params *p)
{

	int norm_prob = 0, norm_prob_con = 0;
	int check_triplet, n_triplet = 0, n_couple_con = 0, n_couple_dis = 0, prob_BinB = 0;
	double J1, J2, u1, u2, h, probJ_old, prob_con_old;

	while (n_triplet < N && n_couple_con < N && n_couple_dis < N)
	{

		check_triplet = extract_random(Juu, ww_con, ww_dis, &u1, &u2, &J1, &J2, p); 

		h = Generate_field(u2, 0);

		if (check_triplet == 2)
		{
			new_triplet_or_couple(J1, J2, u1, h, &n_triplet, &n_couple_con);
			norm_prob++;
		}
		else
		{
			if (check_triplet == 1)
			{
				new_couple_con_or_dis(J1, J2, u1, h, &n_couple_con, &n_couple_dis, &prob_BinB);
				norm_prob_con++;
			}
			else
			{
				new_couple_dis(J1, J2, u1, h, &n_couple_dis);
			}
		}
		
	}

	while (n_triplet < N && n_couple_con < N)
	{
		check_triplet = extract_triplet_or_con(Juu, ww_con, &u1, &u2, &J1, &J2, p);
		h = Generate_field(u2, 0);
		if (check_triplet == 2)
		{
			new_triplet_or_couple(J1, J2, u1, h, &n_triplet, &n_couple_con);
			norm_prob++;
		}
		else
		{
			if (check_triplet == 1)
			{
				just_couple_con(J1, J2, u1, h, &n_couple_con, &prob_BinB);
				norm_prob_con++;
			}
		}
	}
	while (n_triplet < N)
	{
		extract_triplet(Juu, &u1, &u2, &J1, &J2);
		h = Generate_field(u2, 0);
		just_triplet(J1, J2, u1, h, &n_triplet);
		norm_prob++;
	}
	
	copy(Juu, ww_con, ww_dis);

	probJ_old = p->probJ;
	prob_con_old = p->prob_con;
	p->probJ = probJ_old * (double)N / (double)norm_prob; 
	p->prob_con = probJ_old * (1 - (double)N / (double)norm_prob) + (norm_prob_con != 0 ? prob_con_old * (double)prob_BinB / (double)norm_prob_con : 0); 

	p->n_con = n_couple_con;
	p->n_dis = n_couple_dis;
}

void KahanSum_mean_error(double *input, int n, double *mean, double *error)
{

	double sum = 0.0, y, t;
	double c = 0.0; // A running compensation for lost low-order bits.
	int i;

	for (i = 0; i < n; i++)
	{
		y = input[i] - c;  // So far, so good: c is zero.
		t = sum + y;	   // Alas, sum is big, y small, so low-order digits of y are lost.
		c = (t - sum) - y; // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
		sum = t;	   // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
	}			   // Next time around, the lost low part will be added to y in a fresh attempt.

	*mean = sum / (double)n;
	sum = 0.;
	c = 0.;

	for (i = 0; i < n; i++)
	{
		y = input[i] * input[i] - c;
		t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}

	*error = sqrt((sum / ((double)n) - (*mean) * (*mean)) / (n - 1.));
}

void correlations(Juu_type *Juu, double *con, double *error, double probJ)
{ 

	int i, j;
	double h1, h2, J;

	(*con) = 0;

	for (i = 0; i < N; i++)
	{
		h1 = Juu[i].u0 + field;
		//  h=u2+u3+field*pm1;
		for (j = 0; j < z - 1; j++)
		{
			h1 += u[(int)(FRANDOM * N_BP)];
		}
		h2 = Juu[i].uL + field;
		//  h=u2+u3+field*pm1;
		for (j = 0; j < z - 1; j++)
		{
			h2 += u[(int)(FRANDOM * N_BP)];
		}
		J = Juu[i].J;
		(*con) += (double)theta(fabs(J), fabs(h1));
	}

	// KahanSum_mean_error(data,N,con,error);
	(*con) /= (double)N;
	(*error) = sqrt(((*con) - (*con) * (*con)) / (N - 1.));

	(*con) *= probJ;
	(*error) *= probJ;



}

double Generate_field_loop(double u1, double u2)
{ 
	double h;
	int i;

	h = u1 + u2 + field;
	//  h=u2+u3+field*pm1;
	for (i = 0; i < z - 3; i++)
	{
		h += u[(int)(FRANDOM * N_BP)];
	}

	return h;
}



int compose_link(double J1, double J2, double u0, double uL, double h, Juu_type *Juu_temp)
{ 
	double A, B, C, D;
	int check_triplet;

	A = fabs(J1 + J2 + h);
	B = fabs(J1 + J2 - h);
	C = fabs(J1 - J2 + h);
	D = fabs(-J1 + J2 + h);
	if (fabs(h) > fabs(J1) + fabs(J2))
	{
		Juu_temp->u0 = u0 + (A + C - B - D) / 4.;
		Juu_temp->uL = uL + (A + D - C - B) / 4.;
		Juu_temp->J = 0.;
		;
		if (J1 != 0 && J2 != 0)
		{
			check_triplet = 1;
		}
		else
		{
			check_triplet = 0;
		}
	}
	else
	{
		Juu_temp->u0 = u0 + (A + C - B - D) / 4.;
		Juu_temp->uL = uL + (A + D - C - B) / 4.;
		Juu_temp->J = (A + B - C - D) / 4.;
		if (J1 != 0 && J2 != 0)
		{
			check_triplet = 2;
		}
		else
		{
			check_triplet = 1;
		}
	}
	return check_triplet;
}

int loop_correlation(double J2, double J3, double ux2, double ux3, double uy2, double uy3, int compute_l2, int compute_l3)
{

	double h0;
	Juu_type Jloop;
	int corr, j;

	Jloop.J = J2 + J3;
	h0 = field + ux2 + ux3;
	for (j = 0; j < z - 2; j++)
	{
		h0 += u[(int)(FRANDOM * N_BP)];
	}
	corr = theta(fabs(Jloop.J), fabs(h0));

	



	if (compute_l2 != 0) 
	{
		
		h0 = ux2 + field;
		for (j = 0; j < z - 1; j++)
		{
			h0 += u[(int)(FRANDOM * N_BP)];
		}
		corr -= theta(fabs(J2), fabs(h0));
	}

	if (compute_l3 != 0)
	{

		
		h0 = ux3 + field;
		for (j = 0; j < z - 1; j++)
		{
			h0 += u[(int)(FRANDOM * N_BP)];
		}
		corr -= theta(fabs(J3), fabs(h0));
	}

	return corr;
}

void loop(int L, Juu_type *JuuL, ww_type *ww_conL, ww_type *ww_disL, Juu_type *Juu2L, ww_type *ww_con2L, ww_type *ww_dis2L, params paramsL, params params2L, FILE *f)
{
	int ran, corrAA = 0, corr_quadAA = 0, corrAB_1 = 0, corr_quadAB_1 = 0, corrAB_2 = 0, corr_quadAB_2 = 0, corrAC_1 = 0, corr_quadAC_1 = 0, corrAC_2 = 0, corr_quadAC_2 = 0, temp, n;
	double J2, J3, u0, uL, ux2, ux3, uy2, uy3, corr_tot, error_corrAA, error_corrAB, error_corrAC, prob_disL, prob_dis2L, corrAAdouble, corrAB1double, corrAB2double, corrAC1double, corrAC2double;

	for (n = 0; n < (double)N / 4.; n++) 
	{
		
		ran = (int)(FRANDOM * N);
		J2 = JuuL[ran].J;
		ux2 = JuuL[ran].u0;
		uy2 = JuuL[ran].uL;

		ran = (int)(FRANDOM * params2L.n_con);
		ux3 = ww_con2L[ran].u0;
		uy3 = ww_con2L[ran].uL;
		J3 = 0.;

		
		temp = loop_correlation(J2, J3, ux2, ux3, uy2, uy3, 1, 0);

		corrAB_1 += temp;
		corr_quadAB_1 += temp * temp;
	}

	corrAB1double = (double)corrAB_1 * paramsL.probJ * params2L.prob_con / ((double)n); 
	corr_tot = corrAB1double;
	error_corrAB = paramsL.probJ * params2L.prob_con *
		       sqrt(((double)corr_quadAB_1 / ((double)n) - (double)corrAB_1 / ((double)n) * (double)corrAB_1 / ((double)n)) / ((double)n - 1.));

	for (n = 0; n < (double)N / 4.; n++)
	{
		
		ran = (int)(FRANDOM * N);
		J2 = Juu2L[ran].J;
		ux2 = Juu2L[ran].u0;
		uy2 = Juu2L[ran].uL;

		ran = (int)(FRANDOM * paramsL.n_con);
		ux3 = ww_conL[ran].u0;
		uy3 = ww_conL[ran].uL;
		J3 = 0.;

		
		temp = loop_correlation(J2, J3, ux2, ux3, uy2, uy3, 1, 0);

		corrAB_2 += temp;
		corr_quadAB_2 += temp * temp;
	}

	corrAB2double = (double)corrAB_2 * params2L.probJ * paramsL.prob_con / ((double)n);

	corr_tot += (double)corrAB_2 * params2L.probJ * paramsL.prob_con / ((double)n);
	error_corrAB += params2L.probJ * paramsL.prob_con *
			sqrt(((double)corr_quadAB_2 / ((double)n) - (double)corrAB_2 / ((double)n) * (double)corrAB_2 / ((double)n)) / ((double)n - 1.));

	for (n = 0; n < (double)N / 4.; n++)
	{
		
		ran = (int)(FRANDOM * N);
		J2 = JuuL[ran].J;
		ux2 = JuuL[ran].u0;
		uy2 = JuuL[ran].uL;

		ran = (int)(FRANDOM * N);
		ux3 = Juu2L[ran].u0;
		uy3 = Juu2L[ran].uL;
		J3 = Juu2L[ran].J;

		temp = loop_correlation(J2, J3, ux2, ux3, uy2, uy3, 1, 1);
		corrAA += temp;
		corr_quadAA += temp * temp;
	}

	corrAAdouble = (double)corrAA * paramsL.probJ * params2L.probJ / ((double)n);
	corr_tot += corrAAdouble;
	error_corrAA = paramsL.probJ * params2L.probJ *
		       sqrt(((double)corr_quadAA / ((double)n) - (double)corrAA / ((double)n) * (double)corrAA / ((double)n)) / ((double)n - 1.));

	if (paramsL.n_dis > 0)
	
	{
		for (n = 0; n < (double)N / 4.; n++)
		{
			
			ran = (int)(FRANDOM * N);
			J2 = Juu2L[ran].J;
			ux2 = Juu2L[ran].u0;
			uy2 = Juu2L[ran].uL;

			ran = (int)(FRANDOM * paramsL.n_dis);
			ux3 = ww_disL[ran].u0;
			uy3 = ww_disL[ran].uL;
			J3 = 0.;

			
			temp = loop_correlation(J2, J3, ux2, ux3, uy2, uy3, 1, 0);

			corrAC_1 += temp;
			corr_quadAC_1 += temp * temp;
		}
	}
	prob_disL = 1. - paramsL.probJ - paramsL.prob_con;
	corrAC1double = (double)corrAC_1 * params2L.probJ * prob_disL / ((double)n);
	corr_tot += corrAC1double;
	error_corrAC = params2L.probJ * prob_disL *
		       sqrt(((double)corr_quadAC_1 / ((double)n) - (double)corrAC_1 / ((double)n) * (double)corrAC_1 / ((double)n)) / ((double)n - 1.));

	if (params2L.n_dis > 0) 
	{
		for (n = 0; n < (double)N / 4.; n++)
		{
			
			ran = (int)(FRANDOM * N);
			J2 = JuuL[ran].J;
			ux2 = JuuL[ran].u0;
			uy2 = JuuL[ran].uL;

			ran = (int)(FRANDOM * params2L.n_dis);
			ux3 = ww_dis2L[ran].u0;
			uy3 = ww_dis2L[ran].uL;
			J3 = 0.;

			
			temp = loop_correlation(J2, J3, ux2, ux3, uy2, uy3, 1, 0);

			corrAC_2 += temp;
			corr_quadAC_2 += temp * temp;
		}
	}
	prob_dis2L = 1. - params2L.probJ - params2L.prob_con;
	corrAC2double = (double)corrAC_2 * paramsL.probJ * prob_dis2L / ((double)n);
	corr_tot += corrAC2double;
	error_corrAC += paramsL.probJ * prob_dis2L *
			sqrt(((double)corr_quadAC_2 / ((double)n) - (double)corrAC_2 / ((double)n) * (double)corrAC_2 / ((double)n)) / ((double)n - 1.));

	fprintf(f, "%d %g %g %g %g %g ", L, corr_tot, sqrt(error_corrAA * error_corrAA + error_corrAB * error_corrAB + error_corrAC * error_corrAC), corrAAdouble, corrAB1double + corrAB2double, corrAC1double + corrAC2double);
}

int loopvertex_correlation(double J12, double J23, double J13, double u1a, double u2a, double u2b, double u3b, double u3c, double u1c, int compute_la, int compute_lb, int compute_lc)
{

	double h1, h2;
	int corr, j;

	h1 = field + u1a + u1c;
	h2 = field + u2a + u2b;
	for (j = 0; j < z - 2; j++)
	{
		h1 += u[(int)(FRANDOM * N_BP)];
	}
	for (j = 0; j < z - 2; j++)
	{
		h2 += u[(int)(FRANDOM * N_BP)];
	}
	corr = abs( ( theta(fabs(J13),fabs(J12))*theta(fabs(J13)+fabs(J12),fabs(h1))*theta(fabs(h1),fabs(J13)-fabs(J12))*theta( fabs(J23+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))),fabs(h2+sign(J12*h1)*0.5*(fabs(h1)-fabs(J13)+fabs(J12))) )*sign(J23+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))) +  theta(fabs(J13),fabs(J12))*theta(fabs(J13),fabs(J12)+fabs(h1))*theta(fabs(J23+sign(J13)*J12),fabs(h2))*sign(J23+sign(J13)*J12) + theta(fabs(J12),fabs(J13))*theta(fabs(J13)+fabs(J12),fabs(h1))*theta(fabs(h1),fabs(J12)-fabs(J13))*theta( fabs(J23+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))),fabs(h2+sign(J12*h1)*0.5*(fabs(h1)-fabs(J13)+fabs(J12))) )*sign(J23+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))) +  theta(fabs(J12),fabs(J13))*theta(fabs(J12),fabs(J13)+fabs(h1))*theta(fabs(J23+sign(J12)*J13),fabs(h2+sign(J12)*h1))*sign(J23+sign(J12)*J13) )*( theta(fabs(J23),fabs(J12))*theta(fabs(J23)+fabs(J12),fabs(h2))*theta(fabs(h2),fabs(J23)-fabs(J12))*theta(fabs(J13+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))),fabs(h1+sign(J12*h2)*0.5*(fabs(h2)+fabs(J12)-fabs(J23))))*sign(J13+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))) + theta(fabs(J23),fabs(J12))*theta(fabs(J23),fabs(J12)+fabs(h2))*theta(fabs(J13+sign(J23)*J12),fabs(h1))*sign(J13+sign(J23)*J12) + theta(fabs(J12),fabs(J23))*theta(fabs(J23)+fabs(J12),fabs(h2))*theta(fabs(h2),fabs(J12)-fabs(J23))*theta(fabs(J13+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))),fabs(h1+sign(J12*h2)*0.5*(fabs(h2)+fabs(J12)-fabs(J23))))*sign(J13+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))) + theta(fabs(J12),fabs(J23))*theta(fabs(J12)-fabs(J23),fabs(h2))*theta(fabs(J13+sign(J12)*J23),fabs(h1+h2*sign(J12)))*sign(J13+sign(J12)*J23) ) );

	



	if (compute_la != 0 && compute_lb !=0 ) 
	{

		h1 = field + u1a;
		h2 = field + u2a + u2b;
		for (j = 0; j < z - 1; j++)
		{
			h1 += u[(int)(FRANDOM * N_BP)];
		}
		for (j = 0; j < z - 2; j++)
		{
			h2 += u[(int)(FRANDOM * N_BP)];
		}
		corr -= abs( ( theta(fabs(0.),fabs(J12))*theta(fabs(0.)+fabs(J12),fabs(h1))*theta(fabs(h1),fabs(0.)-fabs(J12))*theta( fabs(J23+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h1))),fabs(h2+sign(J12*h1)*0.5*(fabs(h1)-fabs(0.)+fabs(J12))) )*sign(J23+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h1))) +  theta(fabs(0.),fabs(J12))*theta(fabs(0.),fabs(J12)+fabs(h1))*theta(fabs(J23+sign(0.)*J12),fabs(h2))*sign(J23+sign(0.)*J12) + theta(fabs(J12),fabs(0.))*theta(fabs(0.)+fabs(J12),fabs(h1))*theta(fabs(h1),fabs(J12)-fabs(0.))*theta( fabs(J23+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h1))),fabs(h2+sign(J12*h1)*0.5*(fabs(h1)-fabs(0.)+fabs(J12))) )*sign(J23+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h1))) +  theta(fabs(J12),fabs(0.))*theta(fabs(J12),fabs(0.)+fabs(h1))*theta(fabs(J23+sign(J12)*0.),fabs(h2+sign(J12)*h1))*sign(J23+sign(J12)*0.) )*( theta(fabs(J23),fabs(J12))*theta(fabs(J23)+fabs(J12),fabs(h2))*theta(fabs(h2),fabs(J23)-fabs(J12))*theta(fabs(0.+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))),fabs(h1+sign(J12*h2)*0.5*(fabs(h2)+fabs(J12)-fabs(J23))))*sign(0.+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))) + theta(fabs(J23),fabs(J12))*theta(fabs(J23),fabs(J12)+fabs(h2))*theta(fabs(0.+sign(J23)*J12),fabs(h1))*sign(0.+sign(J23)*J12) + theta(fabs(J12),fabs(J23))*theta(fabs(J23)+fabs(J12),fabs(h2))*theta(fabs(h2),fabs(J12)-fabs(J23))*theta(fabs(0.+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))),fabs(h1+sign(J12*h2)*0.5*(fabs(h2)+fabs(J12)-fabs(J23))))*sign(0.+sign(J12*J23)*0.5*(fabs(J12)+fabs(J23)-fabs(h2))) + theta(fabs(J12),fabs(J23))*theta(fabs(J12)-fabs(J23),fabs(h2))*theta(fabs(0.+sign(J12)*J23),fabs(h1+h2*sign(J12)))*sign(0.+sign(J12)*J23) ) );

	}

	if (compute_la != 0 && compute_lc !=0 ) 
	{

		h1 = field + u1a + u1c;
		h2 = field + u2a;
		for (j = 0; j < z - 2; j++)
		{
			h1 += u[(int)(FRANDOM * N_BP)];
		}
		for (j = 0; j < z - 1; j++)
		{
			h2 += u[(int)(FRANDOM * N_BP)];
		}
		corr -= abs( ( theta(fabs(J13),fabs(J12))*theta(fabs(J13)+fabs(J12),fabs(h1))*theta(fabs(h1),fabs(J13)-fabs(J12))*theta( fabs(0.+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))),fabs(h2+sign(J12*h1)*0.5*(fabs(h1)-fabs(J13)+fabs(J12))) )*sign(0.+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))) +  theta(fabs(J13),fabs(J12))*theta(fabs(J13),fabs(J12)+fabs(h1))*theta(fabs(0.+sign(J13)*J12),fabs(h2))*sign(0.+sign(J13)*J12) + theta(fabs(J12),fabs(J13))*theta(fabs(J13)+fabs(J12),fabs(h1))*theta(fabs(h1),fabs(J12)-fabs(J13))*theta( fabs(0.+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))),fabs(h2+sign(J12*h1)*0.5*(fabs(h1)-fabs(J13)+fabs(J12))) )*sign(0.+sign(J12*J13)*0.5*(fabs(J12)+fabs(J13)-fabs(h1))) +  theta(fabs(J12),fabs(J13))*theta(fabs(J12),fabs(J13)+fabs(h1))*theta(fabs(0.+sign(J12)*J13),fabs(h2+sign(J12)*h1))*sign(0.+sign(J12)*J13) )*( theta(fabs(0.),fabs(J12))*theta(fabs(0.)+fabs(J12),fabs(h2))*theta(fabs(h2),fabs(0.)-fabs(J12))*theta(fabs(J13+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h2))),fabs(h1+sign(J12*h2)*0.5*(fabs(h2)+fabs(J12)-fabs(0.))))*sign(J13+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h2))) + theta(fabs(0.),fabs(J12))*theta(fabs(0.),fabs(J12)+fabs(h2))*theta(fabs(J13+sign(0.)*J12),fabs(h1))*sign(J13+sign(0.)*J12) + theta(fabs(J12),fabs(0.))*theta(fabs(0.)+fabs(J12),fabs(h2))*theta(fabs(h2),fabs(J12)-fabs(0.))*theta(fabs(J13+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h2))),fabs(h1+sign(J12*h2)*0.5*(fabs(h2)+fabs(J12)-fabs(0.))))*sign(J13+sign(J12*0.)*0.5*(fabs(J12)+fabs(0.)-fabs(h2))) + theta(fabs(J12),fabs(0.))*theta(fabs(J12)-fabs(0.),fabs(h2))*theta(fabs(J13+sign(J12)*0.),fabs(h1+h2*sign(J12)))*sign(J13+sign(J12)*0.) ) );

	}

	if (compute_lb != 0 && compute_lc != 0) 
	{

		h1 = field + u1c;
		h2 = field + u2b;
		for (j = 0; j < z - 1; j++)
		{
			h1 += u[(int)(FRANDOM * N_BP)];
		}
		for (j = 0; j < z - 1; j++)
		{
			h2 += u[(int)(FRANDOM * N_BP)];
		}
		corr -= abs( ( theta(fabs(J13),fabs(0.))*theta(fabs(J13)+fabs(0.),fabs(h1))*theta(fabs(h1),fabs(J13)-fabs(0.))*theta( fabs(J23+sign(0.*J13)*0.5*(fabs(0.)+fabs(J13)-fabs(h1))),fabs(h2+sign(0.*h1)*0.5*(fabs(h1)-fabs(J13)+fabs(0.))) )*sign(J23+sign(0.*J13)*0.5*(fabs(0.)+fabs(J13)-fabs(h1))) +  theta(fabs(J13),fabs(0.))*theta(fabs(J13),fabs(0.)+fabs(h1))*theta(fabs(J23+sign(J13)*0.),fabs(h2))*sign(J23+sign(J13)*0.) + theta(fabs(0.),fabs(J13))*theta(fabs(J13)+fabs(0.),fabs(h1))*theta(fabs(h1),fabs(0.)-fabs(J13))*theta( fabs(J23+sign(0.*J13)*0.5*(fabs(0.)+fabs(J13)-fabs(h1))),fabs(h2+sign(0.*h1)*0.5*(fabs(h1)-fabs(J13)+fabs(0.))) )*sign(J23+sign(0.*J13)*0.5*(fabs(0.)+fabs(J13)-fabs(h1))) +  theta(fabs(0.),fabs(J13))*theta(fabs(0.),fabs(J13)+fabs(h1))*theta(fabs(J23+sign(0.)*J13),fabs(h2+sign(0.)*h1))*sign(J23+sign(0.)*J13) )*( theta(fabs(J23),fabs(0.))*theta(fabs(J23)+fabs(0.),fabs(h2))*theta(fabs(h2),fabs(J23)-fabs(0.))*theta(fabs(J13+sign(0.*J23)*0.5*(fabs(0.)+fabs(J23)-fabs(h2))),fabs(h1+sign(0.*h2)*0.5*(fabs(h2)+fabs(0.)-fabs(J23))))*sign(J13+sign(0.*J23)*0.5*(fabs(0.)+fabs(J23)-fabs(h2))) + theta(fabs(J23),fabs(0.))*theta(fabs(J23),fabs(0.)+fabs(h2))*theta(fabs(J13+sign(J23)*0.),fabs(h1))*sign(J13+sign(J23)*0.) + theta(fabs(0.),fabs(J23))*theta(fabs(J23)+fabs(0.),fabs(h2))*theta(fabs(h2),fabs(0.)-fabs(J23))*theta(fabs(J13+sign(0.*J23)*0.5*(fabs(0.)+fabs(J23)-fabs(h2))),fabs(h1+sign(0.*h2)*0.5*(fabs(h2)+fabs(0.)-fabs(J23))))*sign(J13+sign(0.*J23)*0.5*(fabs(0.)+fabs(J23)-fabs(h2))) + theta(fabs(0.),fabs(J23))*theta(fabs(0.)-fabs(J23),fabs(h2))*theta(fabs(J13+sign(0.)*J23),fabs(h1+h2*sign(0.)))*sign(J13+sign(0.)*J23) ) );

	}

	return corr;
}

void loopvertex(int L, Juu_type *JuuLa, ww_type *ww_conLa, ww_type *ww_disLa, Juu_type *JuuLb, ww_type *ww_conLb, ww_type *ww_disLb, Juu_type *JuuLc, ww_type *ww_conLc, ww_type *ww_disLc, params paramsLa, params paramsLb, params paramsLc, FILE *f)
{
	int ran, corrAAA = 0, corr_quadAAA = 0, corrAAB = 0, corr_quadAAB = 0, corrBAA = 0, corr_quadBAA = 0, corrABA = 0, corr_quadABA = 0, corrAAC = 0, corr_quadAAC = 0, corrACA = 0, corr_quadACA = 0, corrCAA = 0, corr_quadCAA = 0, temp, n;
	double J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, corr_tot, error_corrAAA, error_corrAAB, error_corrBAA, error_corrABA, error_corrAAC, error_corrACA, error_corrCAA, prob_disLa, prob_disLb, prob_disLc, corrAAAdouble, corrAABdouble, corrBAAdouble, corrABAdouble, corrAACdouble, corrACAdouble, corrCAAdouble;

	for (n = 0; n < (double)N / 4.; n++) 
	{
		ran = (int)(FRANDOM * N);
		J12 = JuuLa[ran].J;
		u1a = JuuLa[ran].u0;
		u2a = JuuLa[ran].uL;

		ran = (int)(FRANDOM * N);
		J23 = JuuLb[ran].J;
		u2b = JuuLb[ran].u0;
		u3b = JuuLb[ran].uL;

		ran = (int)(FRANDOM * N);
		J13 = JuuLc[ran].J;
		u3c = JuuLc[ran].u0;
		u1c = JuuLc[ran].uL;

		
		temp = loopvertex_correlation(J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, 1, 1, 1);

		corrAAA += temp;
		corr_quadAAA += temp * temp;
	}

	corrAAAdouble = (double)corrAAA * paramsLa.probJ * paramsLb.probJ * paramsLc.probJ / ((double)n);
	corr_tot = corrAAAdouble;
	error_corrAAA = paramsLa.probJ * paramsLb.probJ * paramsLc.probJ *
		       sqrt(((double)corr_quadAAA / ((double)n) - (double)corrAAA / ((double)n) * (double)corrAAA / ((double)n)) / ((double)n - 1.));

	for (n = 0; n < (double)N / 4.; n++) 
	{
		ran = (int)(FRANDOM * N);
		J12 = JuuLa[ran].J;
		u1a = JuuLa[ran].u0;
		u2a = JuuLa[ran].uL;

		ran = (int)(FRANDOM * N);
		J23 = JuuLb[ran].J;
		u2b = JuuLb[ran].u0;
		u3b = JuuLb[ran].uL;

		ran = (int)(FRANDOM * paramsLc.n_con);
		J13 = 0.;
		u3c = ww_conLc[ran].u0;
		u1c = ww_conLc[ran].uL;

		
		temp = loopvertex_correlation(J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, 1, 1, 0);

		corrAAB += temp;
		corr_quadAAB += temp * temp;
	}

	corrAABdouble = (double)corrAAB * paramsLa.probJ * paramsLb.probJ * paramsLc.prob_con / ((double)n); 
	corr_tot += corrAABdouble;
	error_corrAAB = paramsLa.probJ * paramsLb.probJ * paramsLc.prob_con *
		       sqrt(((double)corr_quadAAB / ((double)n) - (double)corrAAB / ((double)n) * (double)corrAAB / ((double)n)) / ((double)n - 1.));

	for (n = 0; n < (double)N / 4.; n++) 
	{
		ran = (int)(FRANDOM * paramsLa.n_con);
		J12 = 0.;
		u1a = ww_conLa[ran].u0;
		u2a = ww_conLa[ran].uL;

		ran = (int)(FRANDOM * N);
		J23 = JuuLb[ran].J;
		u2b = JuuLb[ran].u0;
		u3b = JuuLb[ran].uL;

		ran = (int)(FRANDOM * N);
		J13 = JuuLc[ran].J;
		u3c = JuuLc[ran].u0;
		u1c = JuuLc[ran].uL;


		temp = loopvertex_correlation(J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, 0, 1, 1);

		corrBAA += temp;
		corr_quadBAA += temp * temp;
	}

	corrBAAdouble = (double)corrBAA * paramsLa.prob_con * paramsLb.probJ * paramsLc.probJ  / ((double)n); 
	corr_tot += corrBAAdouble;
	error_corrBAA = paramsLa.prob_con * paramsLb.probJ * paramsLc.probJ  *
		       sqrt(((double)corr_quadBAA / ((double)n) - (double)corrBAA / ((double)n) * (double)corrBAA / ((double)n)) / ((double)n - 1.));

	for (n = 0; n < (double)N / 4.; n++)
	{
		ran = (int)(FRANDOM * N);
		J12 = JuuLa[ran].J;
		u1a = JuuLa[ran].u0;
		u2a = JuuLa[ran].uL;

		ran = (int)(FRANDOM * paramsLb.n_con);
		J23 = 0.;
		u2b = ww_conLb[ran].u0;
		u3b = ww_conLb[ran].uL;

		ran = (int)(FRANDOM * N);
		J13 = JuuLc[ran].J;
		u3c = JuuLc[ran].u0;
		u1c = JuuLc[ran].uL;

		
		temp = loopvertex_correlation(J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, 1, 0, 1);

		corrABA += temp;
		corr_quadABA += temp * temp;
	}

	corrABAdouble = (double)corrABA * paramsLa.probJ * paramsLb.prob_con * paramsLc.probJ  / ((double)n);
	corr_tot += corrABAdouble;
	error_corrABA = paramsLa.probJ *  paramsLb.prob_con * paramsLc.probJ  *
		       sqrt(((double)corr_quadABA / ((double)n) - (double)corrABA / ((double)n) * (double)corrABA / ((double)n)) / ((double)n - 1.));

	error_corrAAB += error_corrABA + error_corrBAA;

	if (paramsLc.n_dis > 0)
	
	{
		for (n = 0; n < (double)N / 4.; n++) 
	{
		ran = (int)(FRANDOM * N);
		J12 = JuuLa[ran].J;
		u1a = JuuLa[ran].u0;
		u2a = JuuLa[ran].uL;

		ran = (int)(FRANDOM * N);
		J23 = JuuLb[ran].J;
		u2b = JuuLb[ran].u0;
		u3b = JuuLb[ran].uL;

		ran = (int)(FRANDOM * paramsLc.n_dis);
		J13 = 0.;
		u3c = ww_disLc[ran].u0;
		u1c = ww_disLc[ran].uL;

		
		temp = loopvertex_correlation(J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, 1, 1, 0);

		corrAAC += temp;
		corr_quadAAC += temp * temp;
	}
	}
	prob_disLc = 1. - paramsLc.probJ - paramsLc.prob_con;
	corrAACdouble = (double)corrAAC * paramsLa.probJ * paramsLb.probJ * prob_disLc / ((double)n);
	corr_tot += corrAACdouble;
	error_corrAAC = paramsLa.probJ * paramsLb.probJ * prob_disLc *
		       sqrt(((double)corr_quadAAC / ((double)n) - (double)corrAAC / ((double)n) * (double)corrAAC / ((double)n)) / ((double)n - 1.));

	if (paramsLb.n_dis > 0) 
	{
		for (n = 0; n < (double)N / 4.; n++) 
	{
		ran = (int)(FRANDOM * N);
		J12 = JuuLa[ran].J;
		u1a = JuuLa[ran].u0;
		u2a = JuuLa[ran].uL;

		ran = (int)(FRANDOM * paramsLb.n_dis);
		J23 = 0.;
		u2b = ww_disLb[ran].u0;
		u3b = ww_disLb[ran].uL;

		ran = (int)(FRANDOM * N);
		J13 = JuuLc[ran].J;
		u3c = JuuLc[ran].u0;
		u1c = JuuLc[ran].uL;

		
		temp = loopvertex_correlation(J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, 1, 0, 1);

		corrACA += temp;
		corr_quadACA += temp * temp;
	}
	}
	prob_disLb = 1. - paramsLb.probJ - paramsLb.prob_con;
	corrACAdouble = (double)corrACA * paramsLa.probJ * prob_disLb * paramsLc.probJ / ((double)n);
	corr_tot += corrACAdouble;
	error_corrACA = paramsLa.probJ * prob_disLb * paramsLc.probJ *
		       sqrt(((double)corr_quadACA / ((double)n) - (double)corrACA / ((double)n) * (double)corrACA / ((double)n)) / ((double)n - 1.));

	if (paramsLa.n_dis > 0) 
	{
		for (n = 0; n < (double)N / 4.; n++) 
	{
		ran = (int)(FRANDOM * paramsLa.n_dis);
		J12 = 0.;
		u1a = ww_disLa[ran].u0;
		u2a = ww_disLa[ran].uL;

		ran = (int)(FRANDOM * N);
		J23 = JuuLb[ran].J;
		u2b = JuuLb[ran].u0;
		u3b = JuuLb[ran].uL;

		ran = (int)(FRANDOM * N);
		J13 = JuuLc[ran].J;
		u3c = JuuLc[ran].u0;
		u1c = JuuLc[ran].uL;

		
		temp = loopvertex_correlation(J12, J23, J13, u1a, u2a, u2b, u3b, u3c, u1c, 0, 1, 1);

		corrCAA += temp;
		corr_quadCAA += temp * temp;
	}
	}
	prob_disLa = 1. - paramsLa.probJ - paramsLa.prob_con;
	corrCAAdouble = (double)corrCAA * prob_disLa * paramsLb.probJ  * paramsLc.probJ / ((double)n);
	corr_tot += corrCAAdouble;
	error_corrCAA = prob_disLa * paramsLb.probJ  * paramsLc.probJ *
		       sqrt(((double)corr_quadCAA / ((double)n) - (double)corrCAA / ((double)n) * (double)corrCAA / ((double)n)) / ((double)n - 1.));

	error_corrAAC += error_corrACA + error_corrCAA;

	fprintf(f, "%d %g %g %g %g %g \n", L, -corr_tot, sqrt(error_corrAAA * error_corrAAA + error_corrAAB * error_corrAAB + error_corrAAC * error_corrAAC), corrAAAdouble, corrAABdouble + corrBAAdouble + corrABAdouble, corrAACdouble + corrACAdouble + corrCAAdouble);
}

int main(int argc, char *argv[])
{

	int L = 1, i, iter, n;
	double con, error_con;
	FILE *f;
	char s[60];

	

	Juu_type *JuuLa, *JuuLb, *JuuLc;

	ww_type *ww_conLa, *ww_disLa, *ww_conLb, *ww_disLb, *ww_conLc, *ww_disLc;

	params paramsLa, paramsLb, paramsLc;

	if (argc != 3 && argc != 4)
	{
		fprintf(stderr,
			"usage: %s <fixed degree> <field> \n",
			argv[0]);
		exit(1);
	}
	
	z = atoi(argv[1]);
	field = atof(argv[2]);

	if (argc == 3)
	{ 
		FILE *devran = fopen("/dev/urandom", "r");
		fread(&myrand, 4, 1, devran);
		fclose(devran);
	}
	else
		myrand = atoi(argv[3]); 

	Init_Random();

	ww_conLa = (ww_type *)malloc(N * sizeof(ww_type));
	ww_disLa = (ww_type *)malloc(N * sizeof(ww_type));
	ww_conLb = (ww_type *)malloc(N * sizeof(ww_type));
	ww_disLb = (ww_type *)malloc(N * sizeof(ww_type));
	ww_conLc = (ww_type *)malloc(N * sizeof(ww_type));
	ww_disLc = (ww_type *)malloc(N * sizeof(ww_type));
	new_ww_con = (ww_type *)malloc(N * sizeof(ww_type));
	new_ww_dis = (ww_type *)malloc(N * sizeof(ww_type));
	JuuLa = (Juu_type *)malloc(N * sizeof(Juu_type));
	JuuLb = (Juu_type *)malloc(N * sizeof(Juu_type));
	JuuLc = (Juu_type *)malloc(N * sizeof(Juu_type));
	new_Juu = (Juu_type *)malloc(N * sizeof(Juu_type));

	for (n = N_min; n < N_max; n++)
	{
		sprintf(s, "dati_3pointLL2L/loopvertex_SG_LL2L_%d.dat", n); 
		f = fopen(s, "w");

		
		
		for (i = 0; i < N_BP; i++) 
		{
			u[i] = FRANDOM > 0.5 ? pm1 : 0;
		}
		for (iter = 0; iter < 1000; iter++)
		{ 
			oneStep_BP(z, field);
		}

		Init_triplet(JuuLa, &paramsLa);
		Init_triplet(JuuLb, &paramsLb);
		Init_triplet(JuuLc, &paramsLc);
		L = 1;
		fprintf(f, "# z: %d h: %g N: %d \n", z, field, N);
		fprintf(f, "# 1:L 2:con_loop 3:error 4:con_loopAA 5:con_loopAB 6: con_loopAC \n");
		fflush(f); 

		Lplus1(JuuLc, ww_conLc, ww_disLc, &paramsLc);

		L = 1;
		for (i = 0; i < 8; i++)
		{
			Lplus1(JuuLa, ww_conLa, ww_disLa, &paramsLa);
			Lplus1(JuuLb, ww_conLb, ww_disLb, &paramsLb);
			Lplus1(JuuLc, ww_conLc, ww_disLc, &paramsLc);
            Lplus1(JuuLc, ww_conLc, ww_disLc, &paramsLc);

			L += 1;
			loopvertex(L, JuuLa, ww_conLa, ww_disLa, JuuLb, ww_conLb, ww_disLb, JuuLc, ww_conLc, ww_disLc, paramsLa, paramsLb, paramsLc, f);
			fflush(f);
		}
		fclose(f);
	}


	return 0;
}
