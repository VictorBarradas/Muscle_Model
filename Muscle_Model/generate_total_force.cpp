#include <stdio.h>

const int SAMPLING_RATE = 1024;

double* gen_waveform(double, double, double);
double d_force(double, double, double, double);
double* plotForce(int, double*);
double* gen_h_diff_eq(int);
int* spike_train(double, int);
double active_state(double, double, double, double, double);
double s(double);

int main(){
	int i;
	double *x;
	double *T_list;
	double T_list_1[4000];
	double T_list_2[4000];
	double T_list_3[4000];
	double T_list_4[4000];
	double T_list_5[4000];
	double T_list_6[4000];

	FILE *file;

	x = gen_waveform(1.0, 1.011, 1.0);
	
	T_list = plotForce(5, x);
	for(i = 0; i <= SAMPLING_RATE; i++){
		T_list_1[i] = T_list[i];
	}
	T_list = plotForce(10, x);
	for(i = 0; i <= SAMPLING_RATE; i++){
		T_list_2[i] = T_list[i];
	}
	T_list = plotForce(20, x);
	for(i = 0; i <= SAMPLING_RATE; i++){
		T_list_3[i] = T_list[i];
	}
	T_list = plotForce(50, x);
	for(i = 0; i <= SAMPLING_RATE; i++){
		T_list_4[i] = T_list[i];
	}
	T_list = plotForce(60, x);
	for(i = 0; i <= SAMPLING_RATE; i++){
		T_list_5[i] = T_list[i];
	}
	T_list = plotForce(100, x);
	for(i = 0; i <= SAMPLING_RATE; i++){
		T_list_6[i] = T_list[i];
	}

	file = fopen("test.txt", "w");
	
	for(i = 0; i <= SAMPLING_RATE; i++){
		fprintf(file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i, T_list_1[i], T_list_2[i], T_list_3[i], T_list_4[i], T_list_5[i], T_list_6[i]);
	}

	fclose(file);

	return 0;
}

double s(double x){
	double weighted;

	if (x < 0.5){
		weighted = 0.0;
	}
	else if (x < 1.0){
		weighted = -4.0*x*x + 8.0*x - 3.0; 
	}
	else if (x < 2.0){
		weighted = -x*x + 2.0*x;
	}
	else{
		weighted = 0.0;
	}



	return weighted;
}

double active_state(double x_i0, double x_i1, double x_i2, double y_i1, double y_i2){
	const double b0 = 2299.0;
	const double b1 = -2289.0;
	const double b2 = 0.0;
	const double a0 = 1.0;
	const double a1 = -1.942;
	const double a2 = 0.943;

	double t0;
	double t1;
	double t2;
	double t3;
	double t4;
	double y_i0;

	t0 = b0 * x_i0;
	t1 = b1 * x_i1;
	t2 = b2 * x_i2;
	t3 = a1 * y_i1;
	t4 = a2 * y_i2;

	y_i0 = t0 + t1 + t2 - t3 - t4;

	return y_i0;
}

int* spike_train(double T, int firing_rate){
	int max_n;
	int i;
	static int dtrain[4000];
	
	max_n = T*SAMPLING_RATE;
	for(i = 0; i <= max_n; i++){
		if((i - 1) % int(1.0 / firing_rate*SAMPLING_RATE) == 0){
			dtrain[i] = 1;
		}
		else{
			dtrain[i] = 0;
		}
	}
	return dtrain;
}

//Working
double* gen_h_diff_eq(int firing_rate){
	int *spikes;
	int i;
	int spike_i1;
	int spike_i2;
	static double h[4000];
	double h_i1;
	double h_i2;
	double T;

	spike_i1 = 0;
	spike_i2 = 0;
	h_i1 = 0.0;
	h_i2 = 0.0;
	T = 1.0;
	spikes = spike_train(T, firing_rate);
	
	for(i = 0; i <= T*SAMPLING_RATE; i++){
		h[i] = active_state(*(spikes + i), spike_i1, spike_i2, h_i1, h_i2);
		spike_i2 = spike_i1;
		spike_i1 = *(spikes + i);
		h_i2 = h_i1;
		h_i1 = h[i];
	}

	return h;
}

double* plotForce(int firing_rate, double* x){
	int i;
	double T_i1;
	double t;
	double A;
	static double T[4000];
	double dT[4000];
	double *h;

	T_i1 = 0;
	T[0] = 0.0;
	t = 1.0;
	h = gen_h_diff_eq(firing_rate);

	for(i = 0; i <= t * SAMPLING_RATE; i++){
		//A = *(h + i) * s(*(x + i));
		A = h[i] * s(x[i]);
		dT[i] = d_force(T_i1, 1.0, 0.0, A);
		T[i] = T_i1 + dT[i] * (1.0 / SAMPLING_RATE);
		T_i1 = T[i];
	}

	return T;
}

//Working
double* gen_waveform(double L1, double L2, double T){
	static double x[4000] = {0};
	double dt;
	double ramp_speed;
	double rising_range;
	double x_up[4000] = {0};
	int max_n;
	int i;

	dt = 1.0 / SAMPLING_RATE;
	max_n = (T * SAMPLING_RATE) - 1;
	ramp_speed = 20.0;
	rising_range = (L2 - L1) / ramp_speed / dt;

	for(i = 0; i < max_n; i++){
		if(i <= rising_range){
			x_up[i] = L1 + i*(L2 - L1) / rising_range;
			x[i] = x_up[i];
		}
		/*if(i > 0){*/
		else{
			x[i] = L2;
		}
	}
	x[0] = 1.0;

	//for(i = 0; i < max_n - 1; i++){
	//	//First max_n elements of x are a step-like function. Next (max_n - 1) are
	//	//its derivative (dx)
	//	x[max_n + i] = (x[i + 1] - x[i] + x[max_n - 1])/dt;
	//}

	return x;
}

double d_force(double T_0, double lce, double vel, double A){
	const double L0 = 1.0;
	const double KSE = 136.0;
	const double KPE = 75.0;
	const double B = 50.0;

	double dT_0;

	if (lce > L0){
		dT_0 = KSE / B * (KPE * (lce - L0) + B * vel - (1 + KPE/KSE)*T_0 + A);
	}
	else{
		dT_0 = KSE / B * (KPE * (lce - L0) + B * vel - (1 + KPE/KSE)*T_0 + A);
	}

	return dT_0;

}
