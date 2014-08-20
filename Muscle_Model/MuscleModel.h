#pragma once
#include <stdio.h>

class MuscleModel
{
	int past_spike_1;
	int past_spike_2;
	double current_h;
	double past_h_1;
	double past_h_2;
	double past_T;

	int sample();
	double s(double);
	double active_state(double, double, double, double, double);
	int* spike_train(double, int);
	double* gen_h_diff_eq(int);
	double* plotForce(int, double*);
	double* gen_waveform(double, double, double);
	double d_force(double, double, double, double);

public:
	MuscleModel(void);
	~MuscleModel(void);
	double stepping_model(double, int);
	double T;	
};

