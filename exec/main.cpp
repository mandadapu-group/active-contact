// created by Amaresh Sahu
// amaresh.sahu@berkeley.edu
// 27 Sept. 2018



#include "params.h"
using namespace std;


// functions and descriptions provided below

double calcOmega(double k);
double calcPassiveAnalyticalResult(double k);
double calcPassiveAveragedAnalyticalResult(double k);
double calcAnalyticalResult(double k);
double calcActiveTimeFactor(double t, double t_p, double t_dur);
double calcActivePressureFactor();


/**
 * main method, performs all the
 * analysis and outputs results
 *
 */
int main(int argc, char** argv)
{


	/* *** create data containers *** */
	

	// simulated height (active)
	vector< vector<double> >
		activeHeightData(NUM_MODES, vector<double>(NUM_MODES));

	// simulated height (passive)
	vector< vector<double> >
		passiveHeightData(NUM_MODES, vector<double>(NUM_MODES));

	// passive analytical height
	vector< vector<double> >
		passiveAnalytical(NUM_MODES, vector<double>(NUM_MODES));

	// active analytical height
	vector< vector<double> >
		activeAnalytical(NUM_MODES, vector<double>(NUM_MODES));



	// for active simulations, the pressure
	// factor p_bar dictates the pressure
	// which the active particles exert on
	// the membrane
	double p_bar = calcActivePressureFactor();



	// random number generator for white noise:
	// Gaussian with mean 0, variance 1
	default_random_engine gaussianGenerator;
	normal_distribution<double>
		gaussianDistribution(0.0, 1.0);



	// output header, when required
	
#if SHOW_OUTPUT
	cout << "m\t";
#if INCLUDE_ACTIVE
	cout << "act_sim\t\t\t";
	cout << "act_theo\t";
#else
	cout << "pass_sim\t\t";
	cout << "pass_theo_exact\t";
#if AVERAGE_MODES
	cout << "pass_theo_avg\t";
#endif
#endif
	cout << endl;

#endif



	/* *** generate collision data *** */


#if INCLUDE_ACTIVE

	// calculate the number of collision events
	int num_collisions = (int) (
			NUM_PARTICLES * (
				FINAL_TIME / (
					TRAVERSAL_TIME + TAU_R
					) ) );

	// initialize empty vectors for collision data
	vector<double> collisionTimes(num_collisions);
	vector<double> collisionX(num_collisions);
	vector<double> collisionY(num_collisions);


	// randomly generate collision data:
	// times and locations (x and y)
	for (int i = 0; i < num_collisions; ++i) {

		// times are between 0 and final time
		collisionTimes[i]
			= FINAL_TIME * rand() / RAND_MAX;

		// positions are between -L/2 and L/2
		// in each direction
		collisionX[i]
			= LENGTH * rand() / RAND_MAX - LENGTH/2.0;
		collisionY[i]
			= LENGTH * rand() / RAND_MAX - LENGTH/2.0;
	}

#endif
	


	/* *** main loop through all modes *** */


	// modes in x direction
	for (int m = 0; m < NUM_MODES; ++m) {
		// modes in y direction
		for (int n = 0; n < NUM_MODES; ++n) {


			// ignore the zero mode
			if (m == 0 && n == 0) {
				break;
			}


			// calculate wave vector, relaxation
			// frequency, and appropriate time
			// scale
			
			double k_x = 2. * M_PI * m / LENGTH;
			double k_y = 2. * M_PI * n / LENGTH;
			double k   = sqrt(k_x*k_x + k_y*k_y);
			double omega = calcOmega(k);
			double delta_t = 0.01 / omega;



			// populate the passive analytical
			// solution for all modes
			passiveAnalytical[m][n]
				= calcPassiveAnalyticalResult(k);

			// populate the active analytical
			// solution for all modes
			activeAnalytical[m][n]
				= calcAnalyticalResult(k);



			// the fluctuation spectrum calculation
			// is only relevant at low mode numbers.
			// For higher modes, all activity is lost
			// and we use the theoretical result.

			if (n > NUM_SMALL_MODES || m > NUM_SMALL_MODES) {
				passiveHeightData[m][n] = passiveAnalytical[m][n];
				activeHeightData[m][n] = passiveAnalytical[m][n];
			} else {


				// real and imaginary parts of the
				// Fourier transform of the height,
				// in units of nm^3
				double re_h = 0.0;
				double im_h = 0.0;

				// the cumulative sum of |h|^2 over
				// all times
				double passive_h_squared_total = 0.0;
				double active_h_squared_total = 0.0;


				// calculate the variance in height
				// due to random noise, for a given
				// mode, in units of nm^3
				double height_fluct
					= sqrt(KBT * delta_t / (4.0 * VISCOSITY * k));


				// for a single mode, loop through
				// all times and calculate the
				// height fluctuations


				int num_time_steps = 0;

				for (double t = 0.0; t < FINAL_TIME; t += delta_t) {


					num_time_steps++;


					// evolve the real and imaginary
					// height components according to
					// the discretized Langevin eqn
					re_h = re_h * (1. - omega * delta_t)
						+ height_fluct *
						gaussianDistribution(gaussianGenerator);

					im_h = im_h * (1. - omega * delta_t)
						+ height_fluct *
						gaussianDistribution(gaussianGenerator);


					// passive height fluctuation sum
					passive_h_squared_total += re_h*re_h + im_h*im_h;


					// the evolution of the Fourier
					// height modes is modified in the
					// presence of active particles
#if INCLUDE_ACTIVE


					// loop through all collisions,
					// and calculate the additional
					// contributions

					for (int collIdx = 0; collIdx < num_collisions; ++collIdx) {

						// particle collision time
						double t_p = collisionTimes[collIdx];

						// calculation of phi(t, t_p; tau_R)
						double phi = calcActiveTimeFactor(t, t_p, TAU_R);


						// if there is no temporal component,
						// there is no contribution from the
						// collision at this time

						if (phi > 0) {


							// x and y collision location
							double x_p = collisionX[collIdx];
							double y_p = collisionY[collIdx];


							// the exponential prefactor of the
							// exponential
							double active_pressure
								= phi * p_bar * M_PI * pow(GAUSSIAN_LENGTH,2.) *
								delta_t / 2. / VISCOSITY / k / LENGTH;

							// multiply with the exponential
							active_pressure *=
								exp(-0.5 * pow(GAUSSIAN_LENGTH, 2) * k * k);



							// calculate the real and imaginary parts of
							// e^{-i \bm{\rho}^p \cdot \bm{k}}
							double phase = x_p * k_x + y_p * k_y;
							re_h += active_pressure * cos(phase);
							im_h -= active_pressure * sin(phase);

						}
					}
#endif

					// active height fluctuation sum
					active_h_squared_total += re_h*re_h + im_h*im_h;

				}

				// write the average value of
				// |h|^2 to the data container
				passiveHeightData[m][n]
					= passive_h_squared_total / num_time_steps;

				// write the average value of
				// |h|^2 to the data container
				activeHeightData[m][n]
					= active_h_squared_total / num_time_steps;

			}
		}
	}



	/* *** show output, without mode averaging *** */

#if SHOW_OUTPUT && !AVERAGE_MODES


	for (int m = 0; m < NUM_MODES; ++m) {
		for (int n = 0; n < NUM_MODES; ++n) {

			// do not output the zero mode,
			// as it was not calculated
			if (m == 0 && n == 0) {
				break;
			}

			// output data
			double mode = sqrt(m*m + n*n);
			cout << mode << "\t";
			double passiveHeightVal = passiveHeightData[m][n] / pow(VESICLE_RADIUS, 4);
			double passiveAnalyticalVal = passiveAnalytical[m][n] / pow(VESICLE_RADIUS, 4);
			cout << passiveHeightVal << "\t";
			cout << passiveAnalyticalVal << "\t";
#if INCLUDE_ACTIVE
			double activeHeightVal = activeHeightData[m][n] / pow(VESICLE_RADIUS, 4);
			cout << activeHeightVal << "\t";
#endif
			cout << endl;

		}
	}

#endif


	/* *** show output, with mode averaging *** */

#if SHOW_OUTPUT && AVERAGE_MODES


	// average of calculated height (passive)
	vector<double> passiveHeightDataAvg(NUM_MODES);

	// average of calculated height (active)
	vector<double> activeHeightDataAvg(NUM_MODES);

	// average of analytical passive height
	vector<double> passiveAnalyticalAvg(NUM_MODES);

	// average of analytical active height
	vector<double> activeAnalyticalAvg(NUM_MODES);

	// vector of the exact analytical value,
	// calculated from theory rather than
	// manually integrated
	vector<double> passiveAnalyticalTheory(NUM_MODES);


	// manually average
	for (int m = 0; m < NUM_MODES; ++m) {
		for (int n = 0; n < NUM_MODES; ++n) {
			if (m > 0 || n > 0) {
				// the factor of 1/pi comes from the expression
				// in calculating the fluctuation spectrum at
				// k_y = 0
				passiveHeightDataAvg[m] += passiveHeightData[m][n] / M_PI;
				passiveAnalyticalAvg[m] += passiveAnalytical[m][n] / M_PI;
#if INCLUDE_ACTIVE
				activeHeightDataAvg[m] += activeHeightData[m][n] / M_PI;
				activeAnalyticalAvg[m] += activeAnalytical[m][n]  / M_PI;
#endif
			}
		}
	}


	// the 2 pi / L comes from us integrating over dk_y, not dm
	for (int m = 0; m < NUM_MODES; ++m) {
		passiveHeightDataAvg[m] *= 2.0 * M_PI / LENGTH;
		passiveAnalyticalAvg[m] *= 2.0 * M_PI / LENGTH;
#if INCLUDE_ACTIVE
		activeHeightDataAvg[m]  *= 2.0 * M_PI / LENGTH;
		activeAnalyticalAvg[m]  *= 2.0 * M_PI / LENGTH;
#endif
	}


	// calculate the average from theory
	for (int n = 1; n < NUM_MODES; ++n) {

		double k_y = 2. * M_PI * n / LENGTH;
		passiveAnalyticalTheory[n]
			= calcPassiveAveragedAnalyticalResult(k_y);

	}

	passiveAnalyticalTheory[0] = passiveAnalyticalTheory[1];


	// non-dimensionalize
	for (int m = 0; m < NUM_MODES; ++m) {
		passiveHeightDataAvg[m]    /= pow(VESICLE_RADIUS, 3);
		passiveAnalyticalTheory[m] /= pow(VESICLE_RADIUS, 3);
		passiveAnalyticalAvg[m]    /= pow(VESICLE_RADIUS, 3);
#if INCLUDE_ACTIVE
		activeHeightDataAvg[m]     /= pow(VESICLE_RADIUS, 3);
		activeAnalyticalAvg[m]     /= pow(VESICLE_RADIUS, 3);
#endif
	}


	// output
	for (int m = 2; m < NUM_MODES; ++m) {
		cout << m << "\t";
#if INCLUDE_ACTIVE
		cout << activeHeightDataAvg[m] << "\t\t";
		cout << activeAnalyticalAvg[m] << "\t\t";
#else
		cout << passiveHeightDataAvg[m] << "\t\t";
		cout << passiveAnalyticalTheory[m] << "\t\t";
		cout << passiveAnalyticalAvg[m] << "\t\t";
#endif
		cout << endl;
	}


#endif

	return 0;
}



/**
 * calculate omega(k), in units of 1/us,
 * for k in units of 1/nm
 *
 * -> input: wave vector magnitude
 * <- output: dimensional omega(k)
 *
 */
double calcOmega(double k) {

	double omega = (0.5 * K_BENDING * k*k*k + LAMBDA * k);
	omega /= (4.0 * VISCOSITY);

	return omega;
}



/**
 * calculate the exact, passive,
 * analytical result for the height
 * fluctuation spectrum, WITHOUT
 * averaging over k_y modes
 *
 */
double calcPassiveAnalyticalResult(double k) {

	return KBT / (0.5*K_BENDING*pow(k,4.) + LAMBDA*k*k);
}



/**
 * calculate the exact, passive,
 * analytical result for the height
 * fluctuation spectrum, WITH
 * averaging over k_y modes
 *
 */
double calcPassiveAveragedAnalyticalResult(double k_y) {

	return KBT / 2. / LAMBDA *
			(1./k_y - 1./sqrt(k_y*k_y + 2. * LAMBDA / K_BENDING));
}



/**
 * calculate the exact, active, analytical result
 * from theoretical calculation
 *
 */
double calcAnalyticalResult(double k) {

	double passive = calcPassiveAnalyticalResult(k);
	double p_bar = calcActivePressureFactor();

	double temp = NUM_PARTICLES * TAU_R / (TAU_R + TRAVERSAL_TIME);
	temp *= pow(pow(GAUSSIAN_LENGTH, 2.0) * p_bar / VESICLE_RADIUS * passive / KBT, 2.0);
	temp *= exp(-k*k*pow(PARTICLE_RADIUS,2.));

	return passive + temp;
}



/**
 * calculate the temporal factor of the
 * active forcing (dimensionless).
 *
 * -> input t: time in simulation
 * -> input t_p: time of collision
 * -> input t_dur: collision duration
 * 		(formerly reorientation time)
 *
 * <- output: active time factor,
 *    between 0 and 1
 *
 */
double calcActiveTimeFactor(
		double t, double t_p, double t_dur) {

	// time interval between
	// collision and simulation
	double t_diff = fabs(t - t_p);

	// flat top of trapezoid
	if (t_diff < t_dur / 2.0) {
		return 1.0;
	}
	// slanted sides
	else if (t_diff < t_dur / 2.0 + TAU_P) {
		return 1. - (t_diff - TAU_R/2.0)/TAU_P;
	}
	// zero otherwise
	else {
		return 0.;
	}

}


/**
 * calculate the pressure factor of
 * active forcing, according to the
 * simple choice of 2 * lambda / a
 *
 */
double calcActivePressureFactor() {

	double p_bar;
	p_bar = 2 * LAMBDA / GAUSSIAN_LENGTH;
	return p_bar;

}
