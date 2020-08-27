

#include <iostream>
#include <fstream>
#include <cmath>     
#include "nrutils.hpp"
#include "parameter.hpp"
#include "motif.hpp"

using namespace std;

//-main----------------------------------------------------------------------------------

int main ()
{

  // opening parameter file, parameters can then be read in from "par"
  Parameter par("par_nicotine.par");

  // declaration of variables
  double t; 
  double scan_0_param, scan_0_start, scan_0_end, delta_0_scan;
  double scan_1_param, scan_1_start, scan_1_end, delta_1_scan;
  double scan_2_param, scan_2_start, scan_2_end, delta_2_scan;
  long scan_0_steps, scan_1_steps, scan_2_steps;
  int stage_old;
  

  stage_old = 0;
  // opening file streams
  ofstream final_state_f("nicotine_scan.dat");

    
  scan_0_start  = par.get_param<double>("nicotine_applied1_start");
  scan_0_end    = par.get_param<double>("nicotine_applied1_end");
  scan_0_steps  = par.get_param<long>("nicotine_applied1_steps");
  scan_1_start  = par.get_param<double>("ach_1_start");
  scan_1_end    = par.get_param<double>("ach_1_end");
  scan_1_steps  = par.get_param<long>("ach_1_steps");
  scan_2_start  = par.get_param<double>("glu_2_start");
  scan_2_end    = par.get_param<double>("glu_2_end");
  scan_2_steps  = par.get_param<long>("glu_2_steps");
 
  delta_0_scan=double((scan_0_end-scan_0_start)/scan_0_steps);
  delta_1_scan=double((scan_1_end-scan_1_start)/scan_1_steps);
  delta_2_scan=double((scan_2_end-scan_2_start)/scan_2_steps);
  
 
  // parameter which don't require a reinitialization for each single scan
  const double t_start   = par.get_param<double>("t_start");
  const double t_end     = par.get_param<double>("t_end");
  const double h         = par.get_param<double>("dt");
  const long write_steps = par.get_param<long>("write_steps");
  const double epsilon   = par.get_param<double>("epsilon");
  double dt, dtdid, dtnext;
  long next_write; 
  t = t_start;
  
  // nicotine loop
  scan_0_param = scan_0_start;
  for (int l=0; l<=scan_0_steps; l++) {
	scan_2_param = scan_2_start;
	// Glu loop
	for (int k=0; k<=scan_2_steps; k++) {  
		scan_1_param=scan_1_start;
			
		// ACh loop 
		for (long j=0; j<=scan_1_steps; ++j) {
			// setting up the simulations, initializing variables and routine
			double t   = t_start;
			next_write = write_steps;
			dt         = h;
			dtdid      = h;
			dtnext     = h;
			pre_post_spikes pre_post(par,scan_0_param, scan_1_param,scan_2_param);
			pre_post.write_files(t);  
			
			// temporal integration
			for (long i=1;; ++i) { 
				// performs time step of step-width "dt"
				pre_post.rkqs(t, dt, epsilon, dtdid, dtnext);
				// writes dynamical variables, every "next_write" steps
				if (i>=next_write) {
					pre_post.write_files(t); 
					next_write += write_steps;
				}
				// quits temporal integration if predefined end time is reached
				if ( t >= t_end) {
					break;
				}
				dt = dtnext;
			}
			
			// writes max, min and integrals in file
			final_state_f << scan_0_param << "\t" << scan_1_param << "\t" << scan_2_param << "\t" << pre_post.stage;
			if ( pre_post.da_0 == 0 ) {
				final_state_f  << "\t" <<  pre_post.da_max << "\t" << pre_post.da_min; 
			}
			else {
				final_state_f  << "\t" <<  pre_post.da_max/pre_post.da_0 << "\t" << pre_post.da_min/pre_post.da_0; 
			}
			final_state_f << "\t" << pre_post.gaba_max/pre_post.gaba_0 << "\t" << pre_post.gaba_min/pre_post.gaba_0 << "\t" <<
			pre_post.da_ex_inp_max/pre_post.da_ex_inp_0 << "\t" << pre_post.da_in_inp_max/pre_post.da_in_inp_0 << "\t" << pre_post.da_0 << "\t" << pre_post.gaba_0  << "\t" << pre_post.t_half_decay  << "\t" << pre_post.t_half_rise;
			for (int n = 3 ; n < 10 ; n++ ) {
				final_state_f << "\t" <<  pre_post.yt[n]; 
			} 
			final_state_f << endl;
			// frees parameters
			pre_post.free();
	
			scan_1_param+=delta_1_scan;
		}
		scan_2_param += delta_2_scan;
		
		if (!(scan_2_steps==0)) {
			final_state_f << endl;
		}	
	}
	scan_0_param += delta_0_scan;
  }
  
  final_state_f.close();
  
  return (0);
}
