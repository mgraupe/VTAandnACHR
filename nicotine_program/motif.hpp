
#ifndef motif_H
#define motif_H

#include <fstream>
#include <cmath>
#include <string>

#include "parameter.hpp"
#include "compose.hpp"
#include "nr.hpp"

using namespace NR;

//-pre_post_spikes-----------------------------------

class pre_post_spikes {
public:
  pre_post_spikes(Parameter&, const double&, const double&, const double&);
  ~pre_post_spikes();
  void evaluate_event(const double&, const double&);
  void rkqs(double &, const double, const double, double &, double &);
  void write_files(const double&);
  void free();
  double da, da_b, gaba;
  double *yt;
  double da_old, da_older;
  double da_max, gaba_max, da_min, gaba_min, da_0, gaba_0, da_end_stim, da_bend_stim;
  double da_bmax, da_bmin, da_b0;
  double da_ex_inp_max, da_in_inp_max, da_ex_inp_0, da_in_inp_0, da_bex_inp_max, da_bex_inp_0; 
  double initial_slope, initial_slope_b, t_half_decay, t_half_rise; 
  int stage; 
  
private:   
  void functions(double , double [], double []); 
  void rk4(double [], double [], int , double, double, double []);
  void rkck(const double [], const double [], const double, const int, const double, double []);
  void update(double,double);
  
  double dt, dt_done;
  double nic, nic_0, nic_applied1, nic_applied2, t_app1_start, t_app1_end, t_app2_start, t_app2_end, r,s;
  double ach, ach_0, ach_applied1, ach_applied2;
  double nu_glu, nu_glu_0, nu_glu_applied1, nu_glu_applied2, tau_act_glu;
  double da_rest, gaba_rest;
  double w_g, w_glu, w_ach, I_glu, I_ach, I_da;
  double tau_nic, glu_da, glu_gaba, ach_da, ach_gaba;
  double alpha7, beta2;
  double beta2_act, alpha7_act;
  double n_act_beta2, K_act_beta2, K_d_beta2, n_d_beta2, tau_d_max_beta2, K_tau_d_beta2, n_tau_d_beta2, tau_act_beta2, tau_d_base_beta2, tau_d_beta2;
  double n_act_alpha7, K_act_alpha7, K_d_alpha7, n_d_alpha7, tau_d_max_alpha7, K_tau_d_alpha7, n_tau_d_alpha7, tau_act_alpha7, tau_d_base_alpha7, tau_d_alpha7;
  int choice, phase;
  bool app1, app2, decay, rise;
  double alpha_beta2, alpha_alpha7;

  
  double epsilon;

  int nvar;
  
  double *dydt, *yout, *dy, *yerr, *yscal;

  ofstream motif_nicotine_f;
  //ofstream motif_camkII_f;
};

#endif




