

#include "motif.hpp"


//--pre_post_spikes------------------------------------------------------------
// initialization of the stimulation routine, setting all the parameters and evoking 
// the variables

pre_post_spikes::pre_post_spikes(Parameter& pre_post_spikes_par, const double& scan_pa0, const double& scan_pa1, const double& scan_pa2) {
 
  using String::compose;

  choice      = pre_post_spikes_par.get_param<int>("choice");

  nic_0        = pre_post_spikes_par.get_param<double>("nicotine_rest");
  //nic_applied1 = pre_post_spikes_par.get_param<double>("nicotine_applied1");
  nic_applied1 = scan_pa0;
  nic_applied2 = pre_post_spikes_par.get_param<double>("nicotine_applied2");
  tau_nic      = pre_post_spikes_par.get_param<double>("tau_nicotine");
   
  nu_glu_0       = scan_pa2; //pre_post_spikes_par.get_param<double>("m_pfc_rest");
  //m_pfc_applied1= pre_post_spikes_par.get_param<double>("m_pfc_applied1");
  //m_pfc_applied2= pre_post_spikes_par.get_param<double>("m_pfc_applied2");
    
  ach_0         = scan_pa1; //pre_post_spikes_par.get_param<double>("acetylcholine_rest");
  //ach_applied1  = pre_post_spikes_par.get_param<double>("acetylcholine_applied1");
  //ach_applied2  = pre_post_spikes_par.get_param<double>("acetylcholine_applied2");
  
 
  t_app1_start = pre_post_spikes_par.get_param<double>("t_app1_start");
  t_app1_end   = pre_post_spikes_par.get_param<double>("t_app1_end");
  t_app2_start = pre_post_spikes_par.get_param<double>("t_app2_start");
  t_app2_end   = pre_post_spikes_par.get_param<double>("t_app2_end");
  
  
  alpha7      = pre_post_spikes_par.get_param<double>("alpha7");
  beta2       = pre_post_spikes_par.get_param<double>("beta2");
  w_g         = pre_post_spikes_par.get_param<double>("w_GABA");
  w_glu       = pre_post_spikes_par.get_param<double>("w_Glu");
  w_ach       = pre_post_spikes_par.get_param<double>("w_ACh");
  r           = pre_post_spikes_par.get_param<double>("r_beta2_da");
  s           = pre_post_spikes_par.get_param<double>("s_alpha7_da");
  da_rest     = pre_post_spikes_par.get_param<double>("da_rest");
  gaba_rest     = pre_post_spikes_par.get_param<double>("gaba_rest");
  //m_ppt       = pre_post_spikes_par.get_param<double>("m_ppt");
  glu_da      = pre_post_spikes_par.get_param<double>("glu_DA");
  glu_gaba    = pre_post_spikes_par.get_param<double>("glu_GABA");
  ach_da      = pre_post_spikes_par.get_param<double>("ach_DA");
  ach_gaba    = pre_post_spikes_par.get_param<double>("ach_GABA");
  
  // general initialization
  // desensitization properties
  n_act_beta2 = pre_post_spikes_par.get_param<double>("n_act_beta2");
  K_act_beta2 = pre_post_spikes_par.get_param<double>("K_act_beta2");
  //tau_act_beta2 = pre_post_spikes_par.get_param<double>("tau_act_beta2");
  K_d_beta2   = pre_post_spikes_par.get_param<double>("K_d_beta2");
  n_d_beta2   = pre_post_spikes_par.get_param<double>("n_d_beta2");
  tau_d_max_beta2 = pre_post_spikes_par.get_param<double>("tau_d_max_beta2");
  tau_d_base_beta2 = pre_post_spikes_par.get_param<double>("tau_d_base_beta2");
  K_tau_d_beta2   = pre_post_spikes_par.get_param<double>("K_tau_d_beta2");
  n_tau_d_beta2   = pre_post_spikes_par.get_param<double>("n_tau_d_beta2");
  alpha_beta2     = pre_post_spikes_par.get_param<double>("potency_beta2");

  n_act_alpha7   = pre_post_spikes_par.get_param<double>("n_act_alpha7");
  K_act_alpha7   = pre_post_spikes_par.get_param<double>("K_act_alpha7");
  //tau_act_alpha7 = pre_post_spikes_par.get_param<double>("tau_act_alpha7");
  K_d_alpha7   = pre_post_spikes_par.get_param<double>("K_d_alpha7");
  n_d_alpha7   = pre_post_spikes_par.get_param<double>("n_d_alpha7");
  tau_d_max_alpha7  = pre_post_spikes_par.get_param<double>("tau_d_max_alpha7");
  tau_d_base_alpha7 = pre_post_spikes_par.get_param<double>("tau_d_base_alpha7");
  K_tau_d_alpha7  = pre_post_spikes_par.get_param<double>("K_tau_d_alpha7");
  n_tau_d_alpha7  = pre_post_spikes_par.get_param<double>("n_tau_d_alpha7");
  alpha_alpha7     = pre_post_spikes_par.get_param<double>("potency_alpha7");

  // allocating memory 
  nvar = 10;

  yt   = new double[nvar];
  dydt = new double[nvar];
  yout = new double[nvar];
  dy   = new double[nvar];
  yerr = new double[nvar];
  yscal= new double[nvar];
  
  nic   = nic_0;
  ach   = ach_0;
  nu_glu = nu_glu_0;

  // beta 2
  beta2_act = (pow((alpha_beta2*nic + ach),n_act_beta2)/(pow(K_act_beta2,n_act_beta2) +  pow((alpha_beta2*nic+ach),n_act_beta2)));
  yt[0] = (pow(K_d_beta2,n_d_beta2)/(pow(K_d_beta2,n_d_beta2) + pow(nic,n_d_beta2)));
  // alpha 7
  alpha7_act = (pow((alpha_alpha7*nic + ach),n_act_alpha7)/(pow(K_act_alpha7,n_act_alpha7) +  pow((alpha_alpha7*nic+ach),n_act_alpha7)));
  yt[1] = (pow(K_d_alpha7,n_d_alpha7)/(pow(K_d_alpha7,n_d_alpha7) + pow(nic,n_d_alpha7)));
 
  I_glu =  w_glu*(alpha7*alpha7_act*yt[1] + nu_glu);  
  if (I_glu > w_glu)
    I_glu = w_glu;

  I_ach =  w_ach*beta2*beta2_act*yt[0];
 
  yt[2] = nic;

  da =  da_rest + I_glu*(s*(glu_da + glu_gaba*w_g) - glu_gaba*w_g) + I_ach*(r*(ach_da + ach_gaba*w_g) - ach_gaba*w_g);
  if ( da < 0.) 
	da = 0.;

  gaba = gaba_rest + (1.-s)*glu_gaba*I_glu + (1.-r)*ach_gaba*I_ach; 

 
  da_0          = da_max = da_min = da;
  gaba_0        = gaba_max = gaba_min = gaba;
  da_ex_inp_0   = da_ex_inp_max = s*glu_da*I_glu + r*ach_da*I_ach;
  da_in_inp_0   = da_in_inp_max = w_g*gaba;
  da_old = da_older = da;

  initial_slope = 0.;
  initial_slope_b = 0.;
 
  for (int n=3 ; n<10 ; n++ ) {
	yt[n] = 0.;
  }

  app1 = app2 = false;
  stage = 0;
  decay = rise =  true;
	
  // applied nicotine concentraion and base level ACh and Glu are written in file name
  string s_name_1=compose("nicotine_dynamics_%1_%2_%3.dat",nic_applied1,nu_glu_0,ach_0);
  const char* file_1=s_name_1.c_str();
  motif_nicotine_f.open(file_1);

	
  //motif_nicotine_f.open("nicotine_dynamics.dat");

  //cout.precision(10);

}

//-----------------------------------------------------------------------------
// destructor of the class

pre_post_spikes::~pre_post_spikes() {
  


  motif_nicotine_f.close();
  //motif_spikes_f.close();
  
}



//-----------------------------------------------------------------------------
// calculations performed before the time integration and the time time
// integration

void pre_post_spikes::evaluate_event(const double& t, const double& dt) {
 
  
  
    if ((t > t_app1_start) && (t < t_app1_end)) { 
        app1 = true;
        app2 = false;
        // presynaptic spike
        if ((choice == 0) || (choice == 3) || (choice == 4) || (choice == 6)) 
            nic = nic_applied1;
        if ((choice == 1) || (choice == 4) || (choice == 5) || (choice == 6) )
            nu_glu = nu_glu_applied1;
        if ((choice == 2) || (choice == 3) || (choice == 5) || (choice == 6) )
            ach   = ach_applied1;
    } 
    else if ((t > t_app2_start) && (t < t_app2_end)) {
        app1 = false;
        app2 = true;
        if ((choice == 0) || (choice == 3) || (choice == 4) || (choice == 6)) 
            nic = nic_applied2;
        if ((choice == 1) || (choice == 4) || (choice == 5) || (choice == 6) )
            nu_glu = nu_glu_applied2;
        if ((choice == 2) || (choice == 3) || (choice == 5) || (choice == 6) )
            ach   = ach_applied2;
    }
    else { 
        app1 = false;
        app2 = false;
        nic   = nic_0;
        nu_glu = nu_glu_0;
        ach   = ach_0;
    }
        
  
    // evaluate differential equations
    functions(t,yt,dydt);
    rkck(yt,dydt,t,nvar,dt,yout); 
  
}

//-----------------------------------------------------------------------------
// calculations performed after the time integration 

void pre_post_spikes::update(double t, double ddt) {
  
    beta2_act = (pow((alpha_beta2*yt[2] + ach),n_act_beta2)/(pow(K_act_beta2,n_act_beta2) +  pow((alpha_beta2*yt[2]+ach),n_act_beta2)));
    alpha7_act = (pow((alpha_alpha7*yt[2] + ach),n_act_alpha7)/(pow(K_act_alpha7,n_act_alpha7) +  pow((alpha_alpha7*yt[2]+ach),n_act_alpha7)));
    
    I_glu =  w_glu*(alpha7*alpha7_act*yt[1] + nu_glu);  
              
    if ( I_glu > w_glu ) 
        I_glu = w_glu;

    I_ach =  w_ach*beta2*beta2_act*yt[0];

    da = da_rest + I_glu*(s*(glu_da + glu_gaba*w_g) - glu_gaba*w_g) + I_ach*(r*(ach_da + ach_gaba*w_g) - ach_gaba*w_g);
  	if ( da < 0.) 
		da = 0.;

  	gaba = gaba_rest + (1.-s)*glu_gaba*I_glu + (1.-r)*ach_gaba*I_ach; 


        
    //if (app1) {
        if (da > da_max)
            da_max = da;
        if (da < da_min)
            da_min = da;
        da_end_stim = da;
        da_bend_stim = da;
    
        if (gaba > gaba_max)
            gaba_max = gaba;
		if (gaba < gaba_min)
			gaba_min = gaba;
        if (da_ex_inp_max < (s*glu_da*I_glu + r*ach_da*I_ach) )
            da_ex_inp_max = s*glu_da*I_glu + r*ach_da*I_ach;
        if (da_in_inp_max < w_g*gaba)
            da_in_inp_max = w_g*gaba;

	if ( rise && da/da_0 > ((1.60528 - 1.)/2. + 1.) ) {  // 1.59841
		t_half_rise = t;
		rise = false;
	}

	if ( decay && t > (t_app1_end + 10.) && da < ((da_max-da_0)/2. + da_0) ) {
		t_half_decay = t;
		decay = false;
	}
    
	// excitation (1) - inhibition (2) - excitation (3)
   	if (da > da_0 && stage == 0 ) {
		stage = 1;
	}
	if (da < da_0 && stage == 1 ) {
		stage = 2;
	}
	if (da > da_0 && stage == 2 ) {
		stage = 3;
	}

	// inhibition (4) - excitation (5)
	if (da < da_0 && stage == 0 ) {
		stage = 4;
	} 
	if (da > da_0 && stage == 4 ) {
		stage = 5;
	}
	//cout << (da - da_old) << endl;
	//cout << stage << endl;
	    
    da_older = da_old; 
	da_old = da; 
}
//-----------------------------------------------------------------------------
// deletes the evoked variables from memory

void pre_post_spikes::free() {
  
  delete  [ ] yt;
  delete  [ ] dydt;
  delete  [ ] yout;
  delete  [ ] dy;  
  delete  [ ] yerr;
  delete  [ ] yscal;
}

//-----------------------------------------------------------------------------

// routine with the actual differential equations
void pre_post_spikes::functions(double t, double yt[], double dydt[]) {
  

	tau_d_beta2 = (tau_d_max_beta2*60.*pow(K_tau_d_beta2,n_tau_d_beta2)/(pow(K_tau_d_beta2,n_tau_d_beta2) + pow(yt[2],n_tau_d_beta2))) + tau_d_base_beta2; 
	tau_d_alpha7 = (tau_d_max_alpha7*60.*pow(K_tau_d_alpha7,n_tau_d_alpha7)/(pow(K_tau_d_alpha7,n_tau_d_alpha7) + pow(yt[2],n_tau_d_alpha7))) + tau_d_base_alpha7; 

    // beta 2 desensitization
    dydt[0] = (pow(K_d_beta2,n_d_beta2)/(pow(K_d_beta2,n_d_beta2) + pow(yt[2],n_d_beta2)) - yt[0])/tau_d_beta2;

    // alpha 7 desensitization
    dydt[1] = (pow(K_d_alpha7,n_d_alpha7)/(pow(K_d_alpha7,n_d_alpha7) + pow(yt[2],n_d_alpha7)) - yt[1])/tau_d_alpha7;

	// nicotine dynamics - single exponential with time constant tau_nic
    dydt[2] = (nic - yt[2])/tau_nic;

   	// integration of DA neuron activity, relative from baseline
    dydt[3] = (da - da_0);
	if ( !(da_0 == 0) )
		dydt[3] = dydt[3]/da_0;
    dydt[4] = (gaba - gaba_0)/gaba_0;
	for (int n=1 ; n < 6 ; n++ ) {
		if ( stage == n ) {
			dydt[4+n] = (da - da_0); 
			if ( !(da_0 == 0) ) {
				dydt[4+n] = dydt[4+n]/da_0;
			}
		} 
		else 
			dydt[4+n] = 0.;
     }

    
}

//---------------------------------------------------------------------------------------------------------------------------

// the runge-kutta routine can be run with dynamics step-size and predefined precission. 
// I did the simulations with fixed step size, hence all the commands commented out 
void pre_post_spikes::rkqs(double &t, const double htry, const double eps, double &hdid, double &hnext) {
  
  //const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4, TINY=1.0e-5;
  int i;
  double h;
 
  h=htry;
/*  for (;;) {
    if (pre_spike) {
      yt[1] -= C_pre;
    }
    if (post_spike) {
      yt[2] -= C_post;
    }
    pre_spike = pre_s  = false;
    post_spike= post_s = false;*/
    evaluate_event(t,h);
/*    // scaling used to monitor accuracy 
    for (i=0; i<nvar; i++)
      yscal[i] = fabs(yt[i]) + fabs(dydt[i]*h)+TINY;
    errmax=0.0;
    for (i=0;i<nvar;i++) {
      errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
      //errmax=MAX(errmax,fabs(yerr[i]));
    }
    errmax /= eps;
    if (errmax <= 1.0) break;
    htemp=SAFETY*h*pow(errmax,PSHRNK);
    h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
    tnew=t+h;
    if (tnew == t) nrerror("stepsize underflow in rkqs");
  }*/
//  if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
//  else hnext=5.0*h;
  t += (hdid=h);
  dt_done = h;
  for (i=0;i<nvar;i++) yt[i]=yout[i];
  //pre_s = post_s = false;
  update(t,h);
}

//-----------------------------------------------------------------------------

// routine performs 4-order runge-kutta step
void pre_post_spikes::rkck(const double y[], const double dydx[], const double x, const int nva,
		    const double h, double yout[]) 
{
  static const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
    b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
    b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
    b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
    b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
    c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
    dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
  int i;
  
  //int n=y.size();
  double ak2[nva],ak3[nva],ak4[nva],ak5[nva],ak6[nva],ytemp[nva];
  
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  functions(x+a2*h,ytemp,ak2);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  functions(x+a3*h,ytemp,ak3);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  functions(x+a4*h,ytemp,ak4);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  functions(x+a5*h,ytemp,ak5);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  functions(x+a6*h,ytemp,ak6);
  for (i=0;i<nva;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0;i<nva;i++) {
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  }
}

//---------------------------------------------------------------------------------

// writes dynamics in file 
void pre_post_spikes::write_files(const double& t) {

  motif_nicotine_f << t << "\t" << yt[2] << "\t" << da << "\t" << gaba  << "\t" << da/da_0 << "\t" << gaba/gaba_0 << "\t" << stage; 
  for (int i=0;i<nvar;i++) 
    motif_nicotine_f << "\t" << yt[i];
  motif_nicotine_f << "\t" << I_glu << "\t" << I_ach << "\t" << beta2_act << "\t" << alpha7_act << "\t" <<  app1 << "\t" << app2 << "\t" << nic << "\t" << ach << "\t" << nu_glu << "\t" << (glu_da*I_glu + r*ach_da*I_ach)/da_ex_inp_0 << "\t" << w_g*gaba/da_in_inp_0 <<  endl; 

}

