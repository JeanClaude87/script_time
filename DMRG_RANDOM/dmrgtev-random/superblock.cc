//	superblock.cc			(C) 2011 Fabio Ortolani	fbo 110614
//	==================================================================
#include "timing.hh"
#include "superaction.hh"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
//
//============================================================================
void densityop_parse (size_t, size_t);	       		  	  // [menu.cc]
void ham_parse (size_t);				       	  // [menu.cc]
void hamt_parse (size_t, double);	       		       	  // [menu.cc]
void properties (Block &, Block &); 
void properties (Block &, Block &, 
		 const Action &, const Action &, 
		 double, double);      			    // [superblock.cc]
static void projection (Action &, Block &, Block &,
			const Action &, Block &, Block &);  // [superblock.cc]
void projectionold (Action &, Block &, Block &,
			const Action &, Block &, Block &);  // [superblock.cc]
void evolution  (Block &, Block &);	         	    // [superblock.cc]
void updatebase (Block &, Block &, long);		    // [superblock.cc]
//============================================================================
//static const size_t csize = 					// [action.cc]
//  (sizeof (complex<double>) + sizeof (double) -1) / sizeof (double);
//
extern Barray		blocklft;				  // [dmrg.cc]
extern Barray		blockrgt;				  // [dmrg.cc]
extern size_t		show_hamiltonian;			  // [dmrg.cc]
extern size_t		reflection;				  // [dmrg.cc]
extern size_t		old_sys;				  // [dmrg.cc]
extern size_t		old_uni;				  // [dmrg.cc]
extern size_t		old_sites;				  // [dmrg.cc]
extern	double		n_cutlft;			      // [dmrg.cc]
extern	double		n_cutrgt;				  // [dmrg.cc]
static size_t		sites_old 		= 0;		 
static long		showham 		= 0;
static Superaction	hamaction;
static Superaction      timeaction;
static Superaction      initialaction;
//
double			tensorweight		= 1.e-4;
Apoli		       	hamiltonian;
Apoli                   timehamiltonian;
Aproperty               start_action;
size_t                  ttarget_id;
size_t                  ttarget_index;
double                  timestep                = 0.0;
size_t                  steps_number		= 100;
size_t		        step_divisions          = 10;
size_t                  time_zips               = 1;
vector<Brule>		lfttimerule;
vector<Brule>		rgttimerule;
//
Aarray			super_state;
Qtarget			super_target;
vector<Aproperty>	correlation;
vector<Aproperty>	densityoperator;
vector<double>          super_weight;
//
//	Lanczos default parameters
//
double			zero_energy		= 1.0e-10;
double			zero_norm		= 1.0e-10;
size_t			lanczos_iters		= 30;
size_t			lanczos_repeats		= 500;
size_t			lanczos_strategy	= 1;
//
size_t			initial_guess		= 0;
//
long			purelft 		= -1;
long			purergt			= -1;
Action 			hlft;
Action			hrgt;
double *		energylft		= 0;
double *		energyrgt		= 0;
size_t			edimlft 		= 0;
size_t			edimrgt			= 0;
double 			frustlft		= 0.0;
double			frustrgt		= 0.0;
//
extern double	       	initial_residual;  		   			// [superaction.cc]
extern double	       	initial_value;	                   	// [superaction.cc]
//
vector<double>          renyiarg;
vector<double>		renyientropy;
//
extern bool            	reflect_blocks;
extern bool             reflect_universe;
extern size_t           min_lft;
extern size_t           max_lft;
extern double           n_cutlft;
extern double           n_cutrgt;
extern size_t		min_rgt;
extern size_t		max_rgt;
extern double		truncation;   
extern bool		show_states;
extern bool		show_selection;
extern size_t		show_density;
extern size_t		idop_base;
extern vector<Brule>	lftrule;
extern vector<Brule>	rgtrule;
extern vector<size_t>	symmetry;			          	// [dmrg.cc]
extern Storage		report;					  			// [dmrg.cc]
extern vector<double>	autovalori;
//============================================================================
void trotterlft (Action & result, Action & U, Action & state,
		  Block & system, Block & universe)
{
  Action Blft = system .base ();
  Action bb   = Blft;
  Blft .dagger ();
  Blft *= U;
  Blft *= bb;
  Blft *= state;
  result = Blft;
}
//____________________________________________________________________________
void trotterrgt (Action & result, Action & U, Action & state,
		  Block & system, Block & universe)
{
  Action Ut = U;
  Ut .transpose ();
  Action Brgt = universe .base ();
  Action bb = Brgt;
  bb .dagger ();
  Action aa = state;
  aa *= bb;
  aa *= Ut;
  aa *= Brgt;
  result = aa;
}
//____________________________________________________________________________
void trotterstep (Action & result, Action & U, Action & state,
		  Block & system, Block & universe)
{
  Action Blft, Brgt;
  long sitestates = blocklft [1] .states ();
  //		state :	< a_M, b_N | psi >
  //					M = n. siti system (left)
  //					N = n. siti universe (right)
  Blft = system .base ();	// < a_{M-1}, sigma | a_M> 
  Brgt = universe .base ();	// < lambda, b_{N-1} | b_N>
  Brgt .transpose ();
  multiply (Blft, Blft, state);
  multiply (Blft, Blft, Brgt);
  //	      	Blft :  < a_{M-1} sigma; lambda b_{N-1} | psi >
  size_t na, na1, nb, nb1, ia, ia1, sigma, sigma1, jb, jb1, lambda, lambda1,
    sl, sl1, im, jm;  
  na = Blft .height ();
  nb = Blft .width ();
  na1 = na / sitestates;
  nb1 = nb / sitestates;
  size_t nd = Blft .height () * Blft .width ();
  size_t nu = U .height () * U .width ();
  complex<double> * mm = new complex<double> [2 * nd + nu];
  complex<double> * mr = mm + nd;
  complex<double> * mu = mr + nd;
  Blft .expand (mm);
  U .expand (mu);
  for (jb = 0; jb < nb; jb++) {
    lambda1 = jb % sitestates;
    jb1     = jb / sitestates; 
    for (ia = 0; ia < na; ia++) {
      sigma1 = ia / na1;
      ia1    = ia % na1;
      sl1    = sigma1 + lambda1 * sitestates;
      mr [ia + jb * na] = 0.0;
      for (lambda = 0; lambda < sitestates; lambda++) {
	jm = lambda + jb1 * sitestates;
	for (sigma = 0; sigma < sitestates; sigma++) {
	  sl = sigma + lambda * sitestates;
	  im = ia1 + sigma * na1;
	  //
	  // AAAAA  Controllare se serve un segno statistico per fermioni
	  //
	  complex<double> uu = mu [ sl1 + sl * sitestates * sitestates];
	  complex<double> psi = mm [im + jm * na];
	  mr [ia + jb * na] +=  uu * psi;
	}
      }
    }
  }
  Blft .compress (mr);
  Brgt = system . base ();
  Brgt .dagger ();
  multiply (Blft, Brgt, Blft);
  Brgt = universe .base ();
  Brgt .conjugate ();
  multiply (result, Blft, Brgt);
}
//
//____________________________________________________________________________
void trotter (Block & system, Block & universe)
{
  Aarray phi;
  size_t ket,l;
  size_t sites = system .sites () + universe .sites ();
  cout << "===========> Trotter start " << system .sites () << "+"
       << universe .sites () << " " << system .states () << "x" 
       << universe .states () << endl;
  size_t oddsites = sites % 2;
  //
  if (lfttimerule .size () == 0) lfttimerule = lftrule;
  if (rgttimerule .size () == 0) rgttimerule = rgtrule;
  ket	= 0;
  for (l = 1; l < super_target .subspaces (); l++) 
    if (ttarget_id == super_target .id (l)) ket = l;
  if (ket == 0) {
    cout << "Unidentified ket as initial state " 
	 << name_define(ttarget_id) << endl;
    exit(1);
  }
  if (ttarget_index < super_target .states (ket)) 
    ket = super_target .offset (ket) + ttarget_index;
  else {
    cout << "Initial state target index " << ttarget_index 
	 << " exceeds previously found targets!" << endl;
    exit(1);
  }    
  Action vket = super_state [ket];
  //vket .show ("vket");
  super_state .clear ();
  //
  //    Applies initial operator to initial state
  //
  if (start_action .size () == 0) phi [0] = vket;
  else { 
    start_action .reorder (sites);
    //start_action .show ("ini action");
    hamaction = Superaction (start_action, system, universe, true);
    //hamaction .show ("a");
    biapply (phi [0], hamaction, vket, true);
  }
  double norm2 = multiply (phi [0], phi [0]) .real();
  double norm = sqrt(norm2);
  if (norm < 1.e-08) {
    cout << "Initial evolution state null!" << endl;
    exit(1);
  }
  phi [0] *= (1.0/norm);
  phi [0] .normalize ();
  //cout << "Norm after initial operator " << norm 
  //     << " renormalized to 1" << endl;
  norm2 = norm = 1.0;
  size_t  lft_sites = system .sites ();
  size_t  rgt_sites = sites - lft_sites;
  double  E0        = 0.0;
  size_t  tsteps    =  0;
  if (time_zips < 1) time_zips = 1;
  size_t  zips      = 0;
  long	  zipdirini = -1;
  long    zipdir    = zipdirini;
  double  t, te, tau;
  double  norm0;
  complex<double> survival;
  truncation = 0.0;
  show_states 	 = false;
  show_selection = false;
  show_density	 = 0;
  //
  super_weight .clear ();
  super_weight .push_back (1.0);
  tensorweight = 0.0;
  t  = te = 0.0;
  tau = timestep / time_zips;
  cout << setw (8)  << "time"
       << setw (15) << "<V(t)|V(t)>"
       << setw (15) << "<V(0)|V(0)>"
       << setw (15) << "Zip error"
       << setw (15) << "sites (states)"
       << endl;
  cout << setw (8)  << setprecision (4) << fixed << t
       << setw (15) << setprecision (10) << norm2
       << setw (15) << setprecision (10) << norm2
       << " zip" << setw (3) << zips
       << resetiosflags (ios_base::fixed) << right << showpoint
       << setw (10) << setprecision (3)
       << truncation << right << " " << lft_sites << "+" << rgt_sites
       << " (" << system .states () << "x" << universe .states () << ")"
       << endl;
  phi [2] = phi [0];
  properties (system, universe, phi [0], phi [0], norm2, t);
  while (true) {
    //
    //	check for new sweep
    //
    if ((lft_sites == rgt_sites + oddsites) && (zipdir == zipdirini)) {
      zips++;
      t = te;
      te = t + tau;
    }
    if (zips <= time_zips) {
		if (norm2 < 1.e-06) {
			cout << "Null evolving state!" << endl;
			return;
		}

		if (t < tau) {
			hamt_parse (sites, t);
			hamiltonian += (-E0);
		}
			//
			// AAAAA  TO DO
			// Aggiungere evoluzione dei fattori costanti
			//
			phi [1] = Action ();
			Action U;
			size_t slft = system .sites () - 1;
			// if ((tsteps <= 0) && (zips == 1)) phi [0] .show ("phi");
			//
			if (rgt_sites == 2) {
				U = evolve (hamiltonian, sites - 2, sites - 1, system, universe, tau);
				trotterrgt (phi [1], U, phi [0], system, universe);
				phi [0] = phi [1];
				/*
				if ((tsteps <= 0) && (zips == 1)) {
				  phi [0] .clean ();
				  cout <<sites-2 << " "<< sites-1<< "R "<< zipdir<< " ";
				  phi [0] .show ("Rphi");
				}
				*/
			}      
			if (((zipdir > 0) && (slft % 2 == 0)) ||
			((zipdir < 0) && (slft % 2 == 1)) ||
			(lft_sites == 2) || (rgt_sites == 2)) {
				U = evolve (hamiltonian, slft, slft + 1, system, universe, tau);
				trotterstep (phi [1], U, phi [0], system, universe);
				/*
				if ((tsteps <= 0) && (zips == 1)) {
				  phi [1] .clean ();
				  cout << slft << " " << slft+1 << "  "<<zipdir << " ";
				  phi [1] .show ("Uphi");
				}
				*/
			}
      		else phi [1] = phi [0];
			if (lft_sites == 2) {
				phi [0] = phi [1];
				U = evolve (hamiltonian, 0, 1, system, universe, tau);
				trotterlft (phi [1], U, phi [0], system, universe);
				/*
				if ((tsteps <= 0) && (zips == 1)) {
				  phi [1] .clean ();
				  cout << 0 << " " << 1 << "L " << zipdir << " ";
				  phi [1] .show ("Lphi");
				}
				*/
			}
//
			super_state [0] = phi [1];
//
//	AAA:	the actions arrays (operators) of the 
//		actual system and universe blocks
//		are transferred to oldsystem and 
//		olduniverse blocks (to save memory space)
//		system and universe action lists are now empty
//		(the space structure is the same)
//
      Block oldsystem, olduniverse;
      oldsystem    = system;
      system .action (idop_base,0) = oldsystem .base ();
      olduniverse  = universe;
      universe .action (idop_base,0) = olduniverse .base ();

      size_t actual_zip = zips;
      if (zips == 0) actual_zip = 1;
      min_lft = max_lft = min_rgt = max_rgt = 0;
      //
      for (size_t rule = 0; rule < lfttimerule .size (); rule++) {
	double trule = lfttimerule [rule] .br_time;
	size_t z = lfttimerule [rule] .br_zip;
	size_t s = lfttimerule [rule] .br_sites;
	if ((t >= trule) && (actual_zip >= z) && (system .sites () >= s)) {
	  if (min_lft < lfttimerule [rule] .br_min) 
	    min_lft = lfttimerule [rule] .br_min;
	  if (max_lft < lfttimerule [rule] .br_max) 
	    max_lft = lfttimerule [rule] .br_max;
		n_cutlft = lfttimerule [rule] .br_cut;
	}
      }
      for (size_t rule = 0; rule < rgttimerule .size (); rule++) {
	double trule = rgttimerule [rule] .br_time;
	size_t z = rgttimerule [rule] .br_zip;
	size_t s = rgttimerule [rule] .br_sites;
	if ((t >= trule) && (actual_zip >= z) && (universe .sites () >= s)) {
	  if (min_rgt < rgttimerule [rule] .br_min) 
	    min_rgt = rgttimerule [rule] .br_min;
	  if (max_rgt < rgttimerule [rule] .br_max) 
	    max_rgt = rgttimerule [rule] .br_max;	  
		n_cutrgt = rgttimerule [rule] .br_cut;
	}
      }
      updatebase (system, universe, zipdir);
      if (zipdir < 0) {
	universe .reflectreset ();
	blockrgt [rgt_sites]  .actionclear (0, name_action ());
	blockrgt [rgt_sites] = universe;
	//
	//	define new system and universe (for next step)
	//
	system   = Block (blocklft [lft_sites-2], blocklft [1]);
	if (reflect_universe)
	  universe = Block (blockrgt [1], blockrgt [rgt_sites]);
	else
	  universe = Block (blockrgt [rgt_sites], blockrgt [1]);
      } else {
	//
	//	Update list of blocks (with action transfer)
	//
	blocklft [lft_sites] .actionclear (0, name_action ());
	blocklft [lft_sites] = system;
 	//
 	//	define new system and universe (for next step)
	//
	system   = Block (blocklft [lft_sites], blocklft [1]);
	if (reflect_universe)
	  universe = Block (blockrgt [1], blockrgt [rgt_sites-2]);
	else
	  universe = Block (blockrgt [rgt_sites-2], blockrgt [1]);
      }
      //
      lft_sites += zipdir;
      rgt_sites -= zipdir;
      if ((lft_sites <= 2) || (rgt_sites <= 2)) zipdir = -zipdir;
      //
      phi [0] = Action(system .quantum (), universe .quantum ());
      projection (phi [0], system, universe, phi [1], oldsystem, olduniverse);
      norm2 = multiply (phi [0], phi [0]) .real ();
      norm = sqrt(norm2);
      phi [1] = Action(system .quantum (), universe .quantum ());
      projection (phi [1], system, universe, phi [2], oldsystem, olduniverse); 
      phi [2] = phi [1];
      norm0 = multiply (phi [2], phi [2]) .real ();
      //
      cout << setw (8)  << setprecision (4) << fixed << t
	   << setw (15) << setprecision (10) << norm2
	   << setw (15) << setprecision (10) << norm0
	   << " zip" << setw (3) << zips
	   << "  "   << resetiosflags (ios_base::fixed) << right << showpoint
	   << setw (8) << setprecision (3)
	   << truncation << right << " " << lft_sites << "+" << rgt_sites
	   << " (" << system .states () << "x" << universe .states () << ")"
	   << endl; 
      //
    } // if (zips <= ...
    else { // zips > zipmax:  advance timestep
    //properties every 10 tstep
//    if (tsteps % 10 == 0) {

      double tprop;
      tprop = fmod(t+tau/10,10);
      cout << "Entropy " << setw (14) << setprecision (10) << tprop << "   " << tau/10 << endl;   
      if ( fabs(tprop) <= tau/10) {
        properties (system, universe, phi [0], phi [0], norm2, t); 
        hamiltonian .show ("Hamiltonian polinomial:");
      }
      double ent = phi [0] .entropy ();
      cout << "Entropy " << setw (14) << setprecision (10) << t << "   " << ent << endl;   
      tsteps++;

      if (tsteps >= steps_number) break;
      zips   =  0;
      zipdir =	zipdirini;
      te = t;
    }
  } // while (true)
}
//
//============================================================================
void evolution (Block & system, Block & universe)
{
    Aarray phi;
    Action phihalf;
    size_t ket,l;
    double rungekutta [] =
    { 1.0, 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
    double onethird   [] =
    { 1.0, 31.0/162.0, 14.0/162.0, 14.0/162.0, -5.0/162.0 };
    double twothird   [] =
    { 1.0, 16.0/81.0, 20.0/81.0, 20.0/81.0, -2.0/81.0 };
    double weight     [] =
    { 1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/3.0 };
    //
    size_t sites = system .sites () + universe .sites ();
    cout << "===========> Evolution start " << system .sites () << "+"
    << universe .sites () << " " << system .states () << "x"
    << universe .states () << endl;
    size_t oddsites = sites % 2;
    //
    if (lfttimerule .size () == 0) lfttimerule = lftrule;
    if (rgttimerule .size () == 0) rgttimerule = rgtrule;
    //
    //    Applies initial operator to initial state
    //
    ket	= 0;
    for (l = 1; l < super_target .subspaces (); ++l)
        if (ttarget_id == super_target .id (l)) ket = l;
    if (ket == 0) {
        cout << "Unidentified ket as initial state "
        << name_define(ttarget_id) << endl;
        exit(1);
    }
    if (ttarget_index < super_target .states (ket))
        ket = super_target .offset (ket) + ttarget_index;
    else {
        cout << "Initial state target index " << ttarget_index
        << " exceeds previously found targets!" << endl;
        exit(1);
    }
    Action vket = super_state [ket];
    super_state .clear ();
    if (start_action .size () == 0) phi [0] = vket;
    else {
        start_action .reorder (sites);
        hamaction = Superaction (start_action, system, universe, true);
        biapply (phi [0], hamaction, vket, true);
    }
    super_state .clear ();
    double norm2 = multiply (phi [0], phi [0]) .real();
    double norm = sqrt(norm2);
    if (norm < 1.e-08) {
        cout << "Initial evolution state null!" << endl;
        exit(1);
    }
    phi [0] *= (1.0/norm);
    phi [0] .normalize ();
    cout << "Norm after initial operator " << norm
    << " renormalized to 1" << endl;
    //
    size_t  lft_sites = system .sites ();
    size_t  rgt_sites = sites - lft_sites;
    double  t         = 0.0;
    double  E0        = 0.0;
    size_t  tsteps    = 0;
    size_t  zips      = 0;
    long	  zipdirini = -1;
    long    zipdir    = zipdirini;
    double  tau       = timestep;
    double  tm, te;
    truncation = 0.0;
    show_states 	 = false;
    show_selection = false;
    show_density	 = 0;
    //
    super_weight .clear ();
    for (size_t nw = 0; nw < 4; nw++) super_weight .push_back (weight[nw]);
    phi [6] = phi [0];

    complex<double> phase;
 
    //
    hamt_parse (sites, t);
    hamiltonian += (-E0);
    hamaction = Superaction (hamiltonian, system, universe, true);
    biapply (phi [1], hamaction, phi [0], true);
    phi[1] *= complex<double> (0.,-tau);
    //
    phase = multiply (phi [0], phi [1]) * (complex<double> (0.,-1.));
    double Hsquare = multiply (phi [1], phi [1]) .real ();
    double Esquare = (phase * phase) .real ();
    double variance = Hsquare - Esquare;
    
    cout << "variance " << variance << " HH " << Hsquare << " EE " << Esquare << endl;
    
    while (true) {
        //
        //	check for new sweep
        //
        if ((lft_sites == rgt_sites + oddsites) && (zipdir == zipdirini)) zips++;
        norm2 = multiply (phi [0], phi [0]) .real();
        norm = sqrt(norm2);
        cout << "<phi|phi>( " << setw (10) << setprecision (6) << fixed << t
        << " ) = " << setw (12) << setprecision (10) << norm2
        << " truncation " << resetiosflags (ios_base::fixed)
        << left << showpoint << setw (8) << setprecision (3)
        << truncation << right
        << " zip " << zips << " " << lft_sites << "+" << rgt_sites
        << " " << system .states () << "x" << universe .states ()
        << endl;
        truncation = 0.0;
        if (norm2 < 1.e-06) {
            cout << "Null evolving state!" << endl;
            return;
        }
        
        if (zips <= time_zips) {
            //
            super_state [0] = phi [0];
            for (size_t i = 1; i < 5; ++i) phi [i] = Action ();
            tau = timestep;
            tm = t + tau /2.0;
            te = t + tau;
            //
            hamt_parse (sites, t);
            hamiltonian += (-E0);
            hamaction = Superaction (hamiltonian, system, universe, true);
            biapply (phi [1], hamaction, phi [0], true);
            phi[1] *= complex<double> (0.,-tau);
            //
            hamt_parse  (sites, tm);
            hamiltonian += (-E0);
            hamaction = Superaction (hamiltonian, system, universe, true);
            phihalf = phi[1];
            phihalf *= 0.5;
            phihalf += phi[0]; // phi [0] + 1/2 phi [1] { ~= phi (t + 1/2 tau) }
            phihalf .normalize ();
            biapply (phi[2], hamaction, phihalf, true);
            phi [2] *= complex<double>(0,-tau);
            //
            phihalf  = phi [2];
            phihalf *= 0.5;
            phihalf  += phi [0];
            phihalf .normalize ();
            biapply (phi [3], hamaction, phihalf, true);
            phi [3] *= complex<double> (0, -tau);
            //
            hamt_parse  (sites, te);
            hamiltonian += (-E0);
            hamaction = Superaction (hamiltonian, system, universe, true);
            phihalf =  phi [3];
            phihalf += phi [0];
            phihalf .normalize ();
            biapply (phi [4], hamaction, phihalf, true);
            phi [4]  *= complex<double> (0,-tau);
            //
            phi [5] = phi [0];
            for (size_t i = 1;  i < 5; ++i) {
                phihalf = phi [i];
                phihalf *= rungekutta [i];
                phi [5] += phihalf;
            }
            //
            Action phionethird = phi [0];
            for (size_t i = 1;  i < 5; ++i) {
                phihalf = phi [i];
                phihalf *= onethird [i];
                phionethird += phihalf;
            }
            super_state [1] = phionethird;
            //
            Action phitwothird = phi [0];
            for (size_t i = 1;  i < 5; ++i) {
                phihalf = phi [i];
                phihalf *= twothird [i];
                phitwothird += phihalf;
            }
            super_state [2] = phitwothird;
            //
            //	AAA:	the actions arrays (operators) of the
            //		actual system and universe blocks
            //		are transferred to oldsystem and
            //		olduniverse blocks (to save memory space)
            //		system and universe action lists are now empty
            //		(the space structure is the same)
            //
            Block oldsystem    = system;
            Block olduniverse  = universe;
            //
            super_state [3] = phi [5];
            phi [1] = phionethird;
            phi [2] = phitwothird;
            phi [3] = phi [5];
            //
            //
            //	Set the number of wanted dmrg states according to input rules.
            //	A rule specify the minimum and maximum number of DMRG states
            //	starting from a given sweep (zip) and number of sites.
            //
            size_t actual_zip = zips;
            if (zips == 0) actual_zip = 1;
            min_lft = max_lft = min_rgt = max_rgt = 0;
            //
            for (size_t rule = 0; rule < lfttimerule .size (); rule++) {
                double trule = lfttimerule [rule] .br_time;
                size_t z = lfttimerule [rule] .br_zip;
                size_t s = lfttimerule [rule] .br_sites;
                if ((t >= trule) && (actual_zip >= z) && (system .sites () >= s)) {
                    if (min_lft < lfttimerule [rule] .br_min)
                        min_lft = lfttimerule [rule] .br_min;
                    if (max_lft < lfttimerule [rule] .br_max)
                        max_lft = lfttimerule [rule] .br_max;
                    n_cutlft = lfttimerule [rule] .br_cut;
                }
            }
            for (size_t rule = 0; rule < rgttimerule .size (); rule++) {
                double trule = rgttimerule [rule] .br_time;
                size_t z = rgttimerule [rule] .br_zip;
                size_t s = rgttimerule [rule] .br_sites;
                if ((t >= trule) && (actual_zip >= z) && (universe .sites () >= s)) {
                    if (min_rgt < rgttimerule [rule] .br_min)
                        min_rgt = rgttimerule [rule] .br_min;
                    if (max_rgt < rgttimerule [rule] .br_max)
                        max_rgt = rgttimerule [rule] .br_max;
                    n_cutrgt = rgttimerule [rule] .br_cut;
                }
            }
            updatebase (system, universe, zipdir);
            
            if (zips == time_zips) {
            	double ent = phi [0] .entropy ();
            	cout << "Entropy_inner " << setw (10) << setprecision (6) << fixed
                << t << setw (10) << setprecision (6) << fixed
                << system .sites () << " " << setw (14) << setprecision (10) << ent << endl; }               
            
            
            if (zipdir < 0) {
                universe .reflectreset ();
                blockrgt [rgt_sites]  .actionclear (0, name_action ());
                blockrgt [rgt_sites] = universe;
                //
                //	define new system and universe (for next step)
                //
                system   = Block (blocklft [lft_sites-2], blocklft [1]);
                universe = Block (blockrgt [rgt_sites],   blockrgt [1]);
            } else {
                //
                //	Update list of blocks (with action transfer)
                //
                blocklft [lft_sites] .actionclear (0, name_action ());;
                blocklft [lft_sites] = system;
                //
                //	define new system and universe (for next step)
                //
                system   = Block (blocklft [lft_sites],   blocklft [1]);
                universe = Block (blockrgt [rgt_sites-2], blockrgt [1]);
            }
            /*
             if (zipdir < 0) {
             Action density (universe .quantum (), 0.);
             for (size_t i = 0; i < 4; ++i) {
             phihalf = phi [i];
             phihalf .dagger ();
             phihalf *= phi [i];
             phihalf *= weight [i];
             density += phihalf;
             }
             density .transpose ();
             universe .select_states (density, min_rgt, max_rgt, n_cut);
             //
             //	Update list of blocks (with action transfer)
             //
             universe .reflectreset ();
             blockrgt [rgt_sites]  .actionclear (0, name_action ());
             blockrgt [rgt_sites] = universe;
             //
             //	define new system and universe (for next step)
             //
             system   = Block (blocklft [lft_sites-2], blocklft [1]);
             universe = Block (blockrgt [rgt_sites],   blockrgt [1]);
             } else {
             Action density (system .quantum (), 0.);
             for (size_t i = 0; i < 4; ++i) {
             phihalf     = phi [i];
             phionethird = phi [i];
             phionethird .dagger ();
             phihalf *= phionethird;
             phihalf *= weight [i];
             density += phihalf;
             }
             system .select_states (density, min_lft, max_lft, n_cut);
             //
             //	Update list of blocks (with action transfer)
             //
             blocklft [lft_sites] .actionclear (0, name_action ());;
             blocklft [lft_sites] = system;
             //
             //	define new system and universe (for next step)
             //
             system   = Block (blocklft [lft_sites],   blocklft [1]);
             universe = Block (blockrgt [rgt_sites-2], blockrgt [1]);
             }
             */
            //
            //	Reflection
            //
            if (reflect_blocks) {
                universe .reflectlft ();
                universe .reflectrgt ();
            }
            if (reflect_universe) universe .reflect ();
            //
            phi [4] = Action(system .quantum (), universe .quantum ());
            projection (phi [4], system, universe, phi [0], oldsystem, olduniverse);
            //
            phi [0] = phi [4];
            phi [1] = Action(system .quantum (), universe .quantum ());
            projection (phi [1], system, universe, phi [6], oldsystem, olduniverse);
            //
            phi [6] = phi [1];
            //
            lft_sites += zipdir;
            rgt_sites -= zipdir;
            if ((lft_sites <= 2) || (rgt_sites <= 2)) {
                zipdir = -zipdir;
            
            }
        }
        else {
            //
            // Runge Kutta fine steps
            //
            tau = timestep / step_divisions;
            properties (system, universe, phi [0], phi [0], norm2, t);
            double norm0 = multiply (phi [6], phi [6]) .real ();
            double ent = phi [0] .entropy ();
            cout << "Entropy " << setw (14) << setprecision (10) << ent << endl;
            
            for (size_t l = 1; l < renyientropy .size (); ++l)
                cout << setw (39) << right
                << "Renyi entropy (" << setw (19) << setprecision (6) << fixed
                << renyiarg [l] << ")" <<  setw (19) << setprecision (12)
                << fixed << renyientropy [l] << endl;
            
            for (size_t division = 0; division < step_divisions; division++) {
                for (size_t i = 1; i < 5; ++i) phi [i] = Action ();
                tm = t + tau/2.0;
                te = t + tau;
                hamt_parse (sites, t);
                hamiltonian += (-E0);
                hamaction = Superaction (hamiltonian, system, universe, true);
                biapply (phi [1], hamaction, phi [0], true);
                //
                norm2 = multiply (phi [0], phi [0]) .real ();
                phase = multiply (phi [0], phi [1]) * (complex<double> (0.,-1.));
//                double Hsquare = multiply (phi [1], phi [1]) .real ();
//                double Esquare = (phase * phase) .real ();
//                double variance = Hsquare - Esquare;
                
                complex<double> survival = multiply (phi [6], phi [0]);
                cout << "<phi|phi>( " << setw (10) << setprecision (6) << fixed
                << t
                << " ) = " << setw (10) << setprecision (8) << norm2
                << " <phi|d phi/dt> = " << setw (26) << setprecision (8)
                << phase
                << resetiosflags (ios_base::fixed) << endl;
                cout << "Survival (" << setw (10) << setprecision (6) << fixed
                << t
                << ") = " << setw (10) << setprecision (8)
                << sqrt( (survival .real()) * (survival .real ()) +
                        (survival .imag()) * (survival .imag ()))
                << " <phi0|phi0>   = " << setw (26) << setprecision (8)
                << norm0
                << resetiosflags (ios_base::fixed) << endl;
                phi[1] *= complex<double> (0.,-tau);
                //
                hamt_parse  (sites, tm);
                hamiltonian += (-E0);
                hamaction = Superaction (hamiltonian, system, universe, true);
                phihalf = phi[1];
                phihalf *= 0.5;
                phihalf += phi[0];
                phihalf .normalize ();
                biapply (phi[2], hamaction, phihalf, true);
                phi [2] *= complex<double>(0,-tau);
                //
                phihalf  = phi [2];
                phihalf *= 0.5;
                phihalf  += phi [0];
                phihalf .normalize ();
                biapply (phi [3], hamaction, phihalf, true);
                phi [3] *= complex<double> (0, -tau);
                //
                hamt_parse  (sites, te);
                hamiltonian += (-E0);
                hamaction = Superaction (hamiltonian, system, universe, true);
                phihalf =  phi [3];
                phihalf += phi [0];
                phihalf .normalize ();
                biapply (phi [4], hamaction, phihalf, true);
                phi [4]  *= complex<double> (0,-tau);
                //
                phi [5] = phi [0];
                for (size_t i = 1;  i < 5; ++i) {
                    phihalf = phi [i];
                    phihalf *= rungekutta [i];
                    phi [5] += phihalf;
                }
                t = te;
                phi [0]  = phi [5];
                //properties (system, universe, phi [0], phi [0]);
            }
            tsteps++;
            if (tsteps >= steps_number) break;
            zips   =  0;
            zipdir = -1;
        }
    } // while (true)
    norm2 = multiply (phi [0], phi [0]) .real();
    norm = sqrt(norm2);
    cout << "<phi|phi>( " << setw (10) << setprecision (6) << fixed << t
    << " ) = " << setw (12) << setprecision (10) << norm2
    << " truncation " << resetiosflags (ios_base::fixed)
    << left << showpoint << setw (8) << setprecision (3)
    << truncation << right
    << " zip " << zips << " " << lft_sites << "+" << rgt_sites
    << " " << system .states () << "x" << universe .states ()
    << endl;
    properties (system, universe, phi [0], phi [0], norm2, t);
}
//
//____________________________________________________________________________
void hamiltonian_energies (Block & system, Block & universe, Qtarget & target)
{
  //
  //	Compute wanted eigenstates and eigenvalues of the hamiltonian
  //
  Aarray super;
  size_t sites = system .sites () + universe .sites ();
  //
  //	Update hamiltonian expression (if needed)
  //
  if (sites_old != sites) {
    ham_parse (sites);
    if (sites <= show_hamiltonian || showham)  
      hamiltonian .show ("Hamiltonian polinomial:");
  }
  //
  //	Trasform the hamiltonian polinomial into a sum of
  //	tensor products of operators acting on system and 
  //	universe.  
  //
  hamaction = Superaction (hamiltonian, system, universe, true);
  if (sites <= show_hamiltonian || showham) {
    cout << "Splitted hamiltonian (at " << system .sites () << " sites):" 
	 << endl;
    hamaction .show ();
  }
  showham = 0;
  purelft = hamaction .onlylft ();
  purergt = hamaction .onlyrgt ();
  bool complexham = hamaction .iscomplex ();
  if (purelft >=0) {
    if (energylft) delete [] energylft;
    hlft = hamaction [purelft] .ba_lftop;
    edimlft = hlft .width ();
    energylft = new double [edimlft];
    hlft .eigen (hlft, energylft);
    /*
    cout << setw (16) << left << "elft: " << right;
    for (size_t i = 0; i < edimlft; i++)
      cout << setw (14) << setprecision (8) << fixed << energylft [i] << "  ";
    cout << endl;
    hlft .show ("ulft");
    */
  }
  if (purergt >=0) {
    if (energyrgt) delete [] energyrgt;
    hrgt = hamaction [purergt] .ba_rgtop;
    edimrgt = hrgt .width ();
    energyrgt = new double [edimrgt];
    hrgt .eigen (hrgt, energyrgt);
    /*
    cout << setw (16) << left << "ergt:" << right;
    for (size_t i = 0; i < edimrgt; i++)
      cout << setw (14) << setprecision (8) << fixed << energyrgt [i] << "  ";
    cout << endl;
    */
  }
  double qa = 1.0;
  if (system .sites () % 2) qa = -1.0;
  if ((reflection & 1) && (universe .sites () % 2 == 0)) qa = -qa;
  for (size_t nq = 1; nq < target .subspaces (); nq++) {
    //
    //	Setup initial guess (and result structure)
    //
    Action state (system .quantum (), universe .quantum (), 
		  target .number (nq), complexham, qa);
    state .storage  (state .size ());
    state .scalar   (1.0);
    //
    //	Look if we can project old state as initial guess
    //
    double norm = 1.0;
    //double * m;
    initial_guess = 0;
    /*
    if (old_sites == sites) {
      size_t off = super_target .offset (nq);
      size_t old = super_target .states (nq);
      Action guess = super_state [off];
      for (size_t i = off + 1; i < off + old; i++)  guess += super_state [i];
      //
      Block oldsystem   = Block (blocklft [old_sys - 1], blocklft [1]);
      Block olduniverse = Block (blockrgt [old_uni - 1], blockrgt [1]);
      if (reflection & 2) {
	olduniverse .reflectlft ();
	olduniverse .reflectrgt ();
      }
      if (reflection & 1) olduniverse .reflect ();
      //
      Action projected (state);
      projection (projected, system, universe, guess, oldsystem, olduniverse);
      state += projected;
      m = state .storage ();
      norm = 0.0;
      for (size_t i = csize; i < state .size (); i++) norm += m [i] * m [i];
      if (norm >= zero_norm) {
	initial_guess = 1;
	double dnorm = 1.0 / sqrt (norm);
	norm = ((double) old - norm) / old;
	if (norm < 0.0) norm = 0.0;
	norm = sqrt (norm);
	for (size_t i = csize; i < state .size (); i++) m [i] *= dnorm;
      }
      else norm = 1.0;
    }
    */
    //
    //	Random initial guess
    //
    if (initial_guess == 0)  state .random   ();
    //
    //	Prepare for hamiltonian applications (no explicit releasing)
    //
    size_t oldsize   = state .size ();
    size_t olddim    = state .dimension ();
    size_t oldblocks = state .blocks ();
    (void) biaction (state, hamaction, state);
    cout << "Target " 
	 << setw (20) << left << target .number (nq) .str () << " "
	 << "states " << setw (2) << right << target .states (nq) << "/" 
         << setw (0) << left << olddim << "+" 
	 << state .dimension () - olddim << " " 
	 << "blocks " << setw (3) << right << oldblocks 
	 << "+" << state .blocks () - oldblocks << " "
	 << "statistic " << setw (2) << right << target .statistic (nq) 
	 << setw (0) << endl;
    /*
    size_t insize    = biaction (state, hamaction, state);
    double * sitec = (double *) name_value ("sitecontrol");
    if (sitec && (sites == (size_t) *sitec)) {
      cout << "insize = " << insize << " oldsize = " << oldsize << endl;
      double * mem = state .storage ();
      size_t n1 = csize-1;
      size_t n2 = state .size ();
      n1 = csize   + 1302;
      n2 = csize   + 1750;
      for (size_t i = csize; i < state.size (); i++)	{
	mem [i] = 0.0;
	if ((n1 <= i) && (i < n2))	mem [i] = 1.0;
	}
      Action res;
      Action pp;
      //res .storage (oldsize);
      multiply (pp, state, hamaction [1] .ba_rgtop, true);
      pp .show ("pp");
      biapply (res, hamaction, state);
      mem = res.storage ();
      cout << "full";
      size_t ii = 1;
      for (size_t i = csize; i < res.size (); i++) {
	if (i %20 == csize) cout << endl << setw(3) << ii++ << " ";
	cout << setw(3) << mem[i];
      }
      double al = 0.0;
      for (size_t i = csize; i < res .size (); i++) al += mem [i];
      cout << endl;
      cout << "al = " << al ;
      double bl = multiply (res, state) .real ();
      cout << " <H> = " << bl << " "
	   << state.size () - csize << " " << res .size () - csize << " ";
      res.show ("1shot");
      state .show ("from");
      continue;
    }
    */
    //
    //	Setup thick-restatrt Lanczos algorithm
    //
    //if ((system .sites () == 2) && (universe .sites () == 2)) 
    //  state .show ("Guess");
    Lanczos lanczos;
    lanczos .parameters (state .size (), oldsize,
			 target .states (nq), state .dimension (), 
			 lanczos_strategy, lanczos_iters, lanczos_repeats, 
			 zero_energy, zero_norm);
    size_t found = lanczos .thick (hamaction, state);
    //
    //	Store energies and super states
    //
    target .energy (lanczos .value (), nq, found);
    vector<double> entropy;
    vector<double> frustration;
    //
    //	Output energies, tolerances and extimated errors
    //
    cout << setw (24) << right << "Energy" 
	 << setw (10) << right << "tolerance" 
	 << setw (11) << right << "residual " 
	 << setw (11) << right << "entropy ("
	 << setw (4)  << right << system .sites ()
	 << setw (1)  << ")"
	 << setw (12) << right << "GLIO ("
	 << setw (4) << right << system .sites () << ")" << endl; 
    if (initial_guess) 
      cout << setw (18) << right << "Initial guess  = " 
	   << setw (19) << setprecision (12) << fixed 
	   << right << initial_value << " ("
	   << setw (9) << setprecision (2)  << scientific << norm << " " 
	   << setw (9) << setprecision (2)  << scientific 
	   << initial_residual << ")" << endl;
    for (size_t j = 0; j < lanczos .steps (); j++) {
      if (j < found) {
	state .storage (lanczos .storage (j));
	if ((system .sites () == 2) && (universe .sites () == 2)) 
	  state .show ("State");
	super [target .offset (nq) + j] << state;
	double ent = state .entropy ();
	entropy .push_back (ent);
	frustration .push_back (frustlft);
	state .release ();
      } 
      bool accept, trust;
      double v = lanczos .value     (j);
      double t = lanczos .tolerance (j);
      double e = lanczos .error     (j);
      if (t < 0.0) t = -t;
      if (e < 0.0) e = -e;
      accept = (t < 10.0 * zero_energy) || (e < 10.0 * zero_norm);
      trust  = (t < zero_energy) || (e < zero_norm);
      ofstream file_entropy("entropy.dat",ios::out | ios::app);
      if ((j < found) || accept) {
	//if (trust) 	cout << setw (8) << right << "Trusted ";
	//else	        cout << setw (8) << right << " ";
	cout << "E[" << setw (3) << j << "]= "  
	     << setw (16) << setprecision (10) << fixed << v << " "
	     << setw (9) << setprecision (2)  << scientific << t << " " 
	     << setw (9) << setprecision (2)  << scientific << e << " ";
	if (j < found) { 
	  cout << setw (16) << setprecision (10) << fixed
	       << entropy [j] << " " << setw (16) << setprecision (10)
	       << fixed << frustration [j];
	  file_entropy << setw(3) << system .sites () <<" + "
		       << universe.sites()  << setw (19) << setprecision (12) 
		       << fixed<< entropy [j] <<endl;
	}
	cout << endl;
      }
      file_entropy.close(); 
      if (j < found) {
	if ((renyiarg.size() > 1) && (renyiarg [1] < 0.0)) { 
	  //
	  // Stampa gli autovalori richiesti sul file eigenvalues.dat
	  //
	  ofstream file_eigenvalues("eigenvalues.dat",ios::out | ios::trunc);
	  file_eigenvalues << "system+universe = "<< system.sites ()
			   <<"+"<< universe.sites()<<endl;
	  cout << endl;
	  for (int k = 0; k< ((int)abs(renyiarg[1])); k++) {
	    cout << "eigenvalue[ "<< k << "]= " << scientific 
		 <<autovalori[k]<<endl;
	    file_eigenvalues<< scientific <<autovalori[k]<<endl;
	  }
	  cout << endl;
	  // file_eigenvalues << endl;
	  file_eigenvalues.close();   
	}
	for (size_t l = 1; l < renyientropy .size (); l++)     
	  if (renyiarg [l]>0)
	    cout << setw (39) << right << "Renyi entropy (" 
		 << setw (19) << setprecision (6) << fixed 
		 << renyiarg [l] << ")" <<  setw (19) << setprecision (12) 
		 << fixed << renyientropy [l] << endl;
      }
    }
  }
  //
  //	Remember states and target for next dmrg step 
  //	
  sites_old    = sites;
  super_state  = super;
  super_weight .clear ();
  for (size_t n = 0; n < super_state .size (); n++) 
    super_weight .push_back (1.0);
  super_target = target;
  super_state .release ();  // release memory
  //
  //	reset output options
  //
  cout .setf ((ios_base::fmtflags) 0, ios_base::floatfield);
  cout << setprecision (6);
}
//
//____________________________________________________________________________
static void projection (Action & newstate, 
			Block & newsystem, Block & newuniverse,
			const Action & oldstate,
			Block & oldsystem, Block & olduniverse)
{
  //
  // 	Build projection of oldstate to newstate
  //
  // 	Given a state oldstate on a decomposition (n,m) of n+m sites
  //	(oldsystem,olduniverse):
  //
  //	|psi> = \sum |a_n b_m> <a_n b_m | psi>
  //
  //	compute the projection in a new decomposition (n+1,m-1) or
  //	(n-1,m+1) of n+m sites (newsystem,mewuniverse):
  //
  //	P| psi> = \sum | a_{n+1} b_{m-1}> < a_{n+1} b{m-1} | psi>
  //
  //	or
  //
  //	P| psi> = \sum | a_{n-1} b_{m+1}> < a_{n-1} b{m+1} | psi>
  //
  size_t n  = oldsystem .sites ();
  size_t nn = newsystem .sites ();
  size_t m  = olduniverse .sites ();
  size_t mm = newuniverse .sites ();
  if ((n + m) != (nn + mm)) {
    cout << "Invalid projection!" << endl;
    exit (0);
  }
  //
  Action a, b, c;
  if (nn > n) {
    //
    //	sweep to right
    //
    if (reflect_universe)  {
      //  <a_{n+1}; b_{m-1} | psi> = \sum <a_{n+1} | s_n p> <s_n | a_n>
      //	<b_{m-1} | t_{m-1}> <p t_{m-1} | b_m> <a_n; b_m |psi>
      a = blocklft [n] .base ();
      a .dagger ();	    		// a = <s_n | s_{n-1} p>
      b = oldsystem .base ();		// b = <s_{n-1} p | a_n>
      multiply (c, a, b);	  	// c = <s_n | a_n>
      multiply (a, c, oldstate);  	// a = <s_n; b_m |psi>
      b = olduniverse .base ();		// b = <p t_{m-1} |b_m>
      b .transpose ();
      multiply (newstate, a, b); // newstate = <s_n; p t_{m-1} |psi>
      b = newuniverse .base ();
      b .dagger ();			// b = <b_{m-1} | q t_{m-2}>
      a = blockrgt [m-1] .base ();	// a = <q t_{m-2} | t_{m-1}>
      multiply (c, b, a);		// c = <b_{m-1} | t_{m-1}>
      a = newsystem .base ();
      a .dagger ();			// a = <a_{n+1} | s_n p>
      //
      //	move site index to left
      //
      size_t nd = newstate .height () * newstate .width ();
      complex<double> * mm = new complex<double> [nd];
      newstate .expand (mm);
      newstate = Action (a .domain (), c .domain ());
      newstate .compress (mm);	 // newstate = <s_n p; t_{m-1} |psi>
      delete [] mm;
      //
      c .transpose ();
      multiply (b, newstate, c);	// b = <s_n p; b_{m-1} |psi>
      multiply (newstate, a, b); // newstate = <a_{n+1}; b_{m-1} |psi>
    }
    else {
      cout << "projection with reflect_universe false not implemented!"<< endl;
      exit (0);
    }
  }
  else {
    //
    //	swwep to left
    //
    if (reflect_universe) {
      a = blockrgt [m] .base ();
      a .dagger ();			// a = <t_m | q t_{m-1}>
      b = olduniverse .base ();		// b = <q t_{m-1} | b_m>
      multiply (c, a, b);		// c = <t_m | b_m>
      c .transpose ();
      multiply (a, oldstate, c);	// a = <a_n; t_m | psi>
      b = oldsystem .base ();		// b = <s_{n-1} p | a_n>
      multiply (newstate, b, a); // newstate = <s_{n-1} p; t_m | psi>
      a = newsystem .base ();
      a .dagger ();		 	// a = <a_{n-1} | s_{n-2} p>
      b = blocklft [n-1] .base ();	// b = <s_{n-2} p | s_{n-1}>
      multiply (c, a, b);		// c = <a_{n-1} | s_{n-1}>
      b = newuniverse .base ();		// b = <p t_m | b_{m+1}>
      //
      //  	move site index to right
      //
      size_t nd = newstate .height () * newstate .width ();
      complex<double> * mm = new complex<double> [nd];
      newstate .expand (mm);
      newstate = Action (c .domain (), b .range ());
      newstate .compress (mm);	 // newstate = <s_{n-1}; p t_m |psi>
      delete [] mm;
      b .conjugate ();
      multiply (a, c, newstate);
      multiply (newstate, a, b);
    }
    else {
      cout << "projection with reflect_universe false not implemented!"<< endl;
      exit (0);
    }
  }  
}
//
//____________________________________________________________________________
void projectionold (Action & newstate, 
			   Block & newsystem, Block & newuniverse,
			   const Action & oldstate,
			   Block & oldsystem, Block & olduniverse)
{
  //
  // Build projection of oldstate to newstate
  //
  // newstate [ll pp; rr qq] = < ll pp R(rr qq) | l p R(r q) > 
  //				        oldstate [l p ; r q] 
  //
  // Sweep to right:
  //
  // sites (ll) > sites (l), reflection=0, R(k t)= k t:
  //  newstate [ll pp; rr qq] = < ll | l p > < pp rr | r > 
  //				< qq | q > oldstate [l p ; r q]
  //
  // sites (ll) > sites (l), reflection=1, R(k t)= t R(k)
  //						   (-1)^F(k,t):
  //  newstate [ll pp; rr qq] = < ll | l p > < qq R(rr) | R(r) >
  //     < pp | q > (-1)^F(rr qq) (-1)^F(r q) oldstate [l p ; r q]
  //   = < ll | l p > < pp | q > < rr qq | r >(-1)^F(r q) 
  //					oldstate [l p ; r q]
  //
  // sites (ll) > sites (l), reflection=2, R(k t)= R(k) t:
  //  newstate [ll pp; rr qq] = < ll | l p > < pp R(rr) | R(r) >
  //			 < qq | q > oldstate [l p ; r q]
  //   = < ll | l p >  < rr pp | r > < qq | q > 
  //		       (-1)^F(rr pp) oldstate [l p; r q]
  //
  // sites (ll) > sites (l), reflection=3, R(k t)= t k (-1)^F(k t):
  //  newstate [ll pp; rr qq] = < ll | l p > < qq rr | r >
  //      < pp | q >  (-1)^F(rr qq) (-1)^F(r q) oldstate [l p; r q]
  //
  // Sweep to left:
  //
  // sites (ll) < sites (l), reflection=0, R(k t)= k t:
  //  newstate [ll pp; rr qq] = < ll pp | l > < rr | p r >
  //			 < qq | q >  oldstate [l p ; r q]
  //
  // sites (ll) < sites (l), reflection=1, R(k t)= t R(k) 
  //						   (-1)^F(k,t):
  //  newstate [ll pp; rr qq] = < ll pp | l > < R(rr) | q R(r) > 
  //     < qq | p > (-1)^F(rr qq) (-1)^F(r q) oldstate [l p ; r q]
  //  = < ll pp | l > < qq | p > < rr | r q > (-1)^F(rr qq) 
  //					oldstate [l p ; r q]
  //
  // sites (ll) < sites (l), reflection=2,  R(k t)= R(k) t:
  //  newstate [ll pp; rr qq] = < ll pp | l > < R(rr) | p R(r) >
  //  				< qq | q > oldstate [l p ; r q]
  //   = < ll pp | l > < rr | r p > < qq | q > (-1)^F(r p) 
  //				oldstate [l p ; r q]
  //
  // sites (ll) < sites (l), reflection=3, R(k t)= t k 
  //						   (-1)^F(k t):
  //  newstate [ll pp; rr qq] = <ll pp | l > < qq | p > < rr | q r >
  //                   (-1)^F(rr qq) (-1)^F(r q) oldstate [l p; r q]
  //
  size_t newh = newsystem   .states ();
  size_t neww = newuniverse .states ();
  Storage newfull (newh * neww * sizeof (complex<double>));
  complex<double> * newm = (complex<double> *) newfull .storage ();
  //
  size_t pstates = blocklft [1] .states ();
  //
  size_t newlsites = newsystem   .sites () - 1;
  size_t newrsites = newuniverse .sites () - 1;
  size_t oldlsites = oldsystem   .sites () - 1;
  size_t oldrsites = olduniverse .sites () - 1;
  //
  size_t newlstates = blocklft [newlsites] .states ();
  size_t newrstates = blockrgt [newrsites] .states ();
  //size_t oldlstates = block [oldlsites] .states ();
  size_t oldrstates = blockrgt [oldrsites] .states ();
  //
  Action ubase, rbase;
  // 
  //  On the left (system side) we have always the base 
  //  components (no reflection)
  //
  if (newlsites > oldlsites) {
    //
    // sweep to right
    //
    ubase = blocklft [newlsites] .base ();
    ubase .dagger ();			// ubase [ll; l p] = < ll | l p >
  }
  else {
    //
    //   sweep to left
    //
    ubase = blocklft [oldlsites] .base ();  
    // ubase [ll pp; l] = < ll pp | l >
  }
 //
  //  On the right (universe side) we can have or not 
  //  reflected base components
  //
  if ( (reflection == 0) || (reflection == 3)) {
    //
    // No reflection or double reflection (universe and block): 
    // reflected base
    //
    if (newrsites > oldrsites) {
      //
      //  sweep to left
      //
      rbase = blockrgt [newrsites] .basereflected ();
      rbase .dagger ();                 //  rbase [rr; r q] = < rr | q r > 
    }
    else {
      //
      // sweep to right 
      //
      rbase = blockrgt [oldrsites] .basereflected (); 
      // rbase [rr pp; r] = < pp rr | r >
    }
  }
  if ((reflection == 1) || (reflection == 2)) {
    //
    //   Single reflection (universe or block)
    //
    if (newrsites > oldrsites) {
      //
      //  sweep to left
      //
      rbase = blockrgt [newrsites] .base ();  
      rbase .dagger ();                 // rbase [rr; r p] = < rr | r p >
    }
    else {
      //
      //  sweep to right
      //
      rbase = blockrgt [oldrsites] .base (); 
      // rbase [rr pp; r] = < rr pp | r >
    }
  }
  //
  if (newlsites > oldlsites) {
    //
    //   Sweep to right   
    //
    ubase *= oldstate;
    //
    // ubase [ll ; r q] = < ll | l p > * oldstate [l p; r q] 
    //
    rbase .transpose ();  
    Aarray project;
    olduniverse .contraction (ubase, pstates, 0, rbase, project, 
			      (reflection == 1) || (reflection == 3));
    //
    // project [q] [newl ; newr t] = ubase [newl; oldr q] 
    //		   * rbase [newr t; oldr] (-1)^F(oldr q)
    //
    size_t oldmsize = newlstates * newrstates*pstates;
    Storage oldfull (oldmsize * sizeof (complex<double>));
    complex<double> * oldm = (complex<double> *) oldfull .storage ();
    for (size_t q = 0; q < pstates; q++) {
      for (size_t i = 0; i < oldmsize; i++) oldm [i] = 0.0;
      project [q] .expand (oldm);
      for (size_t t = 0; t < pstates; t++) {
	//
	// define pp qq
	//
	size_t pp = t;
	size_t qq = q;
	if ((reflection == 1) || (reflection == 3)) {
	  pp = q;
	  qq = t;
	}
	size_t tsubs = blockrgt [1] .quantumsub () [t];
	long   tstat = blockrgt [1] .quantumstat (tsubs);	  
	for (size_t rnew = 0 ; rnew < newrstates; rnew++) {
	  size_t rsubs = blockrgt [newrsites] .quantumsub () [rnew];
	  long   rstat = blockrgt [newrsites] .quantumstat (rsubs);
	  // sign (grrr!)
	  double factor = 1.0;
	  if ((reflection == 2) || (reflection == 3)) 
	    if ((rstat < 0) && (tstat < 0)) factor = -1.0;
	  size_t rrqq = newuniverse .tensororder () [rnew + qq * newrstates];
	  size_t rrtt = newuniverse .tensororder () [rnew + t  * newrstates];
	  for (size_t lnew = 0; lnew < newlstates; lnew++) {
	    size_t llpp = newsystem .tensororder () [lnew + pp * newlstates];
	    newm [llpp + rrqq * newh] = factor * oldm [lnew + rrtt*newlstates]; 
	  }
	}
      }
    }
  }
  else {
    //
    //  Sweep to left 
    //
    Aarray project;
    oldsystem .contraction (ubase, 0, pstates, oldstate, project, false);
    //
    // project [p] [ll pp; r q] = < ll pp | l > * oldstate [l p; r q];
    //
    //  Allocate space for single project matrices and rbase components
    //
    size_t oldmsize = newlstates*pstates * oldrstates*pstates;
    Storage oldfull (oldmsize * sizeof (complex<double>)); 
    complex<double> * oldm = (complex<double> *) oldfull .storage ();
    Storage rbasefull (newrstates * oldrstates * pstates * 
		       sizeof (complex<double>));
    complex<double> * rgtbase = (complex<double> *) rbasefull .storage ();
    rbase .expand (rgtbase);
    //
    //   rgtbase [rr ; r q ] = < rr | r q > || < rr | q r >
    //
    for (size_t p = 0; p < pstates; p++) {
      //
      for (size_t i = 0; i < oldmsize; i++) oldm [i] = 0.0;
      project [p] .expand (oldm);
      //
      //   oldm(p) [ll pp; r q] = < ll pp | l > * oldstate [l p; r q]
      //
      size_t psub  = blocklft [1]  .quantumsub (p);
      long   pstat = blocklft [1]  .quantumstat (psub);
      for (size_t q = 0; q < pstates; q++) {
	//
	//  No reflection or block reflection (reflection == 0 | 2):
	//  
	//  newstate [ll pp; rr qq] = oldm(p) [ll pp; r q] * 
	//			      rgtbase [rr; r p ] < qq | q >
	//                            * sign (reflection)
	//
	size_t qq = q;
	size_t rb = p;
	if ((reflection == 1) || (reflection == 3)) {
	  //
	  //  Universe reflection (reflection == 1 | 3):
	  //
	  //  newstate [ll pp; rr qq] = oldm(p) [ll pp; r q ] 
	  //				rgtbase [rr; r q ] <qq | p>
	  //                            * sign (reflection)
	  qq = p;
	  rb = q;
	}
	size_t qsub   = blockrgt [1] .quantumsub (q);
	size_t qqsub  = blockrgt [1] .quantumsub (qq);
	long   qstat  = blockrgt [1] .quantumstat (qsub);
	long   qqstat = blockrgt [1] .quantumstat (qqsub);
	for (size_t rold = 0; rold < oldrstates; rold++) {
	  size_t rq  = olduniverse .tensororder (rold +  q * oldrstates);
	  size_t rbb = olduniverse .tensororder (rold + rb * oldrstates);
	  size_t roldsub  = blockrgt [oldrsites] .quantumsub  (rold);
	  long   roldstat = blockrgt [oldrsites] .quantumstat (roldsub);
	  for (size_t rnew = 0; rnew < newrstates; rnew++) {
	    size_t rrqq = newuniverse .tensororder (rnew + qq * newrstates);
	    size_t rnewsub = blockrgt [newrsites] .quantumsub  (rnew);
	    long rnewstat  = blockrgt [newrsites] .quantumstat (rnewsub);
	    double factor = 1.0;
	    if ((reflection == 1) || (reflection == 3)) 
	      if ((rnewstat < 0) && (qqstat < 0))  factor *= -1.0;
	    if ((reflection == 2) &&
		(roldstat < 0) && (pstat < 0))  factor *= -1.0;
	    if ((reflection == 3) &&
		(roldstat < 0) && (qstat < 0))  factor *= -1.0;
	    for (size_t llpp = 0; llpp < newlstates * pstates; llpp++) {
	      newm [llpp + rrqq * newlstates * pstates] += 
		(factor * 
		 oldm    [llpp + rq  * newlstates*pstates] * 
		 rgtbase [rnew + rbb *newrstates]); 
	    }  
	  }
	}
      }
    }
  }
  newstate .compress (newm);
}
//
//____________________________________________________________________________
double prond (long k, long l, double t)
{
  double res = 0.5;
  double fac = 2.0 / ((l+1.0) * (l+1.0));
  for (long q = 1; q < l; q++) {
    double kq = M_PI * (q + 1.0) / (l + 1.0); 
    double skq = sin (kq * (k+1));
    double eq  = -2.0 * cos (kq);
    for (long p = 0; p < q; p++) {
      double kp = M_PI * (p + 1.0) / (l + 1.0); 
      double skp = sin (kp * (k+1));
      double ep  = -2.0 * cos (kp);
      double ss = (sin (M_PI * (p-q)/2) / sin ((kp-kq)/2) +
		   sin (M_PI * (p+q)/2) / sin ((kp+kq)/2)) * fac;
      res += skp * skq * ss * cos( (ep-eq) * t );
    }
  }
  return res;
}
//
//____________________________________________________________________________
void properties (Block & system, Block & universe,
                 const Action & vbra, const Action & vket,
                 double norm2, double time)
{
    //
    //	Compute properties or correlations
    //
    vector<Aproperty> actual;
    size_t sites = system .sites () + universe .sites ();
    size_t k, l, n, m, prolen;
    long ta, tb;
    ta = timecpu ();
    prolen = 1;
    //
    //	Remember n. of used actions.
    //
    size_t last_action = name_action ();
    //
    //	Look for wanted correlations with actual chain length
    //
    for (n = 0, m = 0; n < correlation .size (); ++n)
        if (correlation [n] .ap_sites == (sites+1)) {
            actual .push_back (correlation [n]);
            //
            //  Maximum length of properties names (for output aesthetic)
            //
            l = name_define (actual [m] .ap_id) .size ();
            if (l > prolen) prolen = l;
            //
            //	Compute property (use hamaction as super action)
            //
            actual [m] .reorder (sites);
            hamaction = Superaction (actual [m], system, universe, true);
            Action result;
            biapply (result, hamaction, vket, true);
            actual [m] .ap_value = multiply (vbra, result)/norm2;
            m++;
            
            
        }
    //
    //	clear new defined actions (to save memory)
    //

            cout << "correlation before " << report.usage() << " time " << timecpustr (timecpu ()) << endl;
    
//    if (name_action () > last_action) {
        system   .actionclear (last_action, name_action ());
        universe .actionclear (last_action, name_action ());
//    }
    
            cout << "correlation after " << report.usage() << " time " << timecpustr (timecpu ()) << endl;
    
    //
    //  for(auto prop : actual) cout << prop.membro << " " ;
    //
    tb = timecpu () - ta;
    long leading = 28 - prolen;
    if (leading < 0) leading = 0;
    if (m == 0) return;
    cout << "Properties evaluation (cpu " << timecpustr (tb) << "): "
    << setw (28) << right << "complex value" << setw (14) << right
    << "variance" << endl;
    for (n = 0; n < m; ++n) {
        size_t id   = actual [n] .ap_id;
        if (id == 0) continue;
        //
        ofstream pout ((name_define (id) + ".prp").c_str (), ios_base::app);
        if (!pout) {
            cout << "Error opening " << name_define (id) + ".prp" << endl;
            exit (0);
        }
        //pout << "# " << name_define (id) << " L=" << sites << endl;
        //
        complex<double> total  = 0.0;
        double 	    total2 = 0.0;
        size_t	    tcount = 0;
        for (l = n; l < m; ++l) {
            //
            if (actual [l] .ap_id != id) continue;
            //
            long index = actual [l] .ap_index;
            //
            complex<double> partial  = 0.0;
            double	      partial2 = 0.0;
            size_t	      pcount   = 0;
            for (k = l; k < m; ++k) {
                if ((actual [k] .ap_id    != id)  ||
                    (actual [k] .ap_index != index)) continue;
                //
                complex<double> & v = actual [k] .ap_value;
                partial  += v;
                partial2 += (v .real () * v .real () + v .imag () * v .imag ());
                pcount++;
                actual [k] .ap_id = 0;
            }
            partial  /= pcount;
            partial2 /= pcount;
            partial2 -= (partial .real () * partial .real () +
                         partial .imag () * partial .imag ());
            if (partial2 < 0.0) partial2 = 0.0;
            partial2 = sqrt (partial2);
            //
            pout << setw (24) << setprecision (16) << fixed << right
            << time << setw (8) << right << index << setw (8) << " "
            << setw (16) << setprecision (16) << fixed << right
            << showpos << partial .real () << " "
            << setw (16) << setprecision (16) << fixed << right
            << partial .imag () << " "
            << noshowpos << endl;
            /*
             if (leading) cout << setw (leading) << "";
             cout << " < " << setw (prolen) << left << name_define (id)
             << " ( " << setw (4) << right << index << " ) > ";
             cout << setw (29) << setprecision (12) << fixed << right
             << showpos << partial << noshowpos;
             if (pcount > 1)
             cout << setw (20) << setprecision (20) << scientific << right
             << partial2;
             cout << setprecision (20) << endl;
             */
            //
            total  += partial;
            total2 += (partial .real () * partial .real () +
                       partial .imag () * partial .imag ());
            tcount++;
        }
        pout << endl;
        //
        if (tcount > 1) {
            cout << "Sum, Mean "  << setw (29) << setprecision (18)
            << fixed << right << showpos << total << noshowpos << ",";
            total  /= tcount;
            total2 /= tcount;
            total2 -= (total .real () * total .real () +
                       total .imag () * total .imag ());
            if (total2 < 0.0) total2 = 0.0;
            total2 = sqrt (total2);
            cout << setw (29) << setprecision (8) << fixed << right
            << showpos << total << noshowpos;
            cout << setw (9) << setprecision (2) << scientific << right
            << total2;
            cout << setprecision (6) << endl;
        }
    }
    //
    //	reset output options
    //
    cout .setf ((ios_base::fmtflags) 0, ios_base::floatfield);
    cout .precision (6);
}
//

//____________________________________________________________________________
void properties (Block & system, Block & universe)
{
  //
  //	Compute properties or correlations 
  //
  vector<Aproperty> actual;
  size_t sites = system .sites () + universe .sites ();
  size_t k, l, n, m, prolen, bralen, ketlen;
  long ta, tb;
  ta = timecpu ();
  bralen = ketlen = 0;
  prolen = 1;
  //
  //	Remember n. of used actions.
  //
  size_t last_action = name_action ();
  //
  //	Look for wanted correlations with actual chain length
  //
  for (n = 0, m = 0; n < correlation .size (); n++) 
    if ((correlation [n] .ap_sites == sites) ||
	(correlation [n] .ap_sites == 0)) {
      size_t bra, bra_id, bra_index, ket, ket_id, ket_index;
      //
      //  look for |ket> and <bra| 
      //
      bra 	= 0;
      bra_id    = correlation [n] .ap_braid;
      bra_index = correlation [n] .ap_braindex;
      if (bra_id) {
	for (l = 1; l < super_target .subspaces (); l++) 
	  if (bra_id == super_target .id (l)) bra = l;
      }
      else bra = 1;
      if (bra) {
	if (bra_index < super_target .states (bra)) 
	  bra = super_target .offset (bra) + bra_index;
	else continue;
      }
      else continue;
      ket	= 0;
      ket_id    = correlation [n] .ap_ketid;
      ket_index = correlation [n] .ap_ketindex;
      if (ket_id) {
	for (l = 1; l < super_target .subspaces (); l++) 
	  if (ket_id == super_target .id (l)) ket = l;
      }
      else ket = 1;
      if (ket) {
	if (ket_index < super_target .states (ket)) 
	  ket = super_target .offset (ket) + ket_index;
	else continue;
      }
      else continue;      
      //
      actual .push_back (correlation [n]);
      actual [m] .ap_sites    = sites;
      actual [m] .ap_braid    = bra;
      actual [m] .ap_ketid    = ket;
      actual [m] .ap_braindex = n;
      actual [m] .ap_ketindex = n;
      //
      //  Maximum length of properties names (for output aesthetic)
      //
      l = 1;
      if (bra_id) l = name_define (bra_id) .size ();
      if (l > bralen) bralen = l;
      l = 1;
      if (ket_id) l = name_define (ket_id) .size ();
      if (l > ketlen) ketlen = l;
      l = name_define (actual [m] .ap_id) .size ();
      if (l > prolen) prolen = l;
      //
      //	Compute property (use hamaction as super action)
      //
      Action & vket = super_state [ket];
      Action & vbra = super_state [bra];
      actual [m] .reorder (sites);
      hamaction = Superaction (actual [m], system, universe, true);
      Action result;
      biapply (result, hamaction, vket, true);
      actual [m] .ap_value = multiply (vbra, result); 
      m++;
    }
  //
  //	clear new defined actions (to save memory)
  //
  if (name_action () > last_action) 
    system .actionclear (last_action, name_action ());
  //
  tb = timecpu () - ta;
  long leading = 25 - bralen - ketlen - prolen;
  if (leading < 0) leading = 0;
  if (m == 0) return;
  cout << "Properties evaluation (cpu " << timecpustr (tb) << "): "
     << setw (28) << right << "complex value" << setw (14) << right 
     << "variance" << endl;
  for (n = 0; n < m; n++) {
    size_t id        = actual [n] .ap_id;
    if (id == 0) continue;
    //
    size_t bra       = actual [n] .ap_braid;
    size_t ket       = actual [n] .ap_ketid;
    size_t nc 	     = actual [n] .ap_braindex;
    //
    size_t bra_id    = correlation [nc] .ap_braid;
    size_t bra_index = correlation [nc] .ap_braindex;
    size_t ket_id    = correlation [nc] .ap_ketid;
    size_t ket_index = correlation [nc] .ap_ketindex;
    complex<double> total  = 0.0;
    double 	    total2 = 0.0;
    size_t	    tcount = 0;
    for (l = n; l < m; l++) {
      // 
      if ((actual [l] .ap_id    != id) ||
	  (actual [l] .ap_braid != bra) ||
	  (actual [l] .ap_ketid != ket)) continue;
      //
      long index = actual [l] .ap_index;
      //
      complex<double> partial  = 0.0;
      double	      partial2 = 0.0;
      size_t	      pcount   = 0;
      for (k = l; k < m; k++) {
	if ((actual [k] .ap_id    != id)  ||
	    (actual [k] .ap_braid != bra) ||
	    (actual [k] .ap_ketid != ket) ||
	    (actual [k] .ap_index != index)) continue;
	//
	complex<double> & v = actual [k] .ap_value;
	partial  += v;
	partial2 += (v .real () * v .real () + v .imag () * v .imag ());
	pcount++;
	actual [k] .ap_id = 0;
      }
      partial  /= pcount;
      partial2 /= pcount;
      partial2 -= (partial .real () * partial .real () +
		   partial .imag () * partial .imag ());
      if (partial2 < 0.0) partial2 = 0.0;
      partial2 = sqrt (partial2);
      //
      if (leading) cout << setw (leading) << "";
      cout << "<";
      if (bra_id) cout << setw (bralen) << left << name_define (bra_id);
      else cout << setw (bralen) << left << "0";
      cout << setw (2) << right << bra_index;
      cout << "| " << setw (prolen) << left << name_define (id)
	   << " (" << setw (4) << right << index << ")|";
      if (ket_id) cout << setw (ketlen) << left << name_define (ket_id);
      else cout << setw (ketlen) << left << "0";
      cout << setw (2) << right << ket_index << "> = ";
      cout << setw (24) << setprecision (8) << fixed << right 
	   << showpos << partial << noshowpos;
      if (pcount > 1) 
	cout << setw (9) << setprecision (2) << scientific << right 
	     << partial2;
      cout << setprecision (6) << endl;
      //
      total  += partial;
      total2 += (partial .real () * partial .real () +
		 partial .imag () * partial .imag ());
      tcount++;
    }
    //			 
    if (tcount > 1) {
      cout << "Sum, Mean "  << setw (29) << setprecision (8) 
	   << fixed << right << showpos << total << noshowpos << ",";
      total  /= tcount;
      total2 /= tcount;
      total2 -= (total .real () * total .real () + 
		 total .imag () * total .imag ());
      if (total2 < 0.0) total2 = 0.0;
      total2 = sqrt (total2);
      cout << setw (29) << setprecision (8) << fixed << right 
	   << showpos << total << noshowpos;
      cout << setw (9) << setprecision (2) << scientific << right 
	   << total2;
      cout << setprecision (6) << endl;
    }
  }
  //
  //	reset output options
  //
  cout .setf ((ios_base::fmtflags) 0, ios_base::floatfield);
  cout .precision (6);
}
//
//____________________________________________________________________________
void	updatebase (Block & system, Block & universe, long tag)
{
  //	
  //	update basis states of system or universe depending on tag.
  //
  //	Parse eventual operators to be applyed to target states
  //
  densityop_parse (system .sites (), universe .sites ());
  Action densitylft (system   .quantum (), 0.0);
  Action densityrgt (universe .quantum (), 0.0);
  Action actual;
  double factor;
  for (size_t n = 0; n < super_state .size (); n++) {
    Action ket = super_state [n];
    Action bra = ket; 
    bra .dagger ();
    if (tag >= 0) {	// system reduced density 
      actual  = ket;
      actual *= bra;
      actual *= super_weight [n];
      densitylft += actual;
    }
    if (tag <= 0) {	// universe reduced density
      actual  = bra;
      actual *= ket;
      actual *= super_weight [n];
      densityrgt += actual;
    }
    //
    //	Apply tensor components of hamiltonian to the target states
    //
    if (tensorweight > 0.0) {
      double weight = 0.0;
      for (size_t j = 0; j < hamaction .size (); j++) {
	Biaction & ba = hamaction [j];
	if (! ba .ba_rgtop .isscalar ()) weight += 1.0;
	if (! ba .ba_lftop .isscalar ()) weight += 1.0;
      }
      //cout << "weight " << weight << endl;
      if (weight > 0.0) weight = tensorweight /weight;
      for (size_t j = 0; j < hamaction .size (); j++) {
	Biaction & ba = hamaction [j];
	if (! ba .ba_rgtop .isscalar ()) {
	  multiply (ket, super_state [n], ba .ba_rgtop, true);
	  factor = multiply (ket, ket) .real ();
	  if (factor  > 1.e-15) { 
	    bra = ket; 
	    bra .dagger ();
	    if (tag >= 0) {
	      actual  = ket;
	      actual *= bra;
	      actual *= (weight/factor);
	      densitylft += actual;
	    }
	    if (tag <= 0) {
	      actual  = bra;
	      actual *= ket;
	      actual *= (weight/factor);
	      densityrgt += actual;
	    }
	  }
	}
	if (! ba .ba_lftop .isscalar ()) {
	  multiply (ket, ba .ba_lftop, super_state [n]);
	  factor = multiply (ket,ket) .real ();
	  if (factor > 1.e-15) {
	    bra = ket; 
	    bra .dagger ();
	    if (tag >= 0) {
	      actual  = ket;
	      actual *= bra;
	      actual *= (weight/factor);
	      densitylft += actual;
	    }
	    if (tag <= 0) {
	      actual  = bra;
	      actual *= ket;
	      actual *= (weight/factor);
	      densityrgt += actual;
	    }
	  }
	}
      }
    }
    //
    //	Apply wanted operators and compute corresponding density matrices
    //
    for (size_t nop = 0; nop < densityoperator .size (); nop++) {
      initialaction = Superaction (densityoperator [nop], 
				   system, universe, true);
      //cout << "stateblocks " << super_state [n] .blocks () 
      //	   << " " << super_state [n] .size () << endl;
      //super_state [n] .show ("stato");
      biapply (ket, initialaction, super_state [n], true);
      factor = multiply (ket, ket) .real ();
      if (factor <= 1.e-15) continue;
      //
      //	temporary cleanup (to do)
      //
      size_t dim = ket .width () * ket .height ();
      Storage st (dim * sizeof (complex<double>));
      complex<double> * m = (complex<double> *) st .storage ();
      ket .expand (m);
      ket .compress (m);
      //cout << "nblocks " << ket .blocks () <<" " << ket .size () <<  endl;
      //ket .show ("cdag");
      bra = ket;
      bra .dagger ();
      if (tag >= 0) {	// system reduced density 
	actual  = ket;
	actual *= bra;
	densitylft += actual;
	//actual .show ("dopsys");
      }
      if (tag <= 0) {	// universe reduced density
	actual  = bra;
	actual *= ket;
	densityrgt += actual;
	//actual .show ("dopuni");
      }
    }
  }
  //
  // 	update basis
  //
  Action sum;
  if (tag >= 0) {
    if (symmetry .size ()) {
      sum = densitylft;
      for (size_t ns = 0; ns < symmetry .size (); ns++) {
	Action sym    = system .symmetry (symmetry [ns]);
	Action actual = densitylft;
	actual *= sym;
	sym .dagger ();
	sym *= actual;
	sum += sym;
      }
      densitylft = sum;
    }
    densitylft .clean ();
//	cout << "LFT " << min_lft <<" " << max_lft << " "<< n_cutlft << endl;
    system .select_states (densitylft, min_lft, max_lft, n_cutlft);
  }
  if (tag <= 0) { 
    densityrgt .transpose ();
    if (symmetry .size ()) {
      sum = densityrgt;
      for (size_t ns = 0; ns < symmetry .size (); ns++) {
	Action sym    = universe .symmetry (symmetry [ns]);
	Action actual = densityrgt;
	actual *= sym;
	sym .dagger ();
	sym *= actual;
	sum += sym;
      }
      densityrgt = sum;
    }
    densityrgt .clean ();

//	cout << "RGT " << min_rgt <<" " << max_rgt << " "<< n_cutrgt << endl;
    universe .select_states (densityrgt, min_rgt, max_rgt, n_cutrgt);
  }
} 
//
//============================================================================

