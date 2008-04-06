#include <../../nrnconf.h>
#include "nrnredef.h"
#include "section.h"
#include "nrnmpiuse.h"
/* stubs for nrnoc. The actual functions are for interviews menus */
nrnallsectionmenu(){ret(0);}
nrnallpointmenu(){ret(0);}
nrnsecmenu(){ret(0);}
nrnglobalmechmenu(){ret(0);}
nrnmechmenu(){ret(0);}
nrnpointmenu(){ret(0);}
make_mechanism(){ret(0);}
make_pointprocess(){ret(0);}
nrnpython() {ret(0);}
/*ARGSUSED*/
void hoc_construct_point(ob) void* ob; {}
nrn_random_play(){}
void* nrn_random_arg(i) int i; { return (void*)0; }
double nrn_random_pick(r) void* r; { return 0.; }
hoc_new_opoint(){}
special_pnt_call(){}
bbs_handle(){}
linmod_alloc(){}
int linmod_extra_eqn_count() { return 0; }
linmod_rhs(){}
linmod_lhs(){}
linmod_update(){}
linmod_init(){}
nrn_solver_prepare(){}
nrn_fihexec(i) int i;{}
nrn_deliver_events(tt) double tt;{}
nrn_record_init(){}
nrn_play_init(){}
fixed_record_continuous(){}
fixed_play_continuous(){}
nrniv_recalc_ptrs(){}

int at_time(te) double te; {
	extern double t, dt;
	double x = te - 1e-11;
	if (x > (t - dt)  && x <= t) {
		return 1;
	}else{
		return 0;
	}
}

net_event(){hoc_execerror("net_event only available in nrniv", (char*)0);}
net_send(){hoc_execerror("net_send only available in nrniv", (char*)0);}
artcell_net_send(){hoc_execerror("net_send only available in nrniv", (char*)0);}
net_move(){hoc_execerror("net_move only available in nrniv", (char*)0);}
artcell_net_move(){hoc_execerror("net_move only available in nrniv", (char*)0);}
nrn_use_daspk(i) int i; { }

#if CVODE
cvode_fadvance(t)double t;{}
cvode_finitialize(){}
cvode_event(x) double x;{}
clear_event_queue(){}
init_net_events(){}
deliver_net_events(){}
hoc_reg_singlechan(){}
_singlechan_declare(){}
_nrn_single_react(){}
#endif

#if defined(CYGWIN)
void* dll_lookup(s) char* s; {return 0;}
void* dll_load(v, s) void* v; char* s; {return 0;}
#endif

void nrn_spike_exchange_init(){}
void nrn_spike_exchange(){}
void nrn_fake_fire(i,x,j)int i, j; double x;{}
void nrn_daq_ao() {}
void nrn_daq_ai() {}
void nrn_daq_scanstart(){}
void nrnpy_reg_mech(int i){}

#if PARANEURON
double* nrn_classicalNodeA(Node* n) {return (double*)0;}
double* nrn_classicalNodeB(Node* n) {return (double*)0;}
#endif
