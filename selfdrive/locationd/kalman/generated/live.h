/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_912188989712334113);
void inv_err_fun(double *nom_x, double *true_x, double *out_6896580442647649017);
void H_mod_fun(double *state, double *out_8002970078773088862);
void f_fun(double *state, double dt, double *out_2545119348829883043);
void F_fun(double *state, double dt, double *out_1665946376783243158);
void h_3(double *state, double *unused, double *out_7121274795261509905);
void H_3(double *state, double *unused, double *out_4835243617977677583);
void h_4(double *state, double *unused, double *out_5575466453095867122);
void H_4(double *state, double *unused, double *out_1884030708145280304);
void h_9(double *state, double *unused, double *out_6837960719668320419);
void H_9(double *state, double *unused, double *out_718902034028390712);
void h_10(double *state, double *unused, double *out_7470992672723713596);
void H_10(double *state, double *unused, double *out_9029471514832757054);
void h_12(double *state, double *unused, double *out_3936529178792141977);
void H_12(double *state, double *unused, double *out_3292588560476036664);
void h_13(double *state, double *unused, double *out_2997693686268920610);
void H_13(double *state, double *unused, double *out_3408344647633756642);
void h_14(double *state, double *unused, double *out_6837960719668320419);
void H_14(double *state, double *unused, double *out_718902034028390712);
void h_19(double *state, double *unused, double *out_6494841857451912197);
void H_19(double *state, double *unused, double *out_8813298631142486660);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);