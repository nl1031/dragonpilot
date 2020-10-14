/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8204215386578902715);
void inv_err_fun(double *nom_x, double *true_x, double *out_2145089199424687675);
void H_mod_fun(double *state, double *out_8064926254512108318);
void f_fun(double *state, double dt, double *out_5875409788587754688);
void F_fun(double *state, double dt, double *out_3415275123736377020);
void h_25(double *state, double *unused, double *out_3239016220435494339);
void H_25(double *state, double *unused, double *out_120949934209931282);
void h_24(double *state, double *unused, double *out_8895723921966059180);
void H_24(double *state, double *unused, double *out_6180294348718162794);
void h_26(double *state, double *unused, double *out_6371035696055707836);
void H_26(double *state, double *unused, double *out_3978496019381186698);
void h_27(double *state, double *unused, double *out_965225713592976720);
void H_27(double *state, double *unused, double *out_2092283034650059750);
void h_29(double *state, double *unused, double *out_8234142898844024910);
void H_29(double *state, double *unused, double *out_2561626517312303074);
void h_28(double *state, double *unused, double *out_8127271147273803940);
void H_28(double *state, double *unused, double *out_9180081764526761856);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
