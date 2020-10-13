/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5936600902336233121);
void inv_err_fun(double *nom_x, double *true_x, double *out_907009349114483511);
void H_mod_fun(double *state, double *out_1617323043114287287);
void f_fun(double *state, double dt, double *out_2358448951961011578);
void F_fun(double *state, double dt, double *out_3534033638387069784);
void h_3(double *state, double *unused, double *out_3609830578751498065);
void H_3(double *state, double *unused, double *out_8349724511825558515);
void h_4(double *state, double *unused, double *out_2901869864305581470);
void H_4(double *state, double *unused, double *out_743938497470033113);
void h_9(double *state, double *unused, double *out_8698656443091399429);
void H_9(double *state, double *unused, double *out_4502690894696551477);
void h_10(double *state, double *unused, double *out_5592994719718198204);
void H_10(double *state, double *unused, double *out_4928365509106764252);
void h_12(double *state, double *unused, double *out_708305998833796345);
void H_12(double *state, double *unused, double *out_4448431500690824261);
void h_13(double *state, double *unused, double *out_6110666387541124501);
void H_13(double *state, double *unused, double *out_2785232682424019824);
void h_14(double *state, double *unused, double *out_8698656443091399429);
void H_14(double *state, double *unused, double *out_4502690894696551477);
void h_19(double *state, double *unused, double *out_1837752462832572563);
void H_19(double *state, double *unused, double *out_9135609480248333551);
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