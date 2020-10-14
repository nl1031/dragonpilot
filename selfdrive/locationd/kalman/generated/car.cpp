
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8204215386578902715) {
   out_8204215386578902715[0] = delta_x[0] + nom_x[0];
   out_8204215386578902715[1] = delta_x[1] + nom_x[1];
   out_8204215386578902715[2] = delta_x[2] + nom_x[2];
   out_8204215386578902715[3] = delta_x[3] + nom_x[3];
   out_8204215386578902715[4] = delta_x[4] + nom_x[4];
   out_8204215386578902715[5] = delta_x[5] + nom_x[5];
   out_8204215386578902715[6] = delta_x[6] + nom_x[6];
   out_8204215386578902715[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2145089199424687675) {
   out_2145089199424687675[0] = -nom_x[0] + true_x[0];
   out_2145089199424687675[1] = -nom_x[1] + true_x[1];
   out_2145089199424687675[2] = -nom_x[2] + true_x[2];
   out_2145089199424687675[3] = -nom_x[3] + true_x[3];
   out_2145089199424687675[4] = -nom_x[4] + true_x[4];
   out_2145089199424687675[5] = -nom_x[5] + true_x[5];
   out_2145089199424687675[6] = -nom_x[6] + true_x[6];
   out_2145089199424687675[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_8064926254512108318) {
   out_8064926254512108318[0] = 1.0;
   out_8064926254512108318[1] = 0.0;
   out_8064926254512108318[2] = 0.0;
   out_8064926254512108318[3] = 0.0;
   out_8064926254512108318[4] = 0.0;
   out_8064926254512108318[5] = 0.0;
   out_8064926254512108318[6] = 0.0;
   out_8064926254512108318[7] = 0.0;
   out_8064926254512108318[8] = 0.0;
   out_8064926254512108318[9] = 1.0;
   out_8064926254512108318[10] = 0.0;
   out_8064926254512108318[11] = 0.0;
   out_8064926254512108318[12] = 0.0;
   out_8064926254512108318[13] = 0.0;
   out_8064926254512108318[14] = 0.0;
   out_8064926254512108318[15] = 0.0;
   out_8064926254512108318[16] = 0.0;
   out_8064926254512108318[17] = 0.0;
   out_8064926254512108318[18] = 1.0;
   out_8064926254512108318[19] = 0.0;
   out_8064926254512108318[20] = 0.0;
   out_8064926254512108318[21] = 0.0;
   out_8064926254512108318[22] = 0.0;
   out_8064926254512108318[23] = 0.0;
   out_8064926254512108318[24] = 0.0;
   out_8064926254512108318[25] = 0.0;
   out_8064926254512108318[26] = 0.0;
   out_8064926254512108318[27] = 1.0;
   out_8064926254512108318[28] = 0.0;
   out_8064926254512108318[29] = 0.0;
   out_8064926254512108318[30] = 0.0;
   out_8064926254512108318[31] = 0.0;
   out_8064926254512108318[32] = 0.0;
   out_8064926254512108318[33] = 0.0;
   out_8064926254512108318[34] = 0.0;
   out_8064926254512108318[35] = 0.0;
   out_8064926254512108318[36] = 1.0;
   out_8064926254512108318[37] = 0.0;
   out_8064926254512108318[38] = 0.0;
   out_8064926254512108318[39] = 0.0;
   out_8064926254512108318[40] = 0.0;
   out_8064926254512108318[41] = 0.0;
   out_8064926254512108318[42] = 0.0;
   out_8064926254512108318[43] = 0.0;
   out_8064926254512108318[44] = 0.0;
   out_8064926254512108318[45] = 1.0;
   out_8064926254512108318[46] = 0.0;
   out_8064926254512108318[47] = 0.0;
   out_8064926254512108318[48] = 0.0;
   out_8064926254512108318[49] = 0.0;
   out_8064926254512108318[50] = 0.0;
   out_8064926254512108318[51] = 0.0;
   out_8064926254512108318[52] = 0.0;
   out_8064926254512108318[53] = 0.0;
   out_8064926254512108318[54] = 1.0;
   out_8064926254512108318[55] = 0.0;
   out_8064926254512108318[56] = 0.0;
   out_8064926254512108318[57] = 0.0;
   out_8064926254512108318[58] = 0.0;
   out_8064926254512108318[59] = 0.0;
   out_8064926254512108318[60] = 0.0;
   out_8064926254512108318[61] = 0.0;
   out_8064926254512108318[62] = 0.0;
   out_8064926254512108318[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5875409788587754688) {
   out_5875409788587754688[0] = state[0];
   out_5875409788587754688[1] = state[1];
   out_5875409788587754688[2] = state[2];
   out_5875409788587754688[3] = state[3];
   out_5875409788587754688[4] = state[4];
   out_5875409788587754688[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5875409788587754688[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5875409788587754688[7] = state[7];
}
void F_fun(double *state, double dt, double *out_3415275123736377020) {
   out_3415275123736377020[0] = 1;
   out_3415275123736377020[1] = 0;
   out_3415275123736377020[2] = 0;
   out_3415275123736377020[3] = 0;
   out_3415275123736377020[4] = 0;
   out_3415275123736377020[5] = 0;
   out_3415275123736377020[6] = 0;
   out_3415275123736377020[7] = 0;
   out_3415275123736377020[8] = 0;
   out_3415275123736377020[9] = 1;
   out_3415275123736377020[10] = 0;
   out_3415275123736377020[11] = 0;
   out_3415275123736377020[12] = 0;
   out_3415275123736377020[13] = 0;
   out_3415275123736377020[14] = 0;
   out_3415275123736377020[15] = 0;
   out_3415275123736377020[16] = 0;
   out_3415275123736377020[17] = 0;
   out_3415275123736377020[18] = 1;
   out_3415275123736377020[19] = 0;
   out_3415275123736377020[20] = 0;
   out_3415275123736377020[21] = 0;
   out_3415275123736377020[22] = 0;
   out_3415275123736377020[23] = 0;
   out_3415275123736377020[24] = 0;
   out_3415275123736377020[25] = 0;
   out_3415275123736377020[26] = 0;
   out_3415275123736377020[27] = 1;
   out_3415275123736377020[28] = 0;
   out_3415275123736377020[29] = 0;
   out_3415275123736377020[30] = 0;
   out_3415275123736377020[31] = 0;
   out_3415275123736377020[32] = 0;
   out_3415275123736377020[33] = 0;
   out_3415275123736377020[34] = 0;
   out_3415275123736377020[35] = 0;
   out_3415275123736377020[36] = 1;
   out_3415275123736377020[37] = 0;
   out_3415275123736377020[38] = 0;
   out_3415275123736377020[39] = 0;
   out_3415275123736377020[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3415275123736377020[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3415275123736377020[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3415275123736377020[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3415275123736377020[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3415275123736377020[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3415275123736377020[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3415275123736377020[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3415275123736377020[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3415275123736377020[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3415275123736377020[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3415275123736377020[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3415275123736377020[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3415275123736377020[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3415275123736377020[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3415275123736377020[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3415275123736377020[56] = 0;
   out_3415275123736377020[57] = 0;
   out_3415275123736377020[58] = 0;
   out_3415275123736377020[59] = 0;
   out_3415275123736377020[60] = 0;
   out_3415275123736377020[61] = 0;
   out_3415275123736377020[62] = 0;
   out_3415275123736377020[63] = 1;
}
void h_25(double *state, double *unused, double *out_3239016220435494339) {
   out_3239016220435494339[0] = state[6];
}
void H_25(double *state, double *unused, double *out_120949934209931282) {
   out_120949934209931282[0] = 0;
   out_120949934209931282[1] = 0;
   out_120949934209931282[2] = 0;
   out_120949934209931282[3] = 0;
   out_120949934209931282[4] = 0;
   out_120949934209931282[5] = 0;
   out_120949934209931282[6] = 1;
   out_120949934209931282[7] = 0;
}
void h_24(double *state, double *unused, double *out_8895723921966059180) {
   out_8895723921966059180[0] = state[4];
   out_8895723921966059180[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6180294348718162794) {
   out_6180294348718162794[0] = 0;
   out_6180294348718162794[1] = 0;
   out_6180294348718162794[2] = 0;
   out_6180294348718162794[3] = 0;
   out_6180294348718162794[4] = 1;
   out_6180294348718162794[5] = 0;
   out_6180294348718162794[6] = 0;
   out_6180294348718162794[7] = 0;
   out_6180294348718162794[8] = 0;
   out_6180294348718162794[9] = 0;
   out_6180294348718162794[10] = 0;
   out_6180294348718162794[11] = 0;
   out_6180294348718162794[12] = 0;
   out_6180294348718162794[13] = 1;
   out_6180294348718162794[14] = 0;
   out_6180294348718162794[15] = 0;
}
void h_26(double *state, double *unused, double *out_6371035696055707836) {
   out_6371035696055707836[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3978496019381186698) {
   out_3978496019381186698[0] = 0;
   out_3978496019381186698[1] = 0;
   out_3978496019381186698[2] = 0;
   out_3978496019381186698[3] = 0;
   out_3978496019381186698[4] = 0;
   out_3978496019381186698[5] = 0;
   out_3978496019381186698[6] = 0;
   out_3978496019381186698[7] = 1;
}
void h_27(double *state, double *unused, double *out_965225713592976720) {
   out_965225713592976720[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2092283034650059750) {
   out_2092283034650059750[0] = 0;
   out_2092283034650059750[1] = 0;
   out_2092283034650059750[2] = 0;
   out_2092283034650059750[3] = 1;
   out_2092283034650059750[4] = 0;
   out_2092283034650059750[5] = 0;
   out_2092283034650059750[6] = 0;
   out_2092283034650059750[7] = 0;
}
void h_29(double *state, double *unused, double *out_8234142898844024910) {
   out_8234142898844024910[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2561626517312303074) {
   out_2561626517312303074[0] = 0;
   out_2561626517312303074[1] = 1;
   out_2561626517312303074[2] = 0;
   out_2561626517312303074[3] = 0;
   out_2561626517312303074[4] = 0;
   out_2561626517312303074[5] = 0;
   out_2561626517312303074[6] = 0;
   out_2561626517312303074[7] = 0;
}
void h_28(double *state, double *unused, double *out_8127271147273803940) {
   out_8127271147273803940[0] = state[5];
   out_8127271147273803940[1] = state[6];
}
void H_28(double *state, double *unused, double *out_9180081764526761856) {
   out_9180081764526761856[0] = 0;
   out_9180081764526761856[1] = 0;
   out_9180081764526761856[2] = 0;
   out_9180081764526761856[3] = 0;
   out_9180081764526761856[4] = 0;
   out_9180081764526761856[5] = 1;
   out_9180081764526761856[6] = 0;
   out_9180081764526761856[7] = 0;
   out_9180081764526761856[8] = 0;
   out_9180081764526761856[9] = 0;
   out_9180081764526761856[10] = 0;
   out_9180081764526761856[11] = 0;
   out_9180081764526761856[12] = 0;
   out_9180081764526761856[13] = 0;
   out_9180081764526761856[14] = 1;
   out_9180081764526761856[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
