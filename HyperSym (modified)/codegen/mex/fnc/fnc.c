/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fnc.c
 *
 * Code generation for function 'fnc'
 *
 */

/* Include files */
#include "fnc.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo j_emlrtRSI = {
    63,    /* lineNo */
    "fnc", /* fcnName */
    "C:\\Users\\vinic\\Dropbox\\PUC\\Artigos\\Modelos Hiperel\xc3\xa1sticos e "
    "Matem\xc3\xa1tica Simb\xc3\xb3lica no MATLAB\\C\xc3\xb3"
    "digos\\NLFEA (modified) + HyperSym\\HyperSym [Modified]\\fnc.m" /* pathName
                                                                      */
};

static emlrtRTEInfo emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m" /* pName
                                                                       */
};

/* Function Definitions */
void fnc(const emlrtStack *sp, const real_T in1[9], const real_T in2[6],
         real_T Svol[6], real_T Dvol[36])
{
  emlrtStack st;
  real_T Svol_tmp;
  real_T b_Svol_tmp;
  real_T c_Svol_tmp;
  real_T d_Svol_tmp;
  real_T e_Svol_tmp;
  real_T f_Svol_tmp;
  real_T g_Svol_tmp;
  real_T h_Svol_tmp;
  real_T i_Svol_tmp;
  real_T j_Svol_tmp;
  real_T k_Svol_tmp;
  real_T t14;
  real_T t15;
  real_T t16;
  real_T t17;
  real_T t18;
  real_T t19;
  real_T t2;
  real_T t20;
  real_T t21;
  real_T t218;
  real_T t219;
  real_T t220;
  real_T t221;
  real_T t222;
  real_T t223;
  real_T t224;
  real_T t225;
  real_T t226;
  real_T t227;
  real_T t229;
  real_T t3;
  real_T t33;
  real_T t34;
  real_T t35;
  real_T t36;
  real_T t37;
  real_T t38;
  real_T t39;
  real_T t4;
  real_T t40;
  real_T t41;
  real_T t42;
  real_T t43;
  real_T t44;
  real_T t47;
  real_T t48;
  real_T t49;
  real_T t5;
  real_T t50;
  real_T t52;
  real_T t53;
  real_T t54;
  real_T t55;
  real_T t6;
  st.prev = sp;
  st.tls = sp->tls;
  /* Yeoh_vol */
  /*     [Svol,Dvol] = Yeoh_vol(IN1,IN2) */
  /*     This function was generated by the Symbolic Math Toolbox version 9.3.
   */
  /*     19-Dec-2024 11:57:52 */
  t2 = in1[3] * in1[3];
  t3 = in1[6] * in1[6];
  t4 = in1[7] * in1[7];
  t5 = in1[3] * in1[6];
  t6 = in1[0] * in1[4];
  t15 = 1.0 / in2[3];
  t16 = 1.0 / in2[4];
  t17 = 1.0 / in2[5];
  t14 = in1[8] * t6;
  t18 = in1[7] * t5 * 2.0;
  t19 = in1[0] * t4;
  t20 = in1[4] * t3;
  t21 = in1[8] * t2;
  t33 = t5 - in1[0] * in1[7];
  t34 = in1[3] * in1[7] - in1[4] * in1[6];
  t35 = in1[6] * in1[7] - in1[3] * in1[8];
  t36 = t2 - t6;
  t37 = t3 - in1[0] * in1[8];
  t38 = t4 - in1[4] * in1[8];
  t39 = t36 * t36;
  t40 = t37 * t37;
  t41 = t38 * t38;
  t42 = t33 * t33;
  t43 = t34 * t34;
  t44 = t35 * t35;
  t47 = 1.0 / ((((t19 + t20) + t21) - t18) - t14);
  st.site = &j_emlrtRSI;
  t5 = (((t14 + t18) - t19) - t20) - t21;
  if (t5 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  t48 = muDoubleScalarSqrt(t5);
  t49 = 1.0 / t48;
  t50 = muDoubleScalarPower(t49, 3.0);
  t52 = (t48 - 1.0) * (t48 - 1.0);
  t53 = muDoubleScalarPower(t48 - 1.0, 3.0);
  t55 = muDoubleScalarPower(t48 - 1.0, 5.0);
  Svol[0] =
      (t15 * t38 * t49 * (t48 - 1.0) * -2.0 - t16 * t38 * t49 * t53 * 4.0) -
      t17 * t38 * t49 * t55 * 6.0;
  Svol_tmp = t15 * t37;
  b_Svol_tmp = t16 * t37;
  c_Svol_tmp = t17 * t37;
  Svol[1] =
      (Svol_tmp * t49 * (t48 - 1.0) * -2.0 - b_Svol_tmp * t49 * t53 * 4.0) -
      c_Svol_tmp * t49 * t55 * 6.0;
  d_Svol_tmp = t15 * t36;
  e_Svol_tmp = t16 * t36;
  f_Svol_tmp = t17 * t36;
  Svol[2] =
      (d_Svol_tmp * t49 * (t48 - 1.0) * -2.0 - e_Svol_tmp * t49 * t53 * 4.0) -
      f_Svol_tmp * t49 * t55 * 6.0;
  t3 = t15 * t35;
  t4 = t16 * t35;
  t14 = t17 * t35;
  Svol[3] = (t3 * t49 * (t48 - 1.0) * 2.0 + t4 * t49 * t53 * 4.0) +
            t14 * t49 * t55 * 6.0;
  g_Svol_tmp = t15 * t33;
  h_Svol_tmp = t16 * t33;
  i_Svol_tmp = t17 * t33;
  Svol[4] =
      (g_Svol_tmp * t49 * (t48 - 1.0) * 2.0 + h_Svol_tmp * t49 * t53 * 4.0) +
      i_Svol_tmp * t49 * t55 * 6.0;
  t229 = t15 * t34;
  j_Svol_tmp = t16 * t34;
  k_Svol_tmp = t17 * t34;
  Svol[5] = (t229 * t49 * (t48 - 1.0) * 2.0 + j_Svol_tmp * t49 * t53 * 4.0) +
            k_Svol_tmp * t49 * t55 * 6.0;
  t54 = t52 * t52;
  t5 = g_Svol_tmp * t36;
  t2 = h_Svol_tmp * t36;
  t6 = i_Svol_tmp * t36;
  t218 = ((((t5 * t47 * 2.0 + t2 * t47 * t52 * 12.0) + t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t5 = t229 * t36;
  t2 = j_Svol_tmp * t36;
  t6 = k_Svol_tmp * t36;
  t219 = ((((t5 * t47 * 2.0 + t2 * t47 * t52 * 12.0) + t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t5 = g_Svol_tmp * t37;
  t2 = h_Svol_tmp * t37;
  t6 = i_Svol_tmp * t37;
  t220 = ((((t5 * t47 * 2.0 + t2 * t47 * t52 * 12.0) + t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t5 = t3 * t37;
  t2 = t4 * t37;
  t6 = t14 * t37;
  t221 = ((((t5 * t47 * 2.0 + t2 * t47 * t52 * 12.0) + t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t5 = t229 * t38;
  t2 = j_Svol_tmp * t38;
  t6 = k_Svol_tmp * t38;
  t222 = ((((t5 * t47 * 2.0 + t2 * t47 * t52 * 12.0) + t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t5 = t3 * t38;
  t2 = t4 * t38;
  t6 = t14 * t38;
  t223 = ((((t5 * t47 * 2.0 + t2 * t47 * t52 * 12.0) + t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t6 = t3 * t36;
  t2 = t4 * t36;
  t5 = t14 * t36;
  t3 = in1[3] * t15 * t49 * (t48 - 1.0);
  t4 = in1[3] * t16 * t49 * t53;
  t14 = in1[3] * t17 * t49 * t55;
  t224 = (((((((t6 * t47 * 2.0 - t3 * 4.0) - t4 * 8.0) - t14 * 12.0) +
             t2 * t47 * t52 * 12.0) +
            t5 * t47 * t54 * 30.0) +
           t6 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t5 * t50 * t55 * 6.0;
  t5 = t229 * t37;
  t2 = j_Svol_tmp * t37;
  t6 = k_Svol_tmp * t37;
  t18 = in1[6] * t15 * t49 * (t48 - 1.0);
  t19 = in1[6] * t16 * t49 * t53;
  t20 = in1[6] * t17 * t49 * t55;
  t225 = (((((((t5 * t47 * 2.0 - t18 * 4.0) - t19 * 8.0) - t20 * 12.0) +
             t2 * t47 * t52 * 12.0) +
            t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t5 = g_Svol_tmp * t38;
  t2 = h_Svol_tmp * t38;
  t6 = i_Svol_tmp * t38;
  t21 = in1[7] * t15 * t49 * (t48 - 1.0);
  t33 = in1[7] * t16 * t49 * t53;
  t36 = in1[7] * t17 * t49 * t55;
  t226 = (((((((t5 * t47 * 2.0 - t21 * 4.0) - t33 * 8.0) - t36 * 12.0) +
             t2 * t47 * t52 * 12.0) +
            t6 * t47 * t54 * 30.0) +
           t5 * t50 * (t48 - 1.0) * 2.0) +
          t2 * t50 * t53 * 4.0) +
         t6 * t50 * t55 * 6.0;
  t2 = g_Svol_tmp * t34;
  t6 = h_Svol_tmp * t34;
  t5 = i_Svol_tmp * t34;
  t227 = (((((((-(t2 * t47 * 2.0) + t3 * 2.0) + t4 * 4.0) + t14 * 6.0) -
             t6 * t47 * t52 * 12.0) -
            t5 * t47 * t54 * 30.0) -
           t2 * t50 * (t48 - 1.0) * 2.0) -
          t6 * t50 * t53 * 4.0) -
         t5 * t50 * t55 * 6.0;
  t6 = g_Svol_tmp * t35;
  t2 = h_Svol_tmp * t35;
  t5 = i_Svol_tmp * t35;
  t34 = (((((((-(t6 * t47 * 2.0) + t18 * 2.0) + t19 * 4.0) + t20 * 6.0) -
            t2 * t47 * t52 * 12.0) -
           t5 * t47 * t54 * 30.0) -
          t6 * t50 * (t48 - 1.0) * 2.0) -
         t2 * t50 * t53 * 4.0) -
        t5 * t50 * t55 * 6.0;
  t6 = t229 * t35;
  t2 = j_Svol_tmp * t35;
  t5 = k_Svol_tmp * t35;
  t229 = (((((((-(t6 * t47 * 2.0) + t21 * 2.0) + t33 * 4.0) + t36 * 6.0) -
             t2 * t47 * t52 * 12.0) -
            t5 * t47 * t54 * 30.0) -
           t6 * t50 * (t48 - 1.0) * 2.0) -
          t2 * t50 * t53 * 4.0) -
         t5 * t50 * t55 * 6.0;
  t2 = d_Svol_tmp * t37;
  t6 = e_Svol_tmp * t37;
  t5 = f_Svol_tmp * t37;
  g_Svol_tmp = in1[0] * t15 * t49 * (t48 - 1.0);
  h_Svol_tmp = in1[0] * t16 * t49 * t53;
  i_Svol_tmp = in1[0] * t17 * t49 * t55;
  t36 = (((((((-(t2 * t47 * 2.0) + g_Svol_tmp * 4.0) + h_Svol_tmp * 8.0) +
             i_Svol_tmp * 12.0) -
            t6 * t47 * t52 * 12.0) -
           t5 * t47 * t54 * 30.0) -
          t2 * t50 * (t48 - 1.0) * 2.0) -
         t6 * t50 * t53 * 4.0) -
        t5 * t50 * t55 * 6.0;
  t6 = d_Svol_tmp * t38;
  t2 = e_Svol_tmp * t38;
  t5 = f_Svol_tmp * t38;
  t20 = in1[4] * t15 * t49 * (t48 - 1.0);
  t21 = in1[4] * t16 * t49 * t53;
  t33 = in1[4] * t17 * t49 * t55;
  t19 = (((((((-(t6 * t47 * 2.0) + t20 * 4.0) + t21 * 8.0) + t33 * 12.0) -
            t2 * t47 * t52 * 12.0) -
           t5 * t47 * t54 * 30.0) -
          t6 * t50 * (t48 - 1.0) * 2.0) -
         t2 * t50 * t53 * 4.0) -
        t5 * t50 * t55 * 6.0;
  t6 = Svol_tmp * t38;
  t2 = b_Svol_tmp * t38;
  t5 = c_Svol_tmp * t38;
  t14 = in1[8] * t15 * t49 * (t48 - 1.0);
  t18 = in1[8] * t16 * t49 * t53;
  t4 = in1[8] * t17 * t49 * t55;
  t2 = (((((((-(t6 * t47 * 2.0) + t14 * 4.0) + t18 * 8.0) + t4 * 12.0) -
           t2 * t47 * t52 * 12.0) -
          t5 * t47 * t54 * 30.0) -
         t6 * t50 * (t48 - 1.0) * 2.0) -
        t2 * t50 * t53 * 4.0) -
       t5 * t50 * t55 * 6.0;
  t6 = t15 * t41;
  t3 = t16 * t41;
  t5 = t17 * t41;
  Dvol[0] = ((((t6 * t47 * -2.0 - t3 * t47 * t52 * 12.0) -
               t6 * t50 * (t48 - 1.0) * 2.0) -
              t5 * t47 * t54 * 30.0) -
             t3 * t50 * t53 * 4.0) -
            t5 * t50 * t55 * 6.0;
  Dvol[1] = t2;
  Dvol[2] = t19;
  Dvol[3] = t223;
  Dvol[4] = t226;
  Dvol[5] = t222;
  Dvol[6] = t2;
  t6 = t15 * t40;
  t3 = t16 * t40;
  t5 = t17 * t40;
  Dvol[7] = ((((t6 * t47 * -2.0 - t3 * t47 * t52 * 12.0) -
               t6 * t50 * (t48 - 1.0) * 2.0) -
              t5 * t47 * t54 * 30.0) -
             t3 * t50 * t53 * 4.0) -
            t5 * t50 * t55 * 6.0;
  Dvol[8] = t36;
  Dvol[9] = t221;
  Dvol[10] = t220;
  Dvol[11] = t225;
  Dvol[12] = t19;
  Dvol[13] = t36;
  t6 = t15 * t39;
  t3 = t16 * t39;
  t5 = t17 * t39;
  Dvol[14] = ((((t6 * t47 * -2.0 - t3 * t47 * t52 * 12.0) -
                t6 * t50 * (t48 - 1.0) * 2.0) -
               t5 * t47 * t54 * 30.0) -
              t3 * t50 * t53 * 4.0) -
             t5 * t50 * t55 * 6.0;
  Dvol[15] = t224;
  Dvol[16] = t218;
  Dvol[17] = t219;
  Dvol[18] = t223;
  Dvol[19] = t221;
  Dvol[20] = t224;
  t6 = t15 * t44;
  t3 = t16 * t44;
  t5 = t17 * t44;
  Dvol[21] = (((((((t6 * t47 * -2.0 - t3 * t47 * t52 * 12.0) -
                   t6 * t50 * (t48 - 1.0) * 2.0) -
                  t5 * t47 * t54 * 30.0) -
                 t3 * t50 * t53 * 4.0) -
                t5 * t50 * t55 * 6.0) -
               t14 * 2.0) -
              t18 * 4.0) -
             t4 * 6.0;
  Dvol[22] = t34;
  Dvol[23] = t229;
  Dvol[24] = t226;
  Dvol[25] = t220;
  Dvol[26] = t218;
  Dvol[27] = t34;
  t6 = t15 * t42;
  t3 = t16 * t42;
  t5 = t17 * t42;
  Dvol[28] = (((((((t6 * t47 * -2.0 - t3 * t47 * t52 * 12.0) -
                   t6 * t50 * (t48 - 1.0) * 2.0) -
                  t5 * t47 * t54 * 30.0) -
                 t3 * t50 * t53 * 4.0) -
                t5 * t50 * t55 * 6.0) -
               g_Svol_tmp * 2.0) -
              h_Svol_tmp * 4.0) -
             i_Svol_tmp * 6.0;
  Dvol[29] = t227;
  Dvol[30] = t222;
  Dvol[31] = t225;
  Dvol[32] = t219;
  Dvol[33] = t229;
  Dvol[34] = t227;
  t6 = t15 * t43;
  t3 = t16 * t43;
  t5 = t17 * t43;
  Dvol[35] = (((((((t6 * t47 * -2.0 - t3 * t47 * t52 * 12.0) -
                   t6 * t50 * (t48 - 1.0) * 2.0) -
                  t5 * t47 * t54 * 30.0) -
                 t3 * t50 * t53 * 4.0) -
                t5 * t50 * t55 * 6.0) -
               t20 * 2.0) -
              t21 * 4.0) -
             t33 * 6.0;
}

/* End of code generation (fnc.c) */
