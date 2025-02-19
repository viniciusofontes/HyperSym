/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fnc_terminate.c
 *
 * Code generation for function 'fnc_terminate'
 *
 */

/* Include files */
#include "fnc_terminate.h"
#include "_coder_fnc_mex.h"
#include "fnc_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void fnc_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void fnc_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (fnc_terminate.c) */
