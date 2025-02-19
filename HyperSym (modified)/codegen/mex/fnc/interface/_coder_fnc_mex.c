/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_fnc_mex.c
 *
 * Code generation for function '_coder_fnc_mex'
 *
 */

/* Include files */
#include "_coder_fnc_mex.h"
#include "_coder_fnc_api.h"
#include "fnc_data.h"
#include "fnc_initialize.h"
#include "fnc_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void fnc_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
                     const mxArray *prhs[2])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        3, "fnc");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 3,
                        "fnc");
  }
  /* Call the function. */
  fnc_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&fnc_atexit);
  /* Module initialization. */
  fnc_initialize();
  /* Dispatch the entry-point. */
  fnc_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  fnc_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_fnc_mex.c) */
