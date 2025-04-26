/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fnc.h
 *
 * Code generation for function 'fnc'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void fnc(const emlrtStack *sp, const real_T in1[9], const real_T in2[6],
         real_T Svol[6], real_T Dvol[36]);

/* End of code generation (fnc.h) */
