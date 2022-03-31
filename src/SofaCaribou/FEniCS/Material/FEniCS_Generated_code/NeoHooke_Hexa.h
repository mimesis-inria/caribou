// This code conforms with the UFC specification version 2018.2.0.dev0
// and was automatically generated by FFCx version 0.3.1.dev0.
//
// This code was generated with the following parameters:
//
//  {'assume_aligned': -1,
//   'epsilon': 1e-14,
//   'output_directory': '../../FEniCS_Generated_code/',
//   'padlen': 1,
//   'profile': False,
//   'scalar_type': 'double',
//   'table_atol': 1e-09,
//   'table_rtol': 1e-06,
//   'tabulate_tensor_void': False,
//   'ufl_file': ['NeoHooke_Hexa.py'],
//   'verbosity': 30,
//   'visualise': False}


#pragma once

#include <ufcx.h>

#ifdef __cplusplus
extern "C" {
#endif

extern ufcx_finite_element element_41ef4b248fae7a876d8b10a862f220d1e7e12bc7;

extern ufcx_finite_element element_36a6e75bd43876b54b708372cec0c664cc09cec0;

extern ufcx_finite_element element_51053d16e9fef09070576a56209b34330c3e2930;

extern ufcx_finite_element element_ef1a1915f47c055271e11d5e56ee2bc57ebf935e;

extern ufcx_dofmap dofmap_41ef4b248fae7a876d8b10a862f220d1e7e12bc7;

extern ufcx_dofmap dofmap_36a6e75bd43876b54b708372cec0c664cc09cec0;

extern ufcx_dofmap dofmap_51053d16e9fef09070576a56209b34330c3e2930;

extern ufcx_dofmap dofmap_ef1a1915f47c055271e11d5e56ee2bc57ebf935e;

extern ufcx_integral integral_0872abcc642237aedf50191768b01dda57127af8;

extern ufcx_integral integral_6d00f7e6dd1add74319cd9d03d5ea6ad2d7d03a9;

extern ufcx_form form_3ffee4268a9a458c71af12bf104d0bc9ee2e76c8;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Hexa_F;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Hexa_F(const char* function_name);

extern ufcx_form form_cf557be0ce3c8b98a6b728de3d39ec3a45a52ab4;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Hexa_J;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Hexa_J(const char* function_name);

#ifdef __cplusplus
}
#endif