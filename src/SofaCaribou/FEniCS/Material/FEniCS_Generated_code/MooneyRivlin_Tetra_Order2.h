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
//   'ufl_file': ['MooneyRivlin_Hexa.py',
//                'MooneyRivlin_Hexa_Order2.py',
//                'MooneyRivlin_Tetra.py',
//                'MooneyRivlin_Tetra_Order2.py'],
//   'verbosity': 30,
//   'visualise': False}


#pragma once

#include <ufcx.h>

#ifdef __cplusplus
extern "C" {
#endif

extern ufcx_finite_element element_827e0731048b4eeeb15d8b84d7fc7208175924a4;

extern ufcx_finite_element element_6fd3c9ff1fb1dbc5751c582f00708d15637ca6f1;

extern ufcx_finite_element element_6a0881dbfb2bde0a242b4240501a9d5aaa7154a6;

extern ufcx_finite_element element_8fcbd9ae6633c6790d5e3c68f72323da1ff6d252;

extern ufcx_dofmap dofmap_827e0731048b4eeeb15d8b84d7fc7208175924a4;

extern ufcx_dofmap dofmap_6fd3c9ff1fb1dbc5751c582f00708d15637ca6f1;

extern ufcx_dofmap dofmap_6a0881dbfb2bde0a242b4240501a9d5aaa7154a6;

extern ufcx_dofmap dofmap_8fcbd9ae6633c6790d5e3c68f72323da1ff6d252;

extern ufcx_integral integral_6e3e85e62daae473a2377742e4b7af41af62b45c;

extern ufcx_integral integral_5ea67f9b81e5b5450c159aba1b49050b1d8401ad;

extern ufcx_integral integral_426ee0244c9057767842c86385bc9fa2169d0ed0;

extern ufcx_form form_951eaf57992fe366bea4a596f5f0cdf26ca51fbb;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_MooneyRivlin_Tetra_Order2_F;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_MooneyRivlin_Tetra_Order2_F(const char* function_name);

extern ufcx_form form_110629da7c8e72d210c0f1da5bec6258d20e1431;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_MooneyRivlin_Tetra_Order2_J;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_MooneyRivlin_Tetra_Order2_J(const char* function_name);

extern ufcx_form form_fa32b752075f9d0db44f959810c569e728756173;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_MooneyRivlin_Tetra_Order2_Pi;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_MooneyRivlin_Tetra_Order2_Pi(const char* function_name);

#ifdef __cplusplus
}
#endif
