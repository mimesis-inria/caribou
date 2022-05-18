// This code conforms with the UFC specification version 2018.2.0.dev0
// and was automatically generated by FFCx version 0.4.3.dev0.
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

extern ufcx_finite_element element_1e6dd21838b5aed150f7a1df56fe135b84d2e70e;

extern ufcx_finite_element element_fa46853f91782f19b6f4c593f59dac7f0e916307;

extern ufcx_finite_element element_00e1a80bbdf056b4a6b4eb1b391c932e14a1ee65;

extern ufcx_finite_element element_0c4e046a76e54e3829221c3403f9b33c105e203b;

extern ufcx_dofmap dofmap_1e6dd21838b5aed150f7a1df56fe135b84d2e70e;

extern ufcx_dofmap dofmap_fa46853f91782f19b6f4c593f59dac7f0e916307;

extern ufcx_dofmap dofmap_00e1a80bbdf056b4a6b4eb1b391c932e14a1ee65;

extern ufcx_dofmap dofmap_0c4e046a76e54e3829221c3403f9b33c105e203b;

extern ufcx_integral integral_2ac95e85388a70078556da50ce6ef3b7590f7e32;

extern ufcx_integral integral_8f96590d15774b49a17c276ac3423edef1680aac;

extern ufcx_integral integral_6269b9f0375985f321a0515f6b9675d18f448abf;

extern ufcx_form form_b2d91c507b91412323d26984512c26287f493212;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_MooneyRivlin_Hexa_F;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_MooneyRivlin_Hexa_F(const char* function_name);

extern ufcx_form form_4a219084b8a702b21ea6681b69da89eb3b8e04c1;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_MooneyRivlin_Hexa_J;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_MooneyRivlin_Hexa_J(const char* function_name);

extern ufcx_form form_5d9041d8603a29c00e402c1879454bc112b35ea8;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_MooneyRivlin_Hexa_Pi;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_MooneyRivlin_Hexa_Pi(const char* function_name);

#ifdef __cplusplus
}
#endif
