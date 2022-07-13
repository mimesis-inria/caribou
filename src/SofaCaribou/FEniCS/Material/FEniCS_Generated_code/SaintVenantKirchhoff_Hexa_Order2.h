// This code conforms with the UFC specification version 2018.2.0.dev0
// and was automatically generated by FFCx version 0.4.3.dev0.
//
// This code was generated with the following parameters:
//
//  {'assume_aligned': -1,
//   'epsilon': 1e-14,
//   'output_directory': '../../FEniCS_Generated_code',
//   'padlen': 1,
//   'profile': False,
//   'scalar_type': 'double',
//   'table_atol': 1e-09,
//   'table_rtol': 1e-06,
//   'tabulate_tensor_void': False,
//   'ufl_file': ['SaintVenantKirchhoff_Hexa_Order2.py'],
//   'verbosity': 30,
//   'visualise': False}


#pragma once

#include <ufcx.h>

#ifdef __cplusplus
extern "C" {
#endif

extern ufcx_finite_element element_f19ec6a3cf30b1196725fb845253bce17bdfb48f;

extern ufcx_finite_element element_49b6ba7218f9f526441b549725eec026a82b9f77;

extern ufcx_finite_element element_119a68fcd15ec88a83d4738ad0c462ee394c1941;

extern ufcx_finite_element element_d7d48a5d5278e8ec3dc7a0b3f220b601d788b79f;

extern ufcx_dofmap dofmap_f19ec6a3cf30b1196725fb845253bce17bdfb48f;

extern ufcx_dofmap dofmap_49b6ba7218f9f526441b549725eec026a82b9f77;

extern ufcx_dofmap dofmap_119a68fcd15ec88a83d4738ad0c462ee394c1941;

extern ufcx_dofmap dofmap_d7d48a5d5278e8ec3dc7a0b3f220b601d788b79f;

extern ufcx_integral integral_34cace2036031317c041c95e66a39cc7fa1a9b89;

extern ufcx_integral integral_240bb06eb4d40cc37701ff4c023b0873babe0663;

extern ufcx_integral integral_5236f7dcb4d4e9c2fad4cfe20caf3a6797169e5a;

extern ufcx_integral integral_85e7b97c1a05f80016cca97889ffee49fb938ba0;

extern ufcx_integral integral_94d1e56c25b184e0f02445a3b3b56fc081ae3f39;

extern ufcx_form form_b37b28731bb9770358fbd6ddd72002ea20d33656;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_SaintVenantKirchhoff_Hexa_Order2_F;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_SaintVenantKirchhoff_Hexa_Order2_F(const char* function_name);

extern ufcx_form form_06e6453dd2866682f0bc67264e5e1151ffbf09fb;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_SaintVenantKirchhoff_Hexa_Order2_J;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_SaintVenantKirchhoff_Hexa_Order2_J(const char* function_name);

extern ufcx_form form_37df8b5d6c32c125d63ef7222ec9673112affc6e;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_SaintVenantKirchhoff_Hexa_Order2_Pi;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_SaintVenantKirchhoff_Hexa_Order2_Pi(const char* function_name);

#ifdef __cplusplus
}
#endif
