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
//   'ufl_file': ['NeoHooke_Tetra.py'],
//   'verbosity': 30,
//   'visualise': False}


#pragma once

#include <ufcx.h>

#ifdef __cplusplus
extern "C" {
#endif

extern ufcx_finite_element element_89b9bb2bf2f707d0b69620f7bd9824d286e16307;

extern ufcx_finite_element element_728dbc04d4646a1f37abc835822781e1c9f53e42;

extern ufcx_dofmap dofmap_89b9bb2bf2f707d0b69620f7bd9824d286e16307;

extern ufcx_dofmap dofmap_728dbc04d4646a1f37abc835822781e1c9f53e42;

extern ufcx_integral integral_14686bb76ad4ae9f60f75351d34e791381e0363a;

extern ufcx_integral integral_1046b4ab5797ae637d804b0eb83b4ff7528d15c9;

extern ufcx_integral integral_fd11d64685f115e557a4befa30cbae4dff5c108c;

extern ufcx_integral integral_b243fbe418bdbb45713a0c4e95a9ccc43ab10d81;

extern ufcx_integral integral_bb7774d738b488d501bb35a56e6e5dd92a5e38f1;

extern ufcx_form form_0cd153591654f37b817c47ca218742006743cad0;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Tetra_F;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Tetra_F(const char* function_name);

extern ufcx_form form_15465d038a09a09c29aafdd0bcf36fdb162b878a;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Tetra_J;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Tetra_J(const char* function_name);

extern ufcx_form form_ed9003a3202cfbc2d7d63fb2f8c036a98b36fc4d;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Tetra_Pi;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Tetra_Pi(const char* function_name);

#ifdef __cplusplus
}
#endif
