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
//   'ufl_file': ['NeoHooke_Tetra_Order2.py'],
//   'verbosity': 30,
//   'visualise': False}


#pragma once

#include <ufcx.h>

#ifdef __cplusplus
extern "C" {
#endif

extern ufcx_finite_element element_016f4d74838397e54d7c6b0fe2c8708e1f8dc831;

extern ufcx_finite_element element_684c4478bdc3465842b882802e6cef80f9c0aa57;

extern ufcx_finite_element element_a0140b23b44cb727ba566623f0839d713a290934;

extern ufcx_finite_element element_ca645ec3d15899d79426639c0b13dcb0eb448682;

extern ufcx_dofmap dofmap_016f4d74838397e54d7c6b0fe2c8708e1f8dc831;

extern ufcx_dofmap dofmap_684c4478bdc3465842b882802e6cef80f9c0aa57;

extern ufcx_dofmap dofmap_a0140b23b44cb727ba566623f0839d713a290934;

extern ufcx_dofmap dofmap_ca645ec3d15899d79426639c0b13dcb0eb448682;

extern ufcx_integral integral_60e1719516a0645da5a232671d22f8e1d211e39e;

extern ufcx_integral integral_a077677632db3364f183cdc9c86773804c3922b1;

extern ufcx_integral integral_a999fed52bbef50742d6cb42e692b55c46125728;

extern ufcx_integral integral_f8c4f940d56452f481040e8f2823ffe0614cd434;

extern ufcx_integral integral_c70c60f391da6d166c73978fc0201ecce87c958b;

extern ufcx_form form_977810de326d0c14ebdaa9ebea9a7258ef295633;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Tetra_Order2_F;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Tetra_Order2_F(const char* function_name);

extern ufcx_form form_4a4ccd2ff9c3fa8a75b4a2d3b1610d94f4da028e;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Tetra_Order2_J;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Tetra_Order2_J(const char* function_name);

extern ufcx_form form_df9a5f818c9006602b042ed36c5b7429a72a8c08;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_NeoHooke_Tetra_Order2_Pi;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_NeoHooke_Tetra_Order2_Pi(const char* function_name);

#ifdef __cplusplus
}
#endif