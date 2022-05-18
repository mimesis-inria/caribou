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
//   'ufl_file': ['Ogden_Hexa.py',
//                'Ogden_Hexa_Order2.py',
//                'Ogden_Tetra.py',
//                'Ogden_Tetra_Order2.py'],
//   'verbosity': 30,
//   'visualise': False}


#pragma once

#include <ufcx.h>

#ifdef __cplusplus
extern "C" {
#endif

extern ufcx_finite_element element_ab1105e07ed213f5e0deb12c29c4391da43095f5;

extern ufcx_finite_element element_c736860d03cdcbe7df10e7a9b33eb640b5ba3d75;

extern ufcx_finite_element element_88104a5f80830c04674a42596f57c3a5e41e9f4d;

extern ufcx_finite_element element_883d0563752c0ad245f0cfab5fa7d07355b5b8d1;

extern ufcx_dofmap dofmap_ab1105e07ed213f5e0deb12c29c4391da43095f5;

extern ufcx_dofmap dofmap_c736860d03cdcbe7df10e7a9b33eb640b5ba3d75;

extern ufcx_dofmap dofmap_88104a5f80830c04674a42596f57c3a5e41e9f4d;

extern ufcx_dofmap dofmap_883d0563752c0ad245f0cfab5fa7d07355b5b8d1;

extern ufcx_integral integral_aed899fb94a834eb2721b126ab236433ef84a0b2;

extern ufcx_integral integral_d0564fb57ad6bcdc38a42585ee442c5a2bbec34e;

extern ufcx_integral integral_c9b8c07843d62a401839b9b72057b9141731dd08;

extern ufcx_form form_0741b54c41f8021cbdb9eb28b3611630bff0da09;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_Ogden_Hexa_F;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_Ogden_Hexa_F(const char* function_name);

extern ufcx_form form_ea82f37706dc91573fafff895fdc3bf007efc84e;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_Ogden_Hexa_J;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_Ogden_Hexa_J(const char* function_name);

extern ufcx_form form_a0489a20b7d11e862f82cea7fb2906faa8e19bdd;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_Ogden_Hexa_Pi;

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_form_Ogden_Hexa_Pi(const char* function_name);

#ifdef __cplusplus
}
#endif
