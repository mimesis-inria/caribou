from distutils.core import setup, Extension

module1 = Extension('Helpers',
                    sources = ['helpers.cpp'])

setup (name = 'Helpers',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1])