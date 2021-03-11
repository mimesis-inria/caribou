### Doxygen generation

The following packages are required:
- texlive-scheme-basic  texlive-epstopdf texlive-newunicodechar

Then, simply run (within `doc/`)
```shell
doxygen 
```

To update the configuration file, run
```shell
doxygen -u Doxyfile
```

### Sphinx

Make sure Doxygen ran well. Then, generate the intersphinx using:

```shell
./generate_doxygen_intersphinx.py
```

Upload the doxygen using:
```shell
./upload_doxygen.sh
```

Go inside the sphinx directory and run
```shell
make html
```