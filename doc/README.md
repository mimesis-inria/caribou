# Building the documentation

To build the documentation, you can either install the required package yourself, 
or use a Docker container. The later can be easily install using
the `doc_builder.dockerfile` file. You can build the image using:
```shell
docker build -t caribou_doc_builder < doc_builder.dockerfile
```

Then, hop into the container using

```shell
docker run --rm -ti -v $CARIBOU_SRC:/w -w /w caribou_doc_builder 
```

Once you are inside the container, simply follow the remaining steps in the next sections.
### Doxygen generation

The following packages are required:
- texlive-scheme-basic  texlive-epstopdf texlive-newunicodechar

Then, run (within `doc/`)
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
python3 generate_doxygen_intersphinx.py
```

Upload the doxygen using:
```shell
./upload_doxygen.sh
```

Go inside the sphinx directory and run
```shell
make html
```