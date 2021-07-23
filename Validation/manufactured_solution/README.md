# Validation against a manufactured solution

Simply start
```shell
docker run -ti --rm -v $PWD:/opt/validation -w /opt/validation jnbrunet:caribou_validation python main.py
```