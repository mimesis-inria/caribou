#/usr/bin/sh

# Change PORT, USER, DOMAIN, PATH
rsync -arvz -e 'ssh -p PORT' --progress --delete  Doxygen/html USER@DOMAIN:PATH
