#!/bin/bash

wget -O 919.zip "https://dl.acm.org/action/downloadSupplement?doi=10.1145%2F2168773.2168781&file=919.zip"
unzip 919.zip
mv 919 matlabsrc/919
rm 919.zip
