#!/bin/sh
#
# To run as an executable on a unix platform, do the following:
# chmod 775 wget_data.bat
# ./wget_data.bat
#
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RSCL_160_pixfrac10.0.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RSCL_70_pixfrac10.0.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RScl-70_5_AFGL190.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RScl-160_8_AFGL190.mod.fits"
