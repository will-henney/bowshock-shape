#!/bin/sh
#
# To run as an executable on a unix platform, do the following:
# chmod 775 wget_data.bat
# ./wget_data.bat
#
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RCAS_160_pixfrac10.0.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RCAS_70_pixfrac10.0.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RCAS-1_160_pixfrac10.0.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RCAS-all_160_pixfrac10.0.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RCAS-1_70_pixfrac10.0.mod.fits"
wget -x "https://irsa.ipac.caltech.edu:443/data/Herschel/MESS/images/RCAS-all_70_pixfrac10.0.mod.fits"
