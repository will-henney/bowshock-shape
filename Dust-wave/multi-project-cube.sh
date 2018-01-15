id=$1
inclinations="07 22 39 61"
azimuths="15 45 75"
for inc in $inclinations; do
    for azi in $azimuths; do
        python project-cube-inc-azi.py ${id}-cube-F $inc $azi
        echo
    done
done
