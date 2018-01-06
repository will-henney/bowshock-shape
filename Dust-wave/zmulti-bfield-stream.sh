angle=$1
for z in $(seq 0.0 0.1 4.5); do
    python dust-bfield-stream.py $angle $z 
    echo
done
