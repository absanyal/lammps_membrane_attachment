mkdir -p runs
cd runs

for i in {1..1}
do
echo "submitting job $i"
    
mkdir -p run_$i
cd run_$i

cp ../../polymer.py .
cp ../../link_dist.py .
cp ../../cylindermath.py .
cp ../../angle.py .
cp ../../e2e_plots.py .

mkdir -p data
mkdir -p plots

python3 polymer.py
time mpirun -n 8 lmp_mpi -i input.lammps > out.run
python3 link_dist.py
python3 e2e_plots.py

echo "--------------------------------"

cd ..

done