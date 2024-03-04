rm runs_measurements/D_single_avg_list.txt
rm runs_measurements/D_single_avg.txt

mkdir -p data
mkdir -p plots
mkdir -p runs_measurements
mkdir -p runs_measurements/sdsam_plots

for i in {1..100}

do
    echo "Running iteration $i"
    bash single_run.sh
    python3 D_single.py
    sed -n '2p' data/D_single_value.txt >> runs_measurements/D_single_avg_list.txt
    python3 avg_D_single.py
    cp plots/single_disp_sampled.pdf runs_measurements/sdsam_plots/sdsam_$i.pdf
    echo "-----------------------------------"
done
