sbatch --array=1-100 job13.sh

for i in {31..37}
do
sbatch --array=1-22 job$i.sh
done

for i in {881043..881060}
do
scancel $i
done

sbatch --array=1-22 job13.sh
sbatch --array=1-22 job14.sh
sbatch --array=1-22 job15.sh
