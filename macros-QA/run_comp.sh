for i in {1..50}
do
	cmd="root -b -q 'CompareToReactionPlane.C(\"/scratch/project_2003583/simO2_outputs/run_5p5TeV_midcent_job${i}/\")'"
	eval $cmd
done
