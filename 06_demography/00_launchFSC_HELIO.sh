#!/bin/bash

#April 6 2020, script modified from Laurent Excoffier script on fsc26 manual
#The script will launch several jobs of fsc26 to estimate demographic parameters from the SFS
#It assumes the following structure of the observed sfs files:
#	scriptDir
#	|
#	|- - - *.est file
#	|- - - *.tpl file
#	|- - - fsc26
#	|- - - targetDir 
#		|
#		- *.obs files

fsc=/apps/fastsimcoal2/2.6.0.3/bin/fsc26
jobcount=0
msgs=conOutputs

#-------- Number of different runs per data set ------
numRuns=6
runBase=1
#-----------------------------

mkdir $msgs 2>/dev/null

#-------- Default run values ------
numSims=200000                 #-n command line option 
numCycles=50                   #-L command line option
#----------------------------------

#-------- Generic Name ------
genericName=HELIO
tplGenericName=HELIO_2orIM
estGenericName=HELIO_2orIM
#----------------------------------

for dirs in $genericName
do
	#Check that dirs is a directory
	if [ -d "$dirs" ]; then
		cd $dirs 
			echo "Main directory : $dirs"
			estFile=$estGenericName.est
			tplFile=$tplGenericName.tpl
			for (( runsDone=$runBase; runsDone<=$numRuns; runsDone++ ))
			do
				runDir="run$runsDone"
				mkdir $runDir 2>/dev/null
				echo "--------------------------------------------------------------------"
				echo ""
				echo "Currrent file: $subDirs $runDir"
				echo ""
				cd $runDir
				#Copying necessary files
				cp ../../$tplFile .
				cp ../../$estFile .
				cp ../*.obs .
				#Renaming files for consistency
				mv $tplFile ${genericName}.tpl
				mv $estFile ${genericName}.est
				let jobcount=jobcount+1
				jobName=${genericName}${jobcount}.slrm
				#Creating bash file on the fly
				( 
				echo "#!/bin/bash"
				echo ""	
				echo "#SBATCH --job-name=${jobName/.slrm/_TDIV2older}"
				echo "#SBATCH --cpus-per-task=1"
                                echo "#SBATCH --ntasks=48"
                                echo "#SBATCH --ntasks-per-node=12"
				echo "#SBATCH --mem=100M"
				echo "#SBATCH --partition=basic"
				echo "#SBATCH --mail-type=ALL"
				echo ""
				echo "module load fastsimcoal2"
				echo "for i in {1..10}"
				echo "do"
				echo "$fsc -t ${genericName}.tpl -n $numSims -e ${genericName}.est -m -M -L $numCycles --foldedSFS -c 48 -B 48 -0"
				echo "cat ./HELIO/HELIO.bestlhoods >> ./best_likelihoods"
				echo "mkdir simulated_sfs_\$i"
                                echo "mv ./HELIO/*txt ./simulated_sfs_\$i/"
                                echo "mv ./HELIO/HELIO_maxL.par ./simulated_sfs_\$i/"
                                echo "done"
				echo ""
				echo "echo \"Job $jobcount terminated\""
				) > $jobName
				chmod 755 $jobName
				echo "Bash file $jobName created"
				sbatch ./$jobName
				#./${jobName}
				cd .. #$runDir
			done
		cd .. #dirs
	fi
done	
