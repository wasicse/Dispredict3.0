#! /bin/bash
# Number of parallel run

start=`date +%s%N` 

echo "Please Enter the number of parallel run with Docker container";
echo "The parallel run should be less than the number of protein sequcnecs in input fasta file."
echo "We have seen that 20 parallel run is the best for 64 core CPU";

read n
rm -rf temp
mkdir -p temp
mkdir -p ./temp/Parallelinputs

../.venv/bin/python ./splitFasta.py $n	

mkdir -p ./temp/ParallelOutputs
mkdir -p ./temp/logoutputs
for ((i=1;i<=$n;i++)); 
do 
	# your-unix-command-here
	echo $i
	docker stop dispredict_$i && docker rm dispredict_$i
done
sleep 5
echo "All containers are stopped"

for ((i=1;i<=$n;i++)); 
do 
	# your-unix-command-here
	echo $i
	docker stop dispredict_$i && docker rm dispredict_$i
	docker run -itd --name dispredict_$i wasicse/dispredict3.0:latest
	docker cp ./temp/Parallelinputs/processedinput_$i.fasta dispredict_$i:/opt/Dispredict3.0/example/sample.fasta
	sleep 2
	docker exec -it dispredict_$i /bin/bash -c "source /opt/Dispredict3.0/.venv/bin/activate && /opt/Dispredict3.0/.venv/bin/python /opt/Dispredict3.0/script/Dispredict3.0.py -f /opt/Dispredict3.0/example/sample.fasta -o /opt/Dispredict3.0/output/" | tee -a ./temp/logoutputs/output_$i.txt &

done
echo "Waiting for the output"
wait


echo "Output is ready"
for ((i=1;i<=$n;i++)); 
do 
	mkdir -p ./temp/ParallelOutputs/output_$i
	docker cp dispredict_$i:/opt/Dispredict3.0/output/ ./temp/ParallelOutputs/output_$i
	docker stop dispredict_$i && docker rm dispredict_$i
done
echo "Output is copied"

end=`date +%s%N` 

echo "Total time in ms for :" 
echo $(( ($end - $start) / 1000000 ))

mkdir -p ./outputs
echo "Merging the output"
for ((i=1;i<=$n;i++)); 
do 
	cat "./temp/ParallelOutputs/output_$i/output/sample_disPred.txt" >> "./outputs/sample_disPred_merged.txt" 
	cat "./temp/ParallelOutputs/output_$i/output/sample_fullydisPred.txt" >> "./outputs/sample_fullydisPred_merged.txt" 
done


