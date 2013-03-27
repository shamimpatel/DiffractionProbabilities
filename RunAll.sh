ProcessArgs=""

for i in $@
do
	if [ $i == -ForceDiffraction ]
	then
		ProcessArgs="-ForceDiffraction"
	fi
done

./DiffractProbMain || exit 1
python ProcessScatterResults.py $ProcessArgs || exit 1