echo "Generating underlying expression level"
python explvl.py
echo "Genearting reads"
./gen_read_new_exp.sh
echo "Merging reads of DE and NON-DE genes"
./merge.sh
