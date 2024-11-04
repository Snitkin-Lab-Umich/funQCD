
cd /opt/interproscan/
python3 setup.py -f interproscan.properties
bash /opt/interproscan/interproscan.sh --input /opt/interproscan/test_all_appl.fasta --output-dir ~/ --tempdir ~/ -f tsv -dp

