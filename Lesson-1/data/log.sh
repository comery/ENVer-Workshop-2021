nucmer --maxmatch -t 2 -l 100 -c 500 -p mummer ref.mito.fa test.mito.fa
dnadiff -d mummer.delta -p mummer
show-aligns mummer.1delta Ref Scaffold1  > demo.alns
