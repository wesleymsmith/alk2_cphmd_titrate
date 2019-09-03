
echo "usage: fix_pdb.bash input.pdb > output.pdb"

reduce -trim $1 > temp.pdb

cat temp.pdb | \
sed "s/ CD  ILE/ CD1 ILE/" | \
sed "s/ HD1 ILE/HD11 ILE/" | \
sed "s/ HD2 ILE/HD12 ILE/" | \
sed "s/ HD3 ILE/HD13 ILE/" | \
sed "s/ HSD / HID /" | \
sed "s/ HSE / HIE /" | \
sed "s/ HT1 / H1  /" | \
sed "s/ HT2 / H2  /" | \
sed "s/ HT3 / H3  /" | \
sed "s/ OT1 / O   /" | \
sed "s/ OT2 / OXT /" 
