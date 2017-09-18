#divide and rule :-) (^_^)
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/1A[0-9]*.pdb 2>&1 > logs/out1A1 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/1A[A-Z]*.pdb 2>&1 > logs/out1A2 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/1B*.pdb 2>&1 > logs/out1B &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/1C*.pdb 2>&1 > logs/out1C &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/1[D-G]*.pdb 2>&1 > logs/out2 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/1[H-Z]*.pdb 2>&1 > logs/out3 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/2*.pdb 2>&1 > logs/out4 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/3*.pdb 2>&1 > logs/out5 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/4*.pdb 2>&1 > logs/out6 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/5*.pdb 2>&1 > logs/out7 &
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/8*.pdb 2>&1 > logs/out8 &

