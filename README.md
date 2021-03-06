# PPII_Classifier
Originally hosted on Dropbox https://www.dropbox.com/sh/698cy79q16ed6da/AAAChkFuAnxIGDUB3Vnp3tR8a?dl=0

## Secondary structure classifier and collator, for PPII secondary structure              
Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore
   Inputs from Sathya Dev Unudurthi, NUS Singapore and project Advisor(s)      
   Uses **vmd**
```
   vmd -dispdev text -e PPII_Classify.vmd -args *.pdb                                     
```
## README for PPII Classifier                                                             #
-----
Tested with VMD 1.8.7. Mileage may vary for other versions..

Usage :-
1. Create an xls directory (Writable) for output files
2. Create xtl and segno directories and put output of xtl and Segno there. (or modify the script to call the external programs)

### Batch mode - output to file
```
vmd -dispdev text -e PPII_Classify.vmd.tcl -args -a pdbs/*.pdb 2>&1 > output_file
```
### Interactive mode
```
vmd -dispdev text -e PPII_Classify.vmd.tcl -args pdbs/*.pdb
```

Tips:

- instead of *.pdb, you can be more specific. 1*.pdb , 1A*.pdb 1ARP.pdb etc.
- insert a `-a`  after -args, then specify pdb files, if you are redirecting output to a file. Significantly reduces file i/o and the output_file size.

Note : must be run from top level directory. to change, edit the Path variable in PPII_Classify.vmd

#### Filesystem structure
```
/
|
-- PPII_Classify.vmd.tcl  - Main wrapper script
-- debug_mesg - dumping scratchpad for some meaningless debugging output.
-- README - this readme.
=> pdbs  (input pdb files )
=> mods ( *.mod.vmd.tcl - files included in vmd script as modules, or as procedures )
    |
    --Preprocessor.mod.vmd.tcl : Extract and calculate average properties.
    --Preprocessor.proc.vmd.tcl - procedures for extracting data, dihedrals etc. called by above module.
    --VMD_read_pdb_ss.proc.vmd.tcl - procedure for extracting secondary structure from pdb file
    --Class_Adzhubei.mod.vmd.tcl : Adzhubei classifier. Implemented algorithm in VMD tcl code. Original definition, rigorous and slow.
    --Class_Eswar.mod.vmd.tcl : Eswar classifier. Implemented algorithm in VMD tcl code
    --Class_Stapeley.mod.vmd.tcl : Stapeley and Creamer classifier. Implemented algorithm in VMD tcl code
    --Class_Segno.mod.vmd.tcl : Segno classifier. Currently only parses existing SEGNO program output. Can be edited to call the external program
    --Class_Xtl.mod.vmd.tcl : XTLSTTR classifier. Currently only parses existing XTL program output. Can be edited to call the external program
=> logs ( messages pertaining to run of each file )
=> segno ( output of segno program. put them here in format <PDB ID>.segno )
=> xtl (output of xtlsttr program. Individual folders)
=> xls   ( output of this program. Tab separated text files with ending .xls . Will open in Excel as spreadsheet )

Format : Column headings self-explanatory.
#	Residue	Chain	Segment	RES	Phi	Psi	Diheco1/Zeta	Diheco2	Alpha	Tau	Dk	Avg_phi	Avg_psi	Avg_Diheco1	Avg_Diheco2	Stride	Crystal	Eswar	Stapeley	Regularity	Adzhubei	Segno	XTLSTTR	#

=> rejects (pdb files which for some reason cannot be handled. Manually dropped here)
```
