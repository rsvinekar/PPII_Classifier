############################################################################################
#   Secondary structure classifier and collator, for PPII secondary structure              #
#   Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore        #
#   Inputs from Sathya Dev Unudurthi, NUS Singapore                                        #
#   Uses vmd                                                                               #
#   vmd -dispdev text -e PPII_Classify.vmd -args *.pdb                                     #
#   This file is included into PPII_Classify.vmd , and not run standalone                  #
#  License: to be decided according to Institute policies                                  #
############################################################################################
#Segno Classifier

# Call external program from here, if necessary. Below code only parses its output

set segnofile "$pdbname.segno"
for { set i 0 } { $i<$length } { incr i } {
		set segno($i) "?"
	}
set segnofile "$Path/segno/$pdbname.segno"

if { [ file exists $segnofile ] } {
	set segno_file [open $segnofile "r"]
	set count_segno 0
	while { [ gets $segno_file line ] != -1 } {
		if { [ regexp {^\s+(\d+)\s+\w\w\w_(\w)\s+(\w)} $line -> c_seg seg_chain seg_val ] } {
		set segno($count_segno) $seg_val
		incr count_segno
		progress_indicate $count_segno
		flush stdout
		}
	}
	set xxSEG 1
} else {
	puts "$segnofile not found"
	set xxSEG 0
}
