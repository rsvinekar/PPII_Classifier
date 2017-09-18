############################################################################################
#   Secondary structure classifier and collator, for PPII secondary structure              #
#   Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore        #
#   Inputs from Sathya Dev Unudurthi, NUS Singapore                                        #
#   Uses vmd                                                                               #
#   vmd -dispdev text -e PPII_Classify.vmd -args *.pdb                                     #
#   This file is included into PPII_Classify.vmd , and not run standalone                  #
#  License: to be decided according to Institute policies                                  #
############################################################################################
#XTLSTTR Classifier

# Call external program from here, if necessary. Below code only parses its output

#puts -nonewline "\nChecking for file <= $lc_pdbname.pd.ln "

set xtlfile "$Path/xtl/$lc_pdbname/$lc_pdbname.pd.ln"
for { set i 0 } { $i<$length } { incr i } {
		set xtl($i) "?"
}
set sel [atomselect $molid "protein or nucleic"]
set chains [$sel get chain]
set i 1
set prev_chain ""
foreach chain $chains {
	if { $chain ne $prev_chain } {
	set connect($i) $chain
	incr i
	}
	set prev_chain $chain
}
if { [ file exists $xtlfile ] } {
	set xtl_file [open $xtlfile "r"]
	set count_xtl 0
	set diff 0
	set cumulative 0
	while { [ gets $xtl_file line ] != -1 } {
		
		if { [ regexp {^\s+(\d+)\s+(\d+)([A-Z]*)\s+([A-Z]+)\s+([0-9]+)\s+([A-Za-z-])} $line -> c_seg c_res c_ins resnam int_num xtl_val ] } {
			if { $c_ins ne {} } {
				set sel2 [atomselect $molid "resid $c_res and chain $connect($c_seg) and insertion $c_ins"]
			} else {
						set sel2 [atomselect $molid "resid $c_res and chain $connect($c_seg)"]
				}
			set count_xtl [lindex [lsort -unique [$sel2 get residue]] 0]
			# We want only one configuration, if 2 or more are possible. ASER,BSER etc. only one is chosen.
			if { ![info exists xtl($count_xtl)] } {
				set nn [$sel2  num]
				puts "stopped at $count_xtl ;$c_seg ; $connect($c_seg) ; $c_res ; '$c_ins' ; $xtl_val ; $resnam ; $int_num $nn $length"
				puts "line is\n$line"
			}
			if { $count_xtl < $length } {
			# XTL program is idiot. It thinks DNA coords are proteins too. 
				if { $xtl($count_xtl) eq "?" } {
						set xtl($count_xtl) $xtl_val
				} else { 	set xtl($count_xtl) "$xtl($count_xtl) $xtl_val" }
				progress_indicate $count_xtl
				flush stdout
			}
		}
	}
	set xxXTL 1
	$sel delete
	$sel2 delete
} else {
	puts "file $lc_pdbname.pd.ln not found"
	set xxXTL 0
}

