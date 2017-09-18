############################################################################################
#   Secondary structure classifier and collator, for PPII secondary structure              #
#   Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore        #
#   Inputs from Sathya Dev Unudurthi, NUS Singapore                                        #
#   Uses vmd                                                                               #
#   vmd -dispdev text -e PPII_Classify.vmd -args *.pdb                                     #
#   This file is included into PPII_Classify.vmd , and not run standalone                  #
#  License: to be decided according to Institute policies                                  #
############################################################################################

proc progress_indicate_stationary { i } {
		upvar auto auto
		if { [info exists auto] && $auto == 1 } {
			;
		} else {
			set display [list "=" "-" "\\" "|" "/" "+"]
			set disvar [expr $i/6]
			puts -nonewline "\b"
			puts -nonewline [lindex $display [expr $disvar%6]]
			flush stdout
		}
}
proc progress_indicate { i } {
		upvar length length auto auto
		set div [expr $length/200]
	#	puts "$length $div"
		if { $div == 0 } { set div 1 }
		progress_indicate_stationary $i
		if { [expr $i%$div] == 0 } {
		#	puts $i
			flush stdout
			if { [expr $i%6] == 0 } { puts -nonewline "\b= " }
		}
	}

proc extract_parameters { index frm res} {
	upvar phi phi psi psi diheco1 diheco1 diheco2 diheco2 alpha alpha tau tau logfile logfile pdbname pdbname local_log local_log length length;
#-----------------------------------------------------Extract indices ----------------------------#
	#reset values of these variables. They will be carried forward by previous runs
	set phi "-"
	set psi "-"
	set diheco1 "-"
	set diheco2 "-"
	set alpha "-"
	set tau "-"
	set thischain [[atomselect $index "protein and residue $res and name CA"] get chain]
	set sel [atomselect $index "protein and name CA and residue $res"]
	set sel_measure [atomselect $index "protein and residue $res and name CA"]
	$sel_measure frame $frm
	set resname [$sel get resname]
	set resid [$sel get resid]
	set insertn [$sel get insertion]
	if { [$sel num] < 1 } { puts -nonewline "\bxx"; return 1 }
	set thischain [lindex [$sel get chain] 0];

	$sel delete;
	puts $local_log  "Line 17 extract_param $pdbname $res $resname $resid"
	flush $local_log 
	set prev [expr $res-1];
	set prev_flag 0
	set next [expr $res+1];
	set sel [atomselect $index "protein and name CA and residue $next"]

	set resnext [$sel get resname]
	set next_flag 0
	set nnxt [expr $res+2];
#chains
#	A number of PDB files have two or more configurations of atoms for the same type. i.e. ASER and BSER may be the same residue, with two C atoms, with occupancies < 1.
#	This means that two configurations of the same residue have been captured. We need to choose one of them
#	lindex makes sure that if a list of two or more atom indices is returned, only the first one will be chosen. If not handled this way, high-resolution crystal structures
# 	with alternate residue configurations fail, with a cryptic "ERROR: Molecule deleted (?) ". VMD's inbuilt phi and psi calculator cannot be used for same reason.

	set sel [atomselect $index "protein and name CA and chain $thischain"]
	if { $prev >= 0 } {
		set prevchain [[atomselect $index "protein and residue $prev and name CA"] get chain]
		#previous residue
		set atmpO  [lindex [[atomselect $index "protein  and residue $prev and name O"  ] get index ] 0]
		set atmpCA [lindex [[atomselect $index "protein  and residue $prev and name CA" ] get index ] 0]
		set atmpN  [lindex [[atomselect $index "protein  and residue $prev and name N"  ] get index ] 0]
		set atmpC  [lindex [[atomselect $index "protein  and residue $prev and name C"  ] get index ] 0]
	} else { set prevchain "-" }
	set atmpO_x [info exists atmpO]
	set atmpCA_x [info exists atmpCA]
	set atmpN_x [info exists atmpN]
	set atmpC_x [info exists atmpC]
	#this residue
	if { $res >= 0 } {
		set atmO  [lindex [[atomselect $index "protein  and residue $res and name O"  ] get index ] 0]
		set atmO_x [info exists atmO]
		set atmCA [lindex [[atomselect $index "protein  and residue $res and name CA" ] get index ] 0]
		set atmCA_x [info exists atmCA]
		set atmN  [lindex [[atomselect $index "protein  and residue $res and name N"  ] get index ] 0]
		set atmN_x [info exists atmN]
		set atmC  [lindex [[atomselect $index "protein  and residue $res and name C"  ] get index ] 0]
		set atmC_x [info exists atmC]
	}
	set atmO_x [info exists atmO]
	set atmCA_x [info exists atmCA]
	set atmN_x [info exists atmN]
	set atmC_x [info exists atmC]
#next residue
	if { $next < $length } {
		set nextchain [[atomselect $index "protein and residue $next and name CA"] get chain]
		set atmnO  [lindex [[atomselect $index "protein  and residue $next and name O"  ] get index ] 0]
		set atmnCA [lindex [[atomselect $index "protein  and residue $next and name CA" ] get index ] 0]
		set atmnN  [lindex [[atomselect $index "protein  and residue $next and name N"  ] get index ] 0]
		set atmnC  [lindex [[atomselect $index "protein  and residue $next and name C"  ] get index ] 0]
	}  else { set nextchain "-" }
	set atmnO_x [info exists atmnO]
	set atmnCA_x [info exists atmnCA]
	set atmnN_x [info exists atmnN]
	set atmnC_x [info exists atmnC]
	if { $nnxt < $length } {
	set nnxtchain [[atomselect $index "protein and residue $nnxt and name CA"] get chain]
		#next next residue
		set atmnnO  [lindex [[atomselect $index "protein  and residue $nnxt and name O"  ] get index ] 0]
		set atmnnCA [lindex [[atomselect $index "protein  and residue $nnxt and name CA" ] get index ] 0]
		set atmnnN  [lindex [[atomselect $index "protein  and residue $nnxt and name N"  ] get index ] 0]
		set atmnnC  [lindex [[atomselect $index "protein  and residue $nnxt and name C"  ] get index ] 0]
	} else { set nnxtchain "-" }
	set atmnnO_x [info exists atmnnO]
	set atmnnCA_x [info exists atmnnCA]
	set atmnnN_x [info exists atmnnN]
	set atmnnC_x [info exists atmnnC]
	if { $prevchain eq $thischain } { set prev_flag TRUE } else { set prev_flag FALSE }
	if { $nextchain eq $thischain } { set next_flag TRUE } else { set next_flag FALSE }
	if { $nnxtchain eq $thischain } { set nnxt_flag TRUE } else { set nnxt_flag FALSE }
#DIHEDRAL MEASURES
#Phi-psi
	if { $prev_flag } {
		if { $atmpC_x && $atmN_x && $atmCA_x && $atmC_x } {
			set phi [ measure dihed [ list $atmpC $atmN $atmCA $atmC ] ]
		}
		if { $atmpO_x && $atmpC_x && $atmC_x && $atmO_x } {
			set diheco1 [ measure dihed [ list  $atmpO $atmpC $atmC $atmO ] ]
		}
		if { $next_flag } {
			if { $atmpO_x && $atmpC_x && $atmnC_x && $atmnO } {
				set diheco2 [ measure dihed [ list $atmpO $atmpC $atmnC $atmnO ] ]
			}
			if { $atmpCA_x && $atmCA_x && $atmnCA } {
				set tau [ measure angle [ list $atmpCA $atmCA $atmnCA ] ]
			}
			if { $nnxt_flag && $atmpCA_x && $atmCA_x && $atmnCA_x && $atmnnCA_x } {
				set alpha [ measure dihed [ list $atmpCA $atmCA $atmnCA $atmnnCA ] ]
			}
		}
	}
	if { $next_flag && $atmN_x && $atmCA_x && $atmC && $atmnN } {
		set psi [ measure dihed [ list $atmN $atmCA $atmC $atmnN ] ]
	}

	$sel_measure delete
	return 0;

}
