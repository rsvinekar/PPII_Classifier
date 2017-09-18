############################################################################################
#   Secondary structure classifier and collator, for PPII secondary structure              #
#   Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore        #
#   Inputs from Sathya Dev Unudurthi, NUS Singapore                                        #
#   Uses vmd                                                                               #
#   vmd -dispdev text -e PPII_Classify.vmd -args *.pdb                                     #
#   This file is included into PPII_Classify.vmd , and not run standalone                  #
#  License: to be decided according to Institute policies                                  #
############################################################################################
set sel [atomselect $molfile "protein"]

set l [$sel num]
set res_begin 0
set sel2 [atomselect $molfile "index 0"]
set prev_res [$sel2 get residue]
set length 1

for { set i 1 } { $i < $l} { incr i } {
	set sel2 [atomselect $molfile "index $i"]
	set this_res [$sel2 get residue]
	if { $this_res > $prev_res } {
		incr length
	}
	set prev_res $this_res
}

set sel [atomselect $molfile "protein and name CA"]
set faux_len [$sel num]

set segid 0
set res_sel [atomselect $molfile "residue 0"]
set chain_ids(0) [ lindex [$res_sel get chain] 0 ]
$res_sel set chain $segid

for { set i 1 } { $i < $length} { incr i } {
	set prev [expr $i-1]
	set sel1 [atomselect $molfile "residue $prev and name N"]
	set sel2 [atomselect $molfile "residue $i and name CA"]
	if { [$sel1 num] > 0 && [$sel2 num] > 0} {
		set peplen [ measure bond [list [lindex [$sel1 get index] 0] [lindex [$sel2 get index] 0]] ]ls
		if { $peplen < 2 } {
		
			incr segid
		}
	} else {
		incr segid
	}

	set res_sel [atomselect $molfile "residue $i"]
	set chain_ids($i) [ lindex [$res_sel get chain] 0 ]
	$res_sel set chain $segid

}

set sel [atomselect $molfile "protein and name CA"]
set Sec(0) [$sel get structure]

if { $faux_len ne $length } {
		puts "Warning!!! Number of CA atoms = $faux_len < length."
		puts "Missing CAs. encountered errors will be marked with x"
		puts $local_log "Warning!!! Number of CA atoms = $faux_len < length."
		puts $local_log "Missing CAs. encountered errors will be marked with x"
}
puts -nonewline "Processing : "
for { set i 0 } { $i < $length } { incr i } {
	progress_indicate $i
	flush stdout
	extract_parameters $molfile 0 $i
#	puts -nonewline "$i list "
#save each of these values in arrays]
	set sel [atomselect $molfile "protein and name CA and residue $i"]
#	puts [$sel list]
	set res($i) [$sel get resname]
#	puts $res($i)
	set phi_s($i) $phi
	set psi_s($i) $psi
	set diheco1_s($i) $diheco1
	set diheco2_s($i) $diheco2
	set alpha_s($i) $alpha
	set tau_s($i) $tau
#set values in average parameters as "-". They will be populated later...
	set avg_phi($i) "-"
	set avg_psi($i) "-"
	set avg_diheco1($i) "-"
	set avg_diheco2($i) "-"
	set D_s($i) "-"
	set Sec($i) [$sel get structure]
	set SegId($i) [$sel get chain]
	set resid($i) [$sel get resid]
	set insertn($i) [$sel get insertion]
	if { ![ regexp {[A-Z]} $insertn($i) -> ] } {
		set insertn($i) ""
	}
	
		
#puts $i
}

for { set i 1 } { $i < [expr $length-1] } { incr i } {
		if { $phi_s([expr $i + 1 ]) ne "-" && $phi_s($i) ne "-" } {
			set delphi  [expr $phi_s([expr $i + 1 ]) - $phi_s($i)]
			set delphi  [expr $delphi*$delphi]
		} else {
			set delphi "-"
		}
		if { $psi_s([expr $i - 1 ]) ne "-" && $psi_s($i) ne "-" } {
			set delpsi  [expr $psi_s([expr $i - 1 ]) - $psi_s($i)]
			set delpsi  [expr $delpsi*$delpsi]
		} else {
			set delpsi "-"
		}
		if { $delphi ne "-" && $delpsi ne "-" } {
			set D_s($i) [expr $delphi+$delpsi]
			set D_s($i) [expr sqrt($D_s($i))]
		} else {
			set D_s($i) "-"
		}
	#	puts $i
}
## Needed. Over a window of 4-residues moving. -->
#      Average phi temp_sum()/temp_num() 0 
#      Average psi 1
#      Average diheco1 2
#      Average diheco2 3
#      D-value 4
for { set i 0 } { $i < [expr $length-3] } { incr i } {
	for { set j 0 } { $j < 5 } { incr j } {
		set temp_num($j) 0
		set temp_sum($j) 0
	}
	for { set j [expr $i] } { $j < [expr $i+3 ]  } { incr j } {
		if { $phi_s([expr $j]) ne "-" } {
			set temp_sum(0) [expr $temp_sum(0)+$phi_s($j)]
			incr temp_num(0)
		}
		if { $psi_s([expr $j+1]) ne "-" } {
			set temp_sum(1) [expr $temp_sum(1)+$psi_s([expr $j+1])]
			incr temp_num(1)
		}
		if { $diheco1_s($j) ne "-" } {
			set temp_sum(2) [expr $temp_sum(2)+$diheco1_s($j)]
			incr temp_num(2)
		}
		if { $diheco2_s($j) ne "-" } {
			set temp_sum(3) [expr $temp_sum(3)+$diheco2_s($j)]
			incr temp_num(3)
		}
	}
	if { $temp_num(0) > 2 } { set avg_phi($i) [expr $temp_sum(0)/$temp_num(0)] }
	if { $temp_num(1) > 2 } { set avg_psi($i) [expr $temp_sum(1)/$temp_num(1)] }
	if { $temp_num(2) > 2 } { set avg_diheco1($i) [expr $temp_sum(2)/$temp_num(2)] }
	if { $temp_num(3) > 2 } { set avg_diheco2($i) [expr $temp_sum(3)/$temp_num(3)] }
}
set xxPRE true


