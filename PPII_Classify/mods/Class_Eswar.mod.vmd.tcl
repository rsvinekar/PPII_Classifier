############################################################################################
#   Secondary structure classifier and collator, for PPII secondary structure              #
#   Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore        #
#   Inputs from Sathya Dev Unudurthi, NUS Singapore                                        #
#   Uses vmd                                                                               #
#   vmd -dispdev text -e PPII_Classify.vmd -args *.pdb                                     #
#   This file is included into PPII_Classify.vmd , and not run standalone                  #
#  License: to be decided according to Institute policies                                  #
############################################################################################
# Stapeley method classifier
for { set i 0 } { $i < $length } { incr i } {
	set Class_NEshwar($i) "-"
}
set bk 0
set matching 0
for { set i 1 } { $i < [expr $length-1] } { incr i } {
	set matching 0
	set m_num 0
progress_indicate $i
	flush stdout
	for { set k $i } { $k < [expr $length-1] } { incr k } {
			#Check condition for average values.
	if { $phi_s($k) eq "-" || $psi_s($k) eq "-" } {
			break
		}
		if { [expr ($phi_s($k) > -180 && $phi_s($k) < -30) && (( $psi_s($k) > 60 && $psi_s($k) < 180 ) || ( $psi_s($k) > -180 && $psi_s($k) < -150 ) ) ] } {
			set m_num [expr $k-$i]
		} else {
			break
		}
	}
	if { $m_num > 3 } {

		for { set k $i } { $k < [expr $i+$m_num ] } { incr k } {
			if { [expr ($phi_s($k) > -90 && $phi_s($k) < -30) && (( $psi_s($k) > 60 && $psi_s($k) < 180 ) || ( $psi_s($k) > -180 && $psi_s($k) < -150 ) ) ] } {
				set Class_NEshwar($k) "P"
			} else {
				set Class_NEshwar($k) "E"
		}		
		
		
		}
		set $i [expr $i+$m_num]
	}
}
set xxESH 1
