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
for { set i $res_begin } { $i < $length } { incr i } {
	set Class_stapeley($i) "-"
}
set bk 0
set matching 0
for { set i $res_begin } { $i < $length } { incr i } {
	set matching 0
	set m_num 4
	set avg_phi_1 0
	set avg_psi_1 0
progress_indicate $i
	flush stdout
	if { $i >= [expr $length-4] } { continue }
	for { set j [expr $i+1] } { $j < [expr $i+4] } { incr j } {
		if { $phi_s($j) ne "-" } {
				set avg_phi_1 [expr $avg_phi_1 + $phi_s($j) ]
		} else {
				puts $logfile "missing phi $j in $pdbname - Class_stapeley"
		}
	}
	for { set j [expr $i] } { $j < [expr $i+3] } { incr j } {
		if { $psi_s($j) ne "-" } {
				set avg_psi_1 [expr $avg_psi_1 + $psi_s($j) ]
		} else {
				puts $logfile "missing psi $j in $pdbname - Class_stapeley"
		}
	}
	
	set avg_phi_1 [expr $avg_phi_1/3]
	set avg_psi_1 [expr $avg_psi_1/3]
	
	if { $avg_phi_1 >= -95 && $avg_phi_1 <= -55 && $avg_psi_1 >= 125 && $avg_psi_1 <= 165 } {
		set matching 1
		set avg_phi_2 $avg_phi_1
		set avg_psi_2 $avg_psi_1
		for { set j 4 } { $j < [expr $length - $i] } { incr j } {
			set avg_phi_1 0
			set avg_psi_1 0
			for { set k $i } { $k < [expr $j+$i] } { incr k } {
				if { $phi_s([expr $k+1]) ne "-" && $avg_phi_1 ne "-" } {
						set avg_phi_1 [expr $avg_phi_1+$phi_s([expr $k+1])]
				} else {
						set avg_phi_1 "-"
				}
				
				if { $psi_s($k) ne "-" && $avg_psi_1 ne "-" } {
						set avg_psi_1 [expr $avg_psi_1+$psi_s($k)]
				} else {
						set avg_psi_1 "-"
				}
			}
			if { $avg_phi_1 ne "-" } { set avg_phi_1 [expr $avg_phi_1/$j] }
			if { $avg_psi_1 ne "-" } { set avg_psi_1 [expr $avg_psi_1/$j] }
			if { $avg_phi_1 >= -95 && $avg_phi_1 <= -55 && $avg_psi_1 >= 125 && $avg_psi_1 <= 165 } {
				set m_num $j
			} else {
				break
			}
		}
	}
	if { $matching eq 1 } {
	#	set avg_phi_1 [expr $avg_phi_1 / [expr $m_num+6] ]
#		set avg_phi_1 [expr $avg_psi_1 / [expr $m_num+6] ]
		set diff 0
		set stddev("phi") 0
		set stddev("psi") 0
		for { set j $i } { $j < [expr $i+$m_num ] } { incr j } {
			if { $avg_phi_1 ne "-" } { set diff_avg("phi") [expr abs($phi_s([expr $j+1]) - $avg_phi_1)]
			} else { set diff_avg("phi") "-" }
			if { $avg_psi_1 ne "-" } { set diff_avg("psi") [expr abs($psi_s($j) - $avg_psi_1)]
			} else { set diff_avg("psi") "-" }
				
			if { $diff_avg("phi") <=20 &&  $diff_avg("psi") <=20 } {
				set diff 1
			}
		#	puts $logfile "came here $pdbname"
			if { $diff_avg("phi") ne "-" } { set stddev("phi") [expr $stddev("phi")+$diff_avg("phi")**2]
			} else { set stddev("phi") 1000000 }
			
			if { $diff_avg("psi") ne "-" } { set stddev("psi") [expr $stddev("psi")+$diff_avg("psi")**2]
			} else { set stddev("phi") 1000000 }
			# Put in a ridiculously high Numerical value for stddev, when - is encountered.
		}
		if { $diff eq 1 && [expr $stddev("phi")/$m_num] <=400 && [expr $stddev("phi")/$m_num] <=400 } {
			for { set j $i } { $j < [expr $i+$m_num ] } { incr j } {
				set Class_stapeley($j) "P"
			}
			set i [expr $i+$m_num]
		}
	}
}
set xxSTA 1
