############################################################################################
#   Secondary structure classifier and collator, for PPII secondary structure              #
#   Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore        #
#   Inputs from Sathya Dev Unudurthi, NUS Singapore                                        #
#   Uses vmd                                                                               #
#   vmd -dispdev text -e PPII_Classify.vmd -args *.pdb                                     #
#   This file is included into PPII_Classify.vmd , and not run standalone                  #
#  License: to be decided according to Institute policies                                  #
############################################################################################
#Adzhubei Classifier

for { set i $res_begin } { $i < $length } { incr i } {
	set Class_Adzhubei($i) "-"
	set Regularity($i) "-"
}
set bk 0
set matching 0
puts $logfile "Length $length"
for { set i $res_begin } { $i < [expr $length-1] } { incr i } {
	progress_indicate $i
	set matching 0
	set m_num 0
	set D_sum 0
	set D_count 0
	for { set j $i } { $j < [expr $i+4] } { incr j } {
		if { $j >= $length } { break }
		if { $D_s($j) eq "-" } {
				continue
		}
		
		set D_sum [expr $D_sum+$D_s($j)]
		incr D_count
		incr m_num
	}
	if { $D_count > 0 } {
		set D_sum [expr $D_sum/$D_count]
		} else { break }
	if { $D_sum <= 55 && $m_num > 3 } {
		set matching 1
		set m_num 4
		set D_sum 0
		set D_count 0
		for { set j 4 } { $j < [expr $length - $i] } { incr j } {
			for { set k $i } { $k < [expr $i+$j] } { incr k } {
				if { $D_s($k) eq "-" } {
						continue
				}
				set D_sum [expr $D_sum + $D_s($k)]
				incr D_count
				progress_indicate_stationary $j
			}
			set D_sum [expr $D_sum/$D_count]
			if { $D_sum < 55 } {
			
			} else {
				break
			}
		}
	}
	if { $matching eq 1 } {
		set pp2 0
		for { set j $i } { $j < [expr $i+$m_num-1 ] } { incr j } {
						set Regularity($j) "R"
				if { $alpha_s($j) eq "-" || $psi_s($j) eq "-" } {
						continue;
				}
				progress_indicate_stationary $j
				if { $alpha_s($j) >= -140 && $alpha_s($j) <= -80 && $psi_s($j) >= 90 && $psi_s($j) <= 180 } {
					set pp2 1
				} else {
					if { [expr $j-$i] < 3 } {
						set pp2 0
					}
					set $m_num [expr $j-$i]
					break
				}
		}
		if { $pp2 eq 1 } {
				for { set j $i } { $j < [expr $i+$m_num ] } { incr j } {
						set Class_Adzhubei($j) "P"
				}
				set $i [expr $i+$m_num]
		}
		
	}
}
set xxADZ 1
