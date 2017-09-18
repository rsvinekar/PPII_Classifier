#!/usr/local/bin/vmd -dispdev text -e
############################################################################################
#   Secondary structure classifier and collator, for PPII secondary structure              #
#   Author : Rithvik Vinekar. Written at Bioinformatics Institute, A*STAR Singapore        #
#   Inputs from Sathya Dev Unudurthi, NUS Singapore                                        #
#   Uses vmd                                                                               #
#   vmd -dispdev text -e PPII_Classify.vmd.tcl -args *.pdb                                 #
#  License: to be decided according to Institute policies                                  #
############################################################################################

set Path "."
#Non-interactive mode. Reduces file IO.
#                          -a flag may be added in automated scripts, which reduces redundant interactive feedback and significantly reduces output file size.
set logfile [open "debug_mesg" "w+"]
source "$Path/mods/Preprocessor.proc.vmd.tcl"
source "$Path/mods/VMD_read_pdb_ss.proc.vmd.tcl"
foreach filenam [lsort $argv] {
	if{ $filenam eq "-a" } {
		set auto 1
		continue
	}
	
	set molfile [mol new $filenam type {pdb} first 0 last -1 step 1 waitfor 1];
	set molid [mol new $filenam type {pdb} first 0 last -1 step 1 waitfor 1];
	set pdbname [string toupper [ lindex [split [lindex [split $filenam "/"] end] "."] 0]]
	mol ssrecalc $molfile
	set local_log [open "$Path/logs/$pdbname.log" "w+"]
	puts $logfile "Opened File $pdbname"
	set file_out "$Path/xls/$pdbname.xls"
	set fileId [open $file_out "w+"]
	flush stdout
	puts "Pdbname : $pdbname"
	set lc_pdbname [string tolower $pdbname]
	
#---------Progress indicators----------------------------
	puts ""
	puts "_______________________________________________ Begin Procedure _______________________________________________________"
	
	source "$Path/mods/Preprocessor.mod.vmd.tcl"
	puts "\b> #"
	flush stdout
	if { [ info exists xxPRE] } {
		puts -nonewline $fileId "#\t"
		puts -nonewline $fileId "Residue\t"
		puts -nonewline $fileId "Chain\t"
		puts -nonewline $fileId "Segment\t"
		puts -nonewline $fileId "RES\t"
		puts -nonewline $fileId "Phi\t"
		puts -nonewline $fileId "Psi\t"
		puts -nonewline $fileId "Diheco1/Zeta\t"
		puts -nonewline $fileId "Diheco2\t"
		puts -nonewline $fileId "Alpha\t"
		puts -nonewline $fileId "Tau\t"
		puts -nonewline $fileId "Dk\t"
		puts -nonewline $fileId "Avg_phi\t"
		puts -nonewline $fileId "Avg_psi\t"
		puts -nonewline $fileId "Avg_Diheco1\t"
		puts -nonewline $fileId "Avg_Diheco2\t"
	}
	#--------------------------------
	puts -nonewline $fileId "Stride\t"
	puts -nonewline $fileId "Crystal\t"
	vmd_use_pdb_ss $molid $filenam
	
	#--------------------------------
	puts -nonewline "Eswar      : "
	source "$Path/mods/Class_Eswar.mod.vmd.tcl"
	puts "\b> #"
	if { [ info exists xxESH] } {
		puts -nonewline $fileId "Eswar\t"
	}
	flush stdout
	#--------------------------------
	puts -nonewline "Stapeley   : "
	source "$Path/mods/Class_Stapeley.mod.vmd.tcl"
	puts "\b> #"
	if { [ info exists xxSTA] } {
		puts -nonewline $fileId "Stapeley\t"
	}
	flush stdout
	#--------------------------------
	puts -nonewline "Adzhubei   : "
	source "$Path/mods/Class_Adzhubei.mod.vmd.tcl"
	puts "\b> #"
	if { [ info exists xxADZ] } {
		puts -nonewline $fileId "Regularity\t"
		puts -nonewline $fileId "Adzhubei\t"
	}
	flush stdout
		#--------------------------------
	puts -nonewline "Segno      : "
	source "$Path/mods/Class_Segno.mod.vmd.tcl"
	puts "\b> #"
	if { [ info exists xxSEG] } {
		puts -nonewline $fileId "Segno\t"
	}
	flush stdout
		#--------------------------------
	puts -nonewline "XTLSTTR    : "
	source "$Path/mods/Class_Xtl.mod.vmd.tcl"
	puts "\b> #"
	if { [ info exists xxXTL] } {
		puts -nonewline $fileId "XTLSTTR\t"
	}
	flush stdout
		#--------------------------------
	puts -nonewline "Writing    : "
	puts $fileId "#"
	for { set i 0 } { $i < $length } { incr i } {
		if { [ info exists xxPRE] } {
			puts -nonewline $fileId "$i\t"
			puts -nonewline $fileId "$resid($i)$insertn($i)\t"
			puts -nonewline $fileId "$chain_ids($i)\t"
			puts -nonewline $fileId "$SegId($i)\t"
			puts -nonewline $fileId "$res($i)\t"
			puts -nonewline $fileId "$phi_s($i)\t"
			puts -nonewline $fileId "$psi_s($i)\t"
			puts -nonewline $fileId "$diheco1_s($i)\t"
			puts -nonewline $fileId "$diheco2_s($i)\t"
			puts -nonewline $fileId "$alpha_s($i)\t"
			puts -nonewline $fileId "$tau_s($i)\t"
			puts -nonewline $fileId "$D_s($i)\t"
			puts -nonewline $fileId "$avg_phi($i)\t"
			puts -nonewline $fileId "$avg_psi($i)\t"
			puts -nonewline $fileId "$avg_diheco1($i)\t"
			puts -nonewline $fileId "$avg_diheco2($i)\t"
		}
		puts -nonewline $fileId "$Sec($i)\t"
		puts -nonewline $fileId "$crystal($i)\t"
		if { [info exists xxESH] } {
			puts -nonewline $fileId "$Class_NEshwar($i)\t"
		}
		if { [ info exists xxSTA] } {
			puts -nonewline $fileId "$Class_stapeley($i)\t"
		}
		if { [ info exists xxADZ] } {
			puts -nonewline $fileId "$Regularity($i)\t"
			puts -nonewline $fileId "$Class_Adzhubei($i)\t"
		}
		if { [ info exists xxSEG] } {
			puts -nonewline $fileId "$segno($i)\t"
		}
		if { [ info exists xxXTL] } {
			puts -nonewline $fileId "$xtl($i)\t"
		}
		puts $fileId "#"
		flush $fileId
		progress_indicate $i
		flush stdout
	}
	puts "\b> #"
	puts "Output    => $file_out"
	puts "____________________________________________________ Done _____________________________________________________________"
	mol delete $molfile;
	mol delete $molid;
	flush $fileId
	close $fileId
	close $local_log
}
puts "____________________________________________________ All Done __________________________________________________________"
close $logfile
quit
		
	
	


