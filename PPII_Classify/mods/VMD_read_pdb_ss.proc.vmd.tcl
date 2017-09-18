#vmd_read_pdb_ss.tcl
# http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/vmd_use_pdb_ss/
#########################################################################################
#
#	REQUIREMENTS: VMD Version 1.0 or better
#
#	DESCRIPTION:
#	Normally VMD uses the program STRIDE in order to determine the
#	secondary structure of molecules.  STRIDE uses standard algorithms,
#	but its conclusions may on occasion differ slightly from
#	what you expect.  This script forces VMD to assign secondary
#	structure based on data in the appropriate record of the PDB file.
#	Thus you can fill in this record to obtain the desired secondary
#	structure representation.  To use this routine, you must
#	call vmd_use_pdb_ss with both the id of the molecule to define
#	secondary structure on and the file where the SS information can
#	be found.
#
#
#	PROCEDURES:
#	vmd_read_pdb_ss -- extracts appropriate fields of a given PDB file
#	vmd_use_pdb_ss -- determines the molecule to which a secondary
#	structure definition is to be applied, as well as the PDB
#	file from which records for the secondary structure are
#	to be obtained via vmd_read_pdb_ss
#
#	EXAMPLE USAGE:
#	Suppose the PDB file corresponding to the top molecule is called
#	PLV.PDB.  It is assumed that this PDB file has secondary structure
#	information contained within it.  Then the following command
#	will define secondary structure according to this information
#		vmd_use_pdb_ss top PLV.PDB
#	If you choose 'Cartoon' as a representation, you will see that
#	the secondary structure has been taken from the PDB file.
#
#
#	DOWNLOAD THE FILE:
#	vmd_use_pdb_ss.tcl
#
#
#	AUTHOR:
#	Andrew Dalke (dalke@ks.uiuc.edu)
#########################################################################################
# read the secondary structure records from a file and return
# the information in the form:
#    {SSType chain1 resid1 chain2 resid2} {... }
proc vmd_read_pdb_ss {pdb_filename} {
    set response {}
    # open the PDB file
    set infile [open $pdb_filename r]
    # read until the end of file
    while {[gets $infile line]} {
        set str [string range $line 0 5]
        if {$str == "HELIX "} {
            set ss helix
            set chain1 [string range $line 19 19]
            set resid1 [string range $line 21 24]
            set chain2 [string range $line 31 31]
            set resid2 [string range $line 33 36]
            lappend response [list $ss $chain1 $resid1 $chain2 $resid2]
            continue
        }
        if {$str == "SHEET "} {
            # get the needed fields
            set ss sheet
            set chain1 [string range $line 21 21]
            set resid1 [string range $line 22 25]
            set chain2 [string range $line 32 32]
            set resid2 [string range $line 33 36]
            lappend response [list $ss $chain1 $resid1 $chain2 $resid2]
            continue
        }

        # also, if read ATOM then there are no more def's
        if {$str == "ATOM  " || $str == ""} {
            break
        }
    }
    # close the file and return the list of info
    close $infile
    return $response
}

#### Now a driver to get info from this routine
# Return 0 if no information available
# Return 1 otherwise
proc vmd_use_pdb_ss {molid pdb_filename} {
    upvar crystal crystal length length
    # get the list of information
    set ssdata [vmd_read_pdb_ss $pdb_filename]
    for { set i 0 } { $i<$length } { incr i } {
		set crystal($i) "?"
	}
    # was there data?
    if {$ssdata == {}} {return 0}

#rlc's modification: first, reset everything to coil as a default
    set sel_all [atomselect $molid "all"]
    $sel_all set structure coil

    # Go through each of the element
    foreach ele $ssdata {
        lassign $ele ss chain1 resid1 chain2 resid2
        # if the chains are " ", don't use them
        if {$chain1 == " "} {
            set seltext "protein and (resid $resid1 to $resid2)"
        } else {
            set seltext "chain $chain1 and (resid $resid1 to $resid2)"
        }
        # get the selection/ make it the right structure/ free it
        set sel [atomselect $molid $seltext]
        $sel set structure $ss
        $sel delete
    }
    for { set i 0 } { $i<$length } { incr i } {
		set sel [atomselect $molid "protein and residue $i"]
		set crystal($i) [lindex [$sel get structure] 0]
	}
    # all done, so return that I did something
    return 1
}
