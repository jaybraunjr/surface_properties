# VMD TCL script to highlight strong interdigitating residues

# Load the structure and trajectory files (replace with actual filenames)
mol new "C:/Users/Jay/Desktop/modules/mlx_v2/6.6_2.gro"
mol addfile "C:/Users/Jay/Desktop/modules/mlx_v2/old_trajs/rep1_skip10.xtc" waitfor all

# List of strong interdigitating residues for each frame
set strong_residues {
    {311 343 350 360 392 562 588 655 730}
    {311 343 350 360 392 468 562 588 655 730}
    {311 350 360 392 468 588 655 730}
}

# Initialize a frame counter
set frame 0

# Loop over each frame and highlight the specified strong residues
foreach frame_residues $strong_residues {
    animate goto $frame
    incr frame

    if {[llength $frame_residues] > 0} {
        set selection [atomselect top "resid [join $frame_residues " or resid "]"]
        
        mol representation VDW
        mol selection "resid [join $frame_residues " or resid "]"
        mol material Opaque
        mol color ColorID 1
        mol addrep top

        $selection delete
    }
}

animate play
