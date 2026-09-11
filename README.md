Explaniation of the files
 - type_1
    -Iring - Simple first mesh
 - type_2 - Follower presures and more complex mesh with
 - type_2v2 - Using TPLOT comand for veryfication with analytical solution
 - type_3_ring - Contact problem
    Referance meshes (all have presure elements on outher edge)
    - IsecR.rev - inner mesh with a 1d elem (doesnt work)
    - IsecS.rev - presure elements on outher edge
    - IsecS2.rev - modified form feap forum - dealited paramieters and aditional data to main file
    - IsecS2L.rev - added a secound longtitude elements
    Main files
    -IsweepR - creating the frame elements ring (dont have contact)
    -IsweepS - creating the shell elements ring (with contact)
    -IsweepS2 - file from feap with consistant surf orientation - surfices point oposive directions
    -IsweepS3 - changes in surfice orientetion to pointing eachouther
    -IsweepSS - creating the inner solid mesh (with contact)

