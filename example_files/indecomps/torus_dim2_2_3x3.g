# Dimension of Module: 7
LoadPackage("qpa");
start_read := Runtime();  # Record the starting time
Q := Quiver(9, [ [1,2, "E0"], [1,4, "E1"], [2,3, "E2"], [2,5, "E3"], [3,6, "E4"], [4,5, "E5"], [4,7, "E6"], [5,6, "E7"], [5,8, "E8"], [6,9, "E9"], [7,8, "E10"], [8,9, "E11"]]);
AQ := PathAlgebra(GF(2), Q);
MQ := RightModuleOverPathAlgebra(AQ, [2,2,0,2,1,0,0,0,0], [ ["E0", [[Z(2)^0,0*Z(2)],[0*Z(2),Z(2)^0]] ], ["E1", [[Z(2)^0,0*Z(2)],[0*Z(2),Z(2)^0]] ], ["E3", [[Z(2)^0],[0*Z(2)]] ], ["E5", [[Z(2)^0],[0*Z(2)]] ] ]);
