#!/bin/bash

run_fits () {
	echo running fits on $1 with aor $2
	echo "exp"
	python Transit_fit_nnbr.py $2 $1 exp 10 true primary
	echo "none"
	python Transit_fit_nnbr.py $2 $1 none 10 true primary
	echo "slant"
	python Transit_fit_nnbr.py $2 $1 slant 10 true primary
}

echo "Generating fits for for WASP-52, WASP-62, WASP-63, WASP-67, WASP-79, WASP-101."

run_fits WASP_52 61517056
run_fits WASP_62 62167040 
run_fits WASP_62 62167552 
run_fits WASP_63 62168064
run_fits WASP_63 62168576
run_fits WASP_67 62169088
run_fits WASP_67 62169600
run_fits WASP_79 62173184
run_fits WASP_79 62173696
run_fits WASP_101 62158336
run_fits WASP_101 62159360