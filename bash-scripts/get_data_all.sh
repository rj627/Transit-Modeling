#!/bin/bash

get_data () {
	echo outputting data to .txt
	python ../get_data.py $1 $2 exp
	python ../get_data.py $1 $2 none
	python ../get_data.py $1 $2 slant
	echo ---------------------------------
}

get_data WASP_52 61517056
get_data WASP_62 62167040 
get_data WASP_62 62167552 
get_data WASP_63 62168064
get_data WASP_63 62168576
get_data WASP_67 62169088
get_data WASP_67 62169600
get_data WASP_79 62173184
get_data WASP_79 62173696
get_data WASP_101 62158336
get_data WASP_101 62159360