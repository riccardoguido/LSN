#!/bin/bash

for T in 2.0 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5
# for T in 2.5 0.5 

do
    echo "Running simulation at TEMP = $T"

    # Update the temperature value in the input file
    sed -i "s/^\(TEMP\s*\).*/\1$T/" ../INPUT/input.dat

    # Run the executable
    ./main.exe

    # Save results
    # Create T_* folder for each TEMP value
	mkdir -p ../T_${T}/CONFIG
	# Move .dat and configuration files from OUTPUT folder to the T_* folder
	mv ../OUTPUT/*.dat ../OUTPUT/*.out ../T_${T}/ 
	mv ../OUTPUT/CONFIG/config.spin ../T_${T}/CONFIG/ 
	# Copy configuration files in INPUT 
	cp ../T_${T}/CONFIG/config.spin ../INPUT/CONFIG/ 
done