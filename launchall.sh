#!/bin/bash
cd memetique
./alternate_launch.sh ../simu_Damien/Simu_naive_2snp_0.25 2 X
mv results/ results_Damien_2snp_0.25/
./alternate_launch.sh ../simu_Damien/Simu_naive_3snp_0.25 3 X
mv results/ results_Damien_3snp_0.25/
cd ../path_relinking
./alternate_launch.sh ../simu_Damien/Simu_naive_3snp_0.25 3 X
mv results/ results_Damien_3snp_0.25/
