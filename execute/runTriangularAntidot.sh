Ne="120"
Ncirc="120"
Nk="35"
latticetype="triangular"
meshfun="antidot"
#
echo "Submitting multiple jobs for a ${latticetype} 2D antidot plasmonic lattice"
#
funname="runDispersionAntidotLattice"
echo "Submitting the function ${funname}"
#
#Running the script
for adivdloop in {"6","5","4","3.5","3","2.5","2","1.5"}; do 
	command="${funname}(${adivdloop},${Ne},round(${Ncirc}*2/${adivdloop}),${Nk},'${latticetype}','${meshfun}',0)"
	logname="logs/dispfine_${latticetype}_${adivdloop}_Ne${Ne}.txt"
	echo "       command = ${command}"
	echo "       adivd = ${adivdloop}"
	nohup matlab -nodesktop -nosplash -noFigureWindows -singleCompThread -nodisplay -r ${command} > ${logname} &
	echo " "
done 