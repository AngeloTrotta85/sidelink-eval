#!/bin/bash

EXEC="/media/angelo/BigLinux/Programs/Eclipse/EclipseCPP/workspaces/Sidelink/sidelink-eval/Release/sidelink-eval"
INFOLDER="test12ott"


for CHAN in 2 3 4 5
do
	for DR in {100..500..100}
	do
		for ALG in 1 2
		do
			PROBLIMIT_SUM=0
			PROBNOLIMIT_SUM=0
			DELAY_SUM=0
			DIST_SUM=0
			
			for RUN in {1..5..1}
			do
				BASICFN="static_pos_10UAV_${CHAN}CH_${DR}DR_R${RUN}"
				FIN_FILENAME="${INFOLDER}/${BASICFN}.txt.lp-DmaxBASE-1000.000000-DmaxUAV-1000.000000-DmaxInterf-1200.000000-OF-${ALG}.out"
				FPOS_FILENAME="${INFOLDER}/${BASICFN}.txt"
				FOUT_FILENAME="${INFOLDER}/${BASICFN}.txt.lp-DmaxBASE-1000.000000-DmaxUAV-1000.000000-DmaxInterf-1200.000000-OF-${ALG}.out.ris"
				FLOG_FILENAME="${INFOLDER}/${BASICFN}.txt.lp-DmaxBASE-1000.000000-DmaxUAV-1000.000000-DmaxInterf-1200.000000-OF-${ALG}.out.log"
				
				#$EXEC -fin ${FIN_FILENAME} -fpos ${FPOS_FILENAME} -fout ${FOUT_FILENAME} > ${FLOG_FILENAME}
				
				RIS=`cat ${FOUT_FILENAME}`
				#VAL=`$(echo $RIS | tr ";" "\n")`
				#VAL=`cat ${FOUT_FILENAME} | tr ";" "\n"`
				#echo $RIS
				#IFS=';'; read -ra VAL <<< "$FOUT_FILENAME"`
				
				A="$(cut -d';' -f1 <<<"$RIS")"
				B="$(cut -d';' -f2 <<<"$RIS")"
				C="$(cut -d';' -f3 <<<"$RIS")"
				D="$(cut -d';' -f11 <<<"$RIS")"
				
				
				PROBLIMIT_SUM=`echo "print(${PROBLIMIT_SUM}+${A})" | python3`
				PROBNOLIMIT_SUM=`echo "print(${PROBNOLIMIT_SUM}+${B})" | python3`
				DELAY_SUM=`echo "print(${DELAY_SUM}+${C})" | python3`
				DIST_SUM=`echo "print(${DIST_SUM}+${D})" | python3`
				
				#echo "${PROBLIMIT_SUM}  ${PROBNOLIMIT_SUM}   ${DELAY_SUM}   "								
			done
			#echo 'print(1/3)' | python3
			PROBLIMIT_SUM=`echo "print(${PROBLIMIT_SUM}/5)" | python3`
			PROBNOLIMIT_SUM=`echo "print(${PROBNOLIMIT_SUM}/5)" | python3`
			DELAY_SUM=`echo "print(${DELAY_SUM}/5)" | python3`
			DIST_SUM=`echo "print(${DIST_SUM}/5)" | python3`
			
			#echo "AVG ${PROBLIMIT_SUM}  ${PROBNOLIMIT_SUM}   ${DELAY_SUM}   "
			
			BASICFN_VARC="${INFOLDER}/static_pos_10UAV_varCH_${DR}DR_A${ALG}.stat"
			BASICFN_VARD="${INFOLDER}/static_pos_10UAV_${CHAN}CH_varDR_A${ALG}.stat"
			
			echo "${CHAN} ${PROBLIMIT_SUM} ${PROBNOLIMIT_SUM} ${DELAY_SUM} ${DIST_SUM}" >> $BASICFN_VARC
			echo "${DR} ${PROBLIMIT_SUM} ${PROBNOLIMIT_SUM} ${DELAY_SUM} ${DIST_SUM}" >> $BASICFN_VARD
		done
	done
done
