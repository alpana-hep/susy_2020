#!/bin/sh
input_Scan=$1
anaExe="analyzeLightBSM"
anaArg=$2 #"signal" //signal models to consider
exeAtWorker="worker_SP.sh"
#dataSetType="T6ttZg"
#dataSetType="T5ttttZg"
dataSetType=$3 #"Summer16v3" // extentsion to find the files
#dataSetType="TChiNG"
#dataSetType="TChiWG"
#dataSetType="T5bbbbZg"
ext="Summer16v3"
#cp rootoutput_newbin/*v18.root . 
hist1=$4
#mkdir /eos/uscms/store/user/bkansal/myProduction/limits_rootout/T5bbbbZg/njet6/${hist1}
filesToTransfer="makeDatacard_SBins.C,higgsAnalysis.tar,FullRun2_WGJets_PhoIdloose_phopt40_BL.root, FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_BL.root, FullRun2_ZNuNu_PhoIdloose_phopt40_BL.root,FullRun2_TTJets_PhoIdloose_phopt40_BL.root,FullRun2_TTGJets_inc_PhoIdloose_phopt40_BL.root,FullRun2_GJets_QCD_PhoIdloose_phopt40_BL.root ,map_crosssection_SMprocess_v1.txt, map_crosssection_SMprocess.txt, ${anaExe}"
#, T5bbbbZg_MassScan.root" 
#PileupHistograms_2018_69mb_pm5.root,PileupHistograms_2016_69mb_pm5.root,PileupHistograms_2017_69mb_pm5.root,T5bbbbZg_MassScan.root"
#,TChiWG_MassScan.root,T5qqqqHg_Summer16v3_MassScan.root,TChiNG_MassScan.root,T5ttttZG_Summer16v3Fast_MassScan.root,T6ttZG_Summer16v3_MassScan.root"
#filesToTransfer="makeDatacard_SBins.C,higgsCombine.tar,GJets_v12.root,QCD_v12.root,TTGJets_v12.root,TTJetsHT_v12.root,WGJetsToLNuG_v12.root,WJetsToLNu_v12.root,ZJetsToNuNu_v12.root,ZGJetsToNuNuG_v12.root,ZGZJ_NuNuG_v12.root,${anaExe},PileupHistograms_0121_69p2mb_pm4p6.root,T5bbbbZg_MassScan.root"

while read -a massP
do 
    echo ${massP[0]} ${massP[1]}
    mass_p=1400
    if [ $(echo $anaArg | grep -c "T6tt") -gt 0 ]
    then
	mass_p= 800

    fi
    echo ${massP[0]} ${mass_p}
    if ((${massP[0]}<$mass_p))
    then
	echo "yes"
	#continue
    else
	echo "no"
	jdl_file="condor_${dataSetType}_${anaArg}_${massP[0]}_${massP[1]}_${hist1}_job.jdl"
	log_prefix="condor_${dataSetType}_${anaArg}_${massP[0]}_${massP[1]}_${hist1}_job"
	echo "universe = vanilla">$jdl_file
	echo "Executable = $exeAtWorker">>$jdl_file
	echo "Should_Transfer_Files = YES">>$jdl_file
	echo "WhenToTransferOutput = ON_EXIT_OR_EVICT">>$jdl_file
	echo "Transfer_Input_Files = ${filesToTransfer}">>$jdl_file
	echo "Output = ${log_prefix}.stdout">>$jdl_file
	echo "Error = ${log_prefix}.stderr">>$jdl_file
	echo "Log = ${log_prefix}.condor">>$jdl_file
	echo "notification = never">>$jdl_file
	echo "x509userproxy = $X509_USER_PROXY">>$jdl_file
	echo "use_x509userproxy = True">>$jdl_file
	echo "request_disk = 1000000" >>$jdl_file
	echo "request_cpus = 1" >>$jdl_file
	echo "request_memory = 4.0GB" >> $jdl_file
	echo "Arguments = ${anaExe} ${massP[0]} ${massP[1]} ${anaArg} ${dataSetType} ${hist1}">>$jdl_file
	echo "Queue">>$jdl_file
	condor_submit $jdl_file
    fi
	
done < $input_Scan
