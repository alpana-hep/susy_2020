python3 BDT_strong_TTJets.py -y all -s1 T5bbbbZg -b1 FullRun2  -m RandS -nt 300 -md 5 -n Equalweight_T5bbbbZg_Phopt40_ST300to1500_MET200_13variables -cuts 'ST<=1500 && ST>300 && MET>200'
python3 BDT_strong_TTJets.py -y all -s1 T5bbbbZg -b1 FullRun2  -m RandS -nt 300 -md 5 -n Equalweight_T5bbbbZg_Phopt40_ST300to1500_MET200_13variables -cuts 'ST>1500 && MET>200'
python3 BDT_strong_TTJets.py -y all -s1 T5bbbbZg -b1 FullRun2  -m RandS -nt 300 -md 5 -n Equalweight_T5bbbbZg_Phopt100_MET200_13variables -cuts 'PhoPt>100 && MET>200'
python3 BDT_strong_TTJets.py -y all -s1 T5bbbbZg -b1 FullRun2  -m RandS -nt 300 -md 5 -n Equalweight_T5bbbbZg_Phopt40to100_MET200_13variables -cuts 'PhoPt>40 && PhoPt<100 && MET>200'
