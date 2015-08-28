
#!/bin/bash
echo "Running all"

channel="mu"
outDir="ratioFiles_${channel}/"

mkdir ${outDir}

Systematics=("MuEff_up" \
"MuEff_down" \
"musmear_up" \
"musmear_down" \
"Btag_up" \
"Btag_down" \
"EleFakeSF_up" \
"EleFakeSF_down" \
"JEC_up" \
"JEC_down" \
"JER_up" \
"JER_down" \
"PU_up" \
"PU_down" \
"QCD_up" \
"QCD_down" \
"otherMC_up" \
"otherMC_down" \
"ZJetsSF_up" \
"ZJetsSF_down" \
"pho_up" \
"pho_down" \
"toppt_up" \
"toppt_down")

# if [ ${channel}=="mu" ]; then
#     Systematics[0]="MuEff_up"
#     Systematics[1]="MuEff_down"
#     Systematics[2]="musmear_up"
#     Systematics[3]="musmear_down"
# fi

for i in `seq 0 23`; do
    echo "python makePlots.py ${channel} ${Systematics[i]} SaveOutput"
    python makePlots.py ${channel} ${Systematics[i]} SaveOutput
    echo "cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_${Systematics[i]}.root"
    cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_${Systematics[i]}.root
    cp templates_barrel_scaled.root ${outDir}templates_barrel_scaled_${Systematics[i]}.root
done

echo "python makePlots.py ${channel} SaveOutput"
python makePlots.py ${channel} SaveOutput

echo "cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_nominal.root"
cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_nominal.root
cp templates_barrel_scaled.root ${outDir}templates_barrel_scaled_nominal.root

cp TTGamma_SF_Lkhood.png TTGamma_SF_Lkhood_${channel}.png
cp Vgamma_SF_Lkhood.png Vgamma_SF_Lkhood_${channel}.png
cp jet_gamma_SF_Lkhood.png jet_gamma_SF_Lkhood_${channel}.png
