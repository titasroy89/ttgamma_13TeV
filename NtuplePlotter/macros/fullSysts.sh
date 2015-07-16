
#!/bin/bash
echo "Running all"

channel="ele"
outDir="ratioFiles/"

mkdir ${outDir}

Systematics=("Btag_up" \
"Btag_down" \
"EleEff_up" \
"EleEff_down" \
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
"elesmear_up" \
"elesmear_down" \
"pho_up" \
"pho_down" \
"toppt_up" \
"toppt_down")

for i in `seq 0 23`; do
    echo "python makePlots.py ${channel} ${Systematics[i]}"
    python makePlots.py ${channel} ${Systematics[i]}
    echo "cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_${Systematics[i]}.root"
    cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_${Systematics[i]}.root
done

echo "python makePlots.py ${channel} > ${outDir}ratio_nominal.txt"
python makePlots.py ${channel} > ${outDir}ratio_nominal.txt

echo "cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_nominal.root"
cp templates_barrel_scaled_afterPhotonM3.root ${outDir}templates_barrel_scaled_afterPhotonM3_nominal.root

