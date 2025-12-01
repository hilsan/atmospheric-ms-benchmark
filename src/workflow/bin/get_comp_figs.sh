#Scale QCxMS and QCxMS2 data
for dir in 0*; do python  $SCRIPTS_NEIMS/../features/scale_intensities.py -i $dir/*2/*2/peaks.csv -o $dir/scaled_qcxms2-xtb.csv;  python  $SCRIPTS_NEIMS/../features/scale_intensities.py -i $dir/*2/*w*/peaks.csv -o $dir/scaled_qcxms2-dft.csv;  python  $SCRIPTS_NEIMS/../features/scale_intensities.py -i $dir/QCxMS/MS-run/result.csv -o $dir/scaled_qcxms.csv; done

#for dir in 0*; do  python $SCRIPTS_NEIMS/../features/ NIST2MSP.py -i $dir/nist.txt; done

#Scale QCxMS data

for dir in 0*; do python  $SCRIPTS_NEIMS/../features/compare_spectra.py --simulated $dir/annotated.sdf --reference $dir/exp.msp --output $dir/exp2neims.png;  python  $SCRIPTS_NEIMS/../features/compare_spectra.py --simulated $dir/scaled_qcxms.csv --reference $dir/exp.msp --output $dir/exp2qcxms.png; python  $SCRIPTS_NEIMS/../features/compare_spectra.py --simulated $dir/scaled_qcxms2-xtb.csv --reference $dir/exp.msp --output $dir/exp2qcxms-xtb.png; python  $SCRIPTS_NEIMS/../features/compare_spectra.py --simulated $dir/scaled_qcxms2-dft.csv --reference $dir/exp.msp --output $dir/exp2qcxms-dft.png; done

#Scale QCxMS data, no max
for dir in 0*; do python $SCRIPTS_NEIMS/../features/compare_spectra_nomax.py --simulated $dir/annotated.sdf --reference $dir/exp.msp --output $dir/exp2neims_nomax.png; python $SCRIPTS_NEIMS/../features/compare_spectra_nomax.py --simulated $dir/scaled_qcxms.csv --reference $dir/exp.msp --output $dir/exp2qcxms_nomax.png; python $SCRIPTS_NEIMS/../features/compare_spectra_nomax.py --simulated $dir/scaled_qcxms2-xtb.csv --reference $dir/exp.msp --output $dir/exp2qcxms-xtb_nomax.png; python $SCRIPTS_NEIMS/../features/compare_spectra_nomax.py --simulated $dir/scaled_qcxms2-dft.csv --reference $dir/exp.msp --output $dir/exp2qcxms-dft_nomax.png; done

mkdir res-allpeaks
python $SCRIPTS_NEIMS/../features/plot_similarity_histograms.py     --dir . --output res-allpeaks --qcxms_flag --qcxms2_flag --qcxms2_dft_flag

mkdir res-nomax
python $SCRIPTS_NEIMS/../features/plot_similarity_histograms.py     --dir . --output res-nomax --qcxms_flag exp2qcxms_nomax_similarity_scores.csv -qcxms2_flag exp2qcxms2_nomax_similarity_scores.csv --qcxms2_dft_flag exp2qcxms2_dft_nomax_similarity_scores.csv --neims exp2neims_nomax_similarity_scores.cs

#4. #Plot all spectra from a calculation 
#python $SCRIPTS_NEIMS/../visualization/plot_stacked_MS.py --experimental $dir/nist.txt -neims $dir/annotated.sdf --qcxms  $dir/scaled_qcxms.csv --#qcxms2_xtb $dir/scaled_qcxms2-xtb.csv --qcxms2_dft $dir/scaled_qcxms2-dft.csv --output $dir/all.png
