# How to generate NEIMS and QCxMS spectra with my scripts. 
    By Hilda Sandstroem, Jan. 2025. 

## Data prep
1. Make a dataset folder in NEIMS
2. Make a csv file with only one column which has the header 'SMILES' and contains all SMILES for the dataset (1 per row)
3. Run the make_sdf_from_smiles_batch.sh from src/file_mgmt folder with SMILES.csv as input (I think this needs the smiles to be in capital letters.

## Neims: 
1. Enter the dataset folder.
2. Run sbatch --array=0-9 submit_batch_neims_prediction.sh in the src/file_mgmt folder.
3. The produced spectrum for each compound will be in the compound folder (named by index) and called annotated_smiles.sdf
4. Run sbatch --array=0-9 submit_NEIMS_xtb.sh

## QCxMS
1. Enter the dataset folder
2. sbatch --array=0-9 submit_GS_MD.sh to run an xtb optimization and an MD run of that optimized GS structure. (1-2 hours)
3. Enter the individual compound folders/QCxMS/MS-run/, then run sbatch batch_hq.sh in order to run all fragmentation MD runs. (5 GB per compound)
        for dir in {0000..0009}; do cd $dir; cd QCxMS/MS-run; sh ../../../../../../src/file_mgmt/batch_hq.sh ; cd ../../../; done
6. Enter the TMP folder, remove all slurm scripts compile results from all fragmentation runs by running compile_freq.sh
8. Go one level above and run getres; ~/PlotMS.v.6.2.0/plotms ; python ../../../../../../src/features/compare_spectra.py --reference ../annotated.sdf --simulated result.csv --output 0001_10ps --similarity "cosine_weighted"
9. for dir in {0000..0001}; do cd $dir/QCxMS/MS-run/; sh ../../../../../../src/submit_frag_parallell.sh; cd ../../../; done

## Visualization of spectra

#Scale QCxMS and QCxMS2 data
for dir in 0*; do cp $dir/*2/*2/peaks.csv $dir/qcxms2-xtb.csv;  cp $dir/*2/*w*/peaks.csv $dir/qcxms2-dft.csv;  cp $dir/QCxMS/MS-run/result.csv $dir/qcxms.csv; done

#for dir in 0*; do  python $SCRIPTS_NEIMS/../features/ NIST2MSP.py -i $dir/nist.txt; done


for dir in 0*/; do python $SCRIPTS_NEIMS/../features/compare_spectra.py --simulated "$dir/annotated.sdf" "$dir/qcxms.csv" "$dir/qcxms2-xtb.csv" "$dir/qcxms2-dft.csv" --reference "$dir/exp.msp" --output "$dir"; done



mkdir res-allpeaks
python $SCRIPTS_NEIMS/../features/plot_similarity_histograms.py     --dir . --output res-allpeaks --qcxms_flag --qcxms2_flag --qcxms2_dft_flag

mkdir res-nomax
python $SCRIPTS_NEIMS/../features/plot_similarity_histograms.py     --dir . --output res-nomax --qcxms_flag exp2qcxms_nomax_similarity_scores.csv -qcxms2_flag exp2qcxms2_nomax_similarity_scores.csv --qcxms2_dft_flag exp2qcxms2_dft_nomax_similarity_scores.csv --neims exp2neims_nomax_similarity_scores.csv

#4. #Plot all spectra from a calculation 
#python $SCRIPTS_NEIMS/../visualization/plot_stacked_MS.py --experimental $dir/nist.txt -neims $dir/annotated.sdf --qcxms  $dir/scaled_qcxms.csv --#qcxms2_xtb $dir/scaled_qcxms2-xtb.csv --qcxms2_dft $dir/scaled_qcxms2-dft.csv --output $dir/all.png