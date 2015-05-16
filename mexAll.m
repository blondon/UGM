% minFunc
fprintf('Compiling minFunc files...\n');
mkdir minFunc/compiled
mex -outdir minFunc/compiled minFunc/mcholC.c
mex -outdir minFunc/compiled minFunc/lbfgsC.c
mex -outdir minFunc/compiled minFunc/lbfgsAddC.c
mex -outdir minFunc/compiled minFunc/lbfgsProdC.c

% UGM
fprintf('Compiling UGM files...\n');
mkdir UGM/compiled
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeEdgeVEC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ExactC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ExactC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ChainC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeClampedPotentialsC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ICMC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_GraphCutC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Sample_GibbsC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_CountBPC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_CountBPC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_MFC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_LBPC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_LBPC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_TRBPC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_TRBPC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_makePotentialsC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_PseudoNLLC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_LogConfigurationPotentialC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_AlphaExpansionC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_AlphaExpansionBetaShrinkC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_NLLC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ChainC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeCRFmapsC.c

% Misc
fprintf('Compiling misc files...\n');
mkdir misc/compiled
mex -outdir misc/compiled misc/sampleDiscrete_cumsumC.c
