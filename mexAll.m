% minFunc
fprintf('Compiling minFunc files...\n');
mex -outdir minFunc/compiled minFunc/mcholC.c
mex -outdir minFunc/compiled minFunc/lbfgsC.c
mex -outdir minFunc/compiled minFunc/lbfgsAddC.c
mex -outdir minFunc/compiled minFunc/lbfgsProdC.c

% KPM
fprintf('Compiling KPM files...\n');
mex -IKPM -outdir KPM/compiled KPM/max_mult.c

% DAGlearn
fprintf('Compiling DAGlearn files...\n');
mex -IDAGlearn/ancestorMatrix -outdir DAGlearn/ancestorMatrix/compiled DAGlearn/ancestorMatrix/ancestorMatrixAddC_InPlace.c
mex -IDAGlearn/ancestorMatrix -outdir DAGlearn/ancestorMatrix/compiled DAGlearn/ancestorMatrix/ancestorMatrixBuildC.c

% L1GeneralOverlapping Group
fprintf('Compiling L1GeneralGroup files...\n');
mex -outdir L1GeneralGroup/compiled L1GeneralGroup/mex/projectRandom2C.c
mex -outdir L1GeneralGroup/compiled L1GeneralGroup/mex/auxGroupLinfProjectC.c
mex -outdir L1GeneralGroup/compiled L1GeneralGroup/mex/auxGroupL2ProjectC.c

% L1GeneralOverlapping Group
fprintf('Compiling L1GeneralOverlapping Group files...\n');
mex -outdir L1GeneralOverlappingGroup/compiled L1GeneralOverlappingGroup/projectNDgroup_DykstraFastC.c

% UGM
fprintf('Compiling UGM files...\n');
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeEdgeVEC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ExactC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ExactC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ChainC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeClampedPotentialsC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ICMC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_GraphCutC.c
mex -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Sample_GibbsC.c
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

% LLM2
fprintf('Compiling LLM2 files...\n');
mex -outdir LLM2/compiled LLM2/mex/LLM2_inferC.c
mex -outdir LLM2/compiled LLM2/mex/LLM2_suffStatC.c
mex -outdir LLM2/compiled LLM2/mex/LLM2_pseudoC.c

% LLM
fprintf('Compiling LLM files...\n');
mex -outdir LLM/compiled LLM/mex/LLM_inferC.c
mex -outdir LLM/compiled LLM/mex/LLM_suffStatC.c
mex -outdir LLM/compiled LLM/mex/LLM_pseudoC.c

% UGMep
fprintf('Compiling UGMep files\n');
mex -outdir UGMep/compiled UGMep/mex/UGMep_Decode_ICMC.c
mex -outdir UGMep/compiled UGMep/mex/UGMep_EnergyC.c
mex -outdir UGMep/compiled UGMep/mex/UGMep_makeClampedEnergyC.c
mex -outdir UGMep/compiled UGMep/mex/UGMep_Decode_GraphCutC.c
mex -outdir UGMep/compiled UGMep/mex/UGMep_Decode_AlphaExpansionC.c
mex -outdir UGMep/compiled UGMep/mex/UGMep_Decode_ExpandShrinkC.c

% Ewout
fprintf('Compiling Ewout files...\n');
mex -IForeign/Ewout -outdir Foreign/Ewout/compiled Foreign/Ewout/projectBlockL1.c Foreign/Ewout/oneProjectorCore.c Foreign/Ewout/heap.c
mex -IForeign/Ewout -outdir Foreign/Ewout/compiled Foreign/Ewout/projectBlockL2.c

% Misc
mex -outdir misc/compiled misc/sampleDiscrete_cumsumC.c

% SAG
mex -outdir SAG/compiled SAG/mex/SGD_logistic.c -largeArrayDims
mex -outdir SAG/compiled SAG/mex/ASGD_logistic.c -largeArrayDims
mex -outdir SAG/compiled SAG/mex/PCD_logistic.c -largeArrayDims
mex -outdir SAG/compiled SAG/mex/DCA_logistic.c -largeArrayDims
mex -outdir SAG/compiled SAG/mex/SAG_logistic.c -largeArrayDims
mex -outdir SAG/compiled SAG/mex/SAGlineSearch_logistic.c -largeArrayDims
mex -outdir SAG/compiled SAG/mex/SAG_LipschitzLS_logistic.c -largeArrayDims
mex -outdir SAG/compiled SAG/mex/SGD_logistic_BLAS.c -largeArrayDims -lmwblas
mex -outdir SAG/compiled SAG/mex/ASGD_logistic_BLAS.c -largeArrayDims -lmwblas
mex -outdir SAG/compiled SAG/mex/PCD_logistic_BLAS.c -largeArrayDims -lmwblas
mex -outdir SAG/compiled SAG/mex/DCA_logistic_BLAS.c -largeArrayDims -lmwblas
mex -outdir SAG/compiled SAG/mex/SAG_logistic_BLAS.c -largeArrayDims -lmwblas
mex -outdir SAG/compiled SAG/mex/SAGlineSearch_logistic_BLAS.c -largeArrayDims -lmwblas
mex -outdir SAG/compiled SAG/mex/SAG_LipschitzLS_logistic_BLAS.c -largeArrayDims -lmwblas

