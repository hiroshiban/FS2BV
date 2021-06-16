**************************************************
README.txt on importing FreeSurfer files into BrainVoyager

Created    : "2017-08-27 16:24:46 ban"
Last Update: "2021-06-16 10:06:50 ban"
**************************************************

To import FreeSurfer recon_all's outputs (individual auto segmentation results e.g. *.mgz) into BrainVoyager,
please use

- ConvertFreeSurferAnnotation2BVpoi   : Converts FreeSurfer surface annotations to BrainVoyager POIs
- ConvertFreeSurferParcellation2BVvoi : Converts FreeSurfer MGZ parcellations to BrainVoyager VOIs
- ConvertFreeSurferMGZ2VMR            : Converts FreeSurer MGZ T1/ROI files to BrainVoayer VMRs
- ConvertFreeSurferRibbon2BL2VMR      : Converts FreeSurfer ribbon.mgz to BrainVoyager *_{LH|RH}_BL2.vmr
- ConvertFreeSurferSurf2SRF           : Converts FreeSurer surface files to BrainVoayer SRFs
- ImportFreeSurfer2BrainVoyager       : Imports FreeSurfer-processed files into BrainVoyager


To combine already generated BrainVoyager VMR with the FreeSurfer segmentation results for the surface
reconstructions on BrainVoyager, please try

- MaskVMRbyFreeSurferSegmentation       : for general masking purposes. Any *.mgz segmentation result can be used as a
                                          mask (by default, the parameters are tuned to process wm.seg.mgz as a mask)
- MaskVMRbyFreeSurferSegmentation_ribbon: Specific for applying a mask using the white (and gray) matter segmentation
                                          result in ribbon.mgz.
                                          Generally, for surface reconstructions, MaskVMRbyFreeSurferSegmentation_ribbon
                                          gives the better results.

if you want to import the FreeSurfer's gray/white matter segmentation results .

Furthermore, I have prepared additional tools listed below to import the ROI etc defined by the external software into BrainVoyager.
Please also try these functions if you are interested.

- ConvertNiftiRoi2BVvoi_ProbThres : Converts NII-format ROI probability map to BrainVoayer VOIs
                                    with thresholding the map values
- ConvertNiftiRoi2BVvoi_Labels    : Converts NII-format ROI probability map to BrainVoayer VOIs
                                    using the label lookuptable corresponding to the map ID
- ExtractFSLroi                   : Extracts specific value(s) from NII based on XML database
- ExtractFSLroiDirect             : Extracts specific value(s) from NII directly for ROI generations
- ConvertFSLroi2BVvoi             : Converts FSL NII ROIs to BrainVoyager VOIs
- ConvertSPMroi2BVvoi             : Converts SPM NII ROIs to BrainVoyager VOIs
- ConvertsAALroi2BVvoi            : Converts SPM AAL antomical tempolate (NII) ROIs to BrainVoyager VOIs

To find the VOIs in specific XYZs of TAL/MNI coords, please use the function below,

- GetAreaNameFromAtlasVOI         : Returns area candidates, in which the input XYZ coordinate(s)
                                    is(are) belonging to, based on the pre-defined VOI atlases.
