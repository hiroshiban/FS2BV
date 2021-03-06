function m = read_siemens_header(fname)


%
% read_siemens_header.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

fid = fopen(fname, 'r', 'b');
if fid < 0
	disp(['can''t open file ' fname]);
	return
end

fseek(fid, 0, 'bof');  m.h_G08_Ide_StudyDate_Year = fread(fid, 1, 'long');
fseek(fid, 4, 'bof');  m.h_G08_Ide_StudyDate_Month = fread(fid, 1, 'long');
fseek(fid, 8, 'bof');  m.h_G08_Ide_StudyDate_Day = fread(fid, 1, 'long');
fseek(fid, 12, 'bof');  m.h_G08_Ide_AcquisitionDate_Year = fread(fid, 1, 'long');
fseek(fid, 16, 'bof');  m.h_G08_Ide_AcquisitionDate_Month = fread(fid, 1, 'long');
fseek(fid, 20, 'bof');  m.h_G08_Ide_AcquisitionDate_Day = fread(fid, 1, 'long');
fseek(fid, 24, 'bof');  m.h_G08_Ide_ImageDate_Year = fread(fid, 1, 'long');
fseek(fid, 28, 'bof');  m.h_G08_Ide_ImageDate_Month = fread(fid, 1, 'long');
fseek(fid, 32, 'bof');  m.h_G08_Ide_ImageDate_Day = fread(fid, 1, 'long');
fseek(fid, 36, 'bof');  m.h_G08_Ide_StudyTime_Hour = fread(fid, 1, 'long');
fseek(fid, 40, 'bof');  m.h_G08_Ide_StudyTime_Minute = fread(fid, 1, 'long');
fseek(fid, 44, 'bof');  m.h_G08_Ide_StudyTime_Second = fread(fid, 1, 'long');
fseek(fid, 48, 'bof');  m.h_G08_Ide_StudyTime_Fraction = fread(fid, 1, 'long');
fseek(fid, 52, 'bof');  m.h_G08_Ide_AcquisitionTime_Hour = fread(fid, 1, 'long');
fseek(fid, 56, 'bof');  m.h_G08_Ide_AcquisitionTime_Minute = fread(fid, 1, 'long');
fseek(fid, 60, 'bof');  m.h_G08_Ide_AcquisitionTime_Second = fread(fid, 1, 'long');
fseek(fid, 64, 'bof');  m.h_G08_Ide_AcquisitionTime_Fraction = fread(fid, 1, 'long');
fseek(fid, 68, 'bof');  m.h_G08_Ide_ImageTime_Hour = fread(fid, 1, 'long');
fseek(fid, 72, 'bof');  m.h_G08_Ide_ImageTime_Minute = fread(fid, 1, 'long');
fseek(fid, 76, 'bof');  m.h_G08_Ide_ImageTime_Second = fread(fid, 1, 'long');
fseek(fid, 80, 'bof');  m.h_G08_Ide_ImageTime_Fraction = fread(fid, 1, 'long');
fseek(fid, 92, 'bof');  m.h_G08_Ide_Modality = fread(fid, 1, 'int');
fseek(fid, 96, 'bof');  t = fread(fid, 9, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_Manufacturer = char(t);
fseek(fid, 105, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_InstitutionID = char(t);
fseek(fid, 132, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_ReferringPhysician = char(t);
fseek(fid, 159, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_StationID = char(t);
fseek(fid, 186, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_ProcedureDescription_1 = char(t);
fseek(fid, 213, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_ProcedureDescription_2 = char(t);
fseek(fid, 240, 'bof');  t = fread(fid, 41, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_AdmittingDiagnosis = char(t);
fseek(fid, 281, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_ManufacturerModel = char(t);
fseek(fid, 308, 'bof');  t = fread(fid, 25, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_AttendingPhysician = char(t);
fseek(fid, 333, 'bof');  t = fread(fid, 25, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_Radiologist = char(t);
fseek(fid, 358, 'bof');  t = fread(fid, 25, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G08_Ide_OperatorIdentification = char(t);
fseek(fid, 396, 'bof');  m.h_G09_Ide_NumberOfMeasurements = fread(fid, 1, 'long');
fseek(fid, 404, 'bof');  m.h_G09_Ide_EvaluationMask = fread(fid, 1, 'long');
fseek(fid, 408, 'bof');  m.h_G09_Ide_Gap1212 = fread(fid, 1, 'long');
fseek(fid, 412, 'bof');  m.h_G09_Ide_LastMoveDate_Year = fread(fid, 1, 'long');
fseek(fid, 416, 'bof');  m.h_G09_Ide_LastMoveDate_Month = fread(fid, 1, 'long');
fseek(fid, 420, 'bof');  m.h_G09_Ide_LastMoveDate_Day = fread(fid, 1, 'long');
fseek(fid, 424, 'bof');  m.h_G09_Ide_LastMoveTime_Hour = fread(fid, 1, 'long');
fseek(fid, 428, 'bof');  m.h_G09_Ide_LastMoveTime_Minute = fread(fid, 1, 'long');
fseek(fid, 432, 'bof');  m.h_G09_Ide_LastMoveTime_Second = fread(fid, 1, 'long');
fseek(fid, 436, 'bof');  m.h_G09_Ide_LastMoveTime_Fraction = fread(fid, 1, 'long');
fseek(fid, 440, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_GeneratorIdentificationLabel = char(t);
fseek(fid, 467, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_GantryIdentificationLabel = char(t);
fseek(fid, 494, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_XRayTubeIdentificationLabel = char(t);
fseek(fid, 521, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_DetectorIdentificationLabel = char(t);
fseek(fid, 548, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_DASIdentificationLabel = char(t);
fseek(fid, 575, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_SMIIdentificationLabel = char(t);
fseek(fid, 602, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_CPUIdentificationLabel = char(t);
fseek(fid, 629, 'bof');  t = fread(fid, 9, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_HeaderVersion = char(t);
fseek(fid, 638, 'bof');  t = fread(fid, 25, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_OperatorIdentification = char(t);
fseek(fid, 663, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_Organ = char(t);
fseek(fid, 690, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_CaseIdentification = char(t);
fseek(fid, 717, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G09_Ide_RequestIdentification = char(t);
fseek(fid, 744, 'bof');  m.h_G09_Ide_LastEvaluation = fread(fid, 1, 'long');
fseek(fid, 748, 'bof');  m.h_G09_Ide_NoWithinScan = fread(fid, 1, 'long');
fseek(fid, 768, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G10_Pat_PatientName = char(t);
fseek(fid, 795, 'bof');  t = fread(fid, 13, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G10_Pat_PatientId = char(t);
fseek(fid, 808, 'bof');  m.h_G10_Pat_PatientBirthdate_Year = fread(fid, 1, 'long');
fseek(fid, 812, 'bof');  m.h_G10_Pat_PatientBirthdate_Month = fread(fid, 1, 'long');
fseek(fid, 816, 'bof');  m.h_G10_Pat_PatientBirthdate_Day = fread(fid, 1, 'long');
fseek(fid, 820, 'bof');  m.h_G10_Pat_PatientSex = fread(fid, 1, 'int');
fseek(fid, 824, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G10_Pat_PatientMaidenName = char(t);
fseek(fid, 851, 'bof');  t = fread(fid, 5, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G10_Pat_PatientAge = char(t);
fseek(fid, 856, 'bof');  m.h_G10_Pat_PatientSize = fread(fid, 1, 'double');
fseek(fid, 864, 'bof');  m.h_G10_Pat_PatientWeight = fread(fid, 1, 'long');
fseek(fid, 868, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G10_Pat_FrameOfReferenceUID = char(t);
fseek(fid, 933, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G10_Pat_SOPInstanceUID = char(t);
fseek(fid, 998, 'bof');  t = fread(fid, 17, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G10_Pat_AccessionNumber = char(t);
fseek(fid, 1024, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G11_Pat_Organ = char(t);
fseek(fid, 1052, 'bof');  m.h_G11_Pat_RegistrationDate_Year = fread(fid, 1, 'long');
fseek(fid, 1056, 'bof');  m.h_G11_Pat_RegistrationDate_Month = fread(fid, 1, 'long');
fseek(fid, 1060, 'bof');  m.h_G11_Pat_RegistrationDate_Day = fread(fid, 1, 'long');
fseek(fid, 1064, 'bof');  m.h_G11_Pat_RegistrationTime_Hour = fread(fid, 1, 'long');
fseek(fid, 1068, 'bof');  m.h_G11_Pat_RegistrationTime_Minute = fread(fid, 1, 'long');
fseek(fid, 1072, 'bof');  m.h_G11_Pat_RegistrationTime_Second = fread(fid, 1, 'long');
fseek(fid, 1076, 'bof');  m.h_G11_Pat_RegistrationTime_Fraction = fread(fid, 1, 'long');
fseek(fid, 1080, 'bof');  m.h_G11_Pat_UsedPatientWeight = fread(fid, 1, 'long');
fseek(fid, 1084, 'bof');  m.h_G11_Pat_OrganCode = fread(fid, 1, 'long');
fseek(fid, 1088, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G11_Pat_CaseIdentification = char(t);
fseek(fid, 1115, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G11_Pat_RequestIdentification = char(t);
fseek(fid, 1152, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G13_PatMod_ModifyingPhysician = char(t);
fseek(fid, 1180, 'bof');  m.h_G13_PatMod_ModificationDate_Year = fread(fid, 1, 'long');
fseek(fid, 1184, 'bof');  m.h_G13_PatMod_ModificationDate_Month = fread(fid, 1, 'long');
fseek(fid, 1188, 'bof');  m.h_G13_PatMod_ModificationDate_Day = fread(fid, 1, 'long');
fseek(fid, 1192, 'bof');  m.h_G13_PatMod_ModificationTime_Hour = fread(fid, 1, 'long');
fseek(fid, 1196, 'bof');  m.h_G13_PatMod_ModificationTime_Minute = fread(fid, 1, 'long');
fseek(fid, 1200, 'bof');  m.h_G13_PatMod_ModificationTime_Second = fread(fid, 1, 'long');
fseek(fid, 1204, 'bof');  m.h_G13_PatMod_ModificationTime_Fraction = fread(fid, 1, 'long');
fseek(fid, 1208, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G13_PatMod_PatientName_DICOM = char(t);
fseek(fid, 1273, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G13_PatMod_PatientId_DICOM = char(t);
fseek(fid, 1340, 'bof');  m.h_G13_PatMod_PatientBirthdate_Year = fread(fid, 1, 'long');
fseek(fid, 1344, 'bof');  m.h_G13_PatMod_PatientBirthdate_Month = fread(fid, 1, 'long');
fseek(fid, 1348, 'bof');  m.h_G13_PatMod_PatientBirthdate_Day = fread(fid, 1, 'long');
fseek(fid, 1352, 'bof');  m.h_G13_PatMod_PatientWeight = fread(fid, 1, 'long');
fseek(fid, 1356, 'bof');  m.h_G13_PatMod_PatientSex = fread(fid, 1, 'int');
fseek(fid, 1360, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G13_PatMod_Comment1 = char(t);
fseek(fid, 1387, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G13_PatMod_Comment2 = char(t);
fseek(fid, 1544, 'bof');  m.h_G18_Acq_SliceThickness = fread(fid, 1, 'double');
fseek(fid, 1552, 'bof');  m.h_G18_Acq_GeneratorVoltage = fread(fid, 1, 'long');
fseek(fid, 1556, 'bof');  m.h_G18_Acq_GeneratorVoltageDual = fread(fid, 1, 'long');
fseek(fid, 1560, 'bof');  m.h_G18_Acq_RepetitionTime = fread(fid, 1, 'double');
fseek(fid, 1568, 'bof');  m.h_G18_Acq_EchoTime = fread(fid, 1, 'double');
fseek(fid, 1576, 'bof');  m.h_G18_Acq_InversionTime = fread(fid, 1, 'double');
fseek(fid, 1584, 'bof');  m.h_G18_Acq_NumberOfAverages = fread(fid, 1, 'long');
fseek(fid, 1592, 'bof');  m.h_G18_Acq_ImagingFrequency = fread(fid, 1, 'double');
fseek(fid, 1604, 'bof');  m.h_G18_Acq_EchoNumber = fread(fid, 1, 'long');
fseek(fid, 1608, 'bof');  m.h_G18_Acq_DataCollectionDiameter = fread(fid, 1, 'long');
fseek(fid, 1612, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_DeviceSerialNumber = char(t);
fseek(fid, 1639, 'bof');  t = fread(fid, 9, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_SoftwareVersion = char(t);
fseek(fid, 1648, 'bof');  m.h_G18_Acq_DistanceSourceToDetector = fread(fid, 1, 'long');
fseek(fid, 1652, 'bof');  m.h_G18_Acq_DistanceSourceToPatient = fread(fid, 1, 'long');
fseek(fid, 1656, 'bof');  m.h_G18_Acq_GantryTilt = fread(fid, 1, 'long');
fseek(fid, 1660, 'bof');  m.h_G18_Acq_TableHeight = fread(fid, 1, 'long');
fseek(fid, 1668, 'bof');  m.h_G18_Acq_ExposureTime = fread(fid, 1, 'long');
fseek(fid, 1672, 'bof');  m.h_G18_Acq_Exposure = fread(fid, 1, 'long');
fseek(fid, 1676, 'bof');  t = fread(fid, 13, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_FilterIdLabel = char(t);
fseek(fid, 1696, 'bof');  m.h_G18_Acq_GeneratorPower = fread(fid, 1, 'double');
fseek(fid, 1704, 'bof');  m.h_G18_Acq_FocalSpot = fread(fid, 1, 'double');
fseek(fid, 1712, 'bof');  m.h_G18_Acq_CalibrationDate_Year = fread(fid, 1, 'long');
fseek(fid, 1716, 'bof');  m.h_G18_Acq_CalibrationDate_Month = fread(fid, 1, 'long');
fseek(fid, 1720, 'bof');  m.h_G18_Acq_CalibrationDate_Day = fread(fid, 1, 'long');
fseek(fid, 1724, 'bof');  m.h_G18_Acq_CalibrationTime_Hour = fread(fid, 1, 'long');
fseek(fid, 1728, 'bof');  m.h_G18_Acq_CalibrationTime_Minute = fread(fid, 1, 'long');
fseek(fid, 1732, 'bof');  m.h_G18_Acq_CalibrationTime_Second = fread(fid, 1, 'long');
fseek(fid, 1736, 'bof');  m.h_G18_Acq_CalibrationTime_Fraction = fread(fid, 1, 'long');
fseek(fid, 1740, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_ConvolutionKernel = char(t);
fseek(fid, 1767, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_ReceivingCoil = char(t);
fseek(fid, 1794, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_Gap1251 = char(t);
fseek(fid, 1828, 'bof');  t = fread(fid, 9, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_ImagedNucleus = char(t);
fseek(fid, 1837, 'bof');  t = fread(fid, 53, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G18_Acq_ProcDescr_new = char(t);
fseek(fid, 1920, 'bof');  m.h_G19_Acq1_CM_NetFrequency = fread(fid, 1, 'long');
fseek(fid, 1940, 'bof');  m.h_G19_Acq1_CM_Gap1040 = fread(fid, 1, 'long');
fseek(fid, 1944, 'bof');  m.h_G19_Acq1_CM_NoiseLevel = fread(fid, 1, 'long');
fseek(fid, 1948, 'bof');  m.h_G19_Acq1_CM_NumberOfDataBytes = fread(fid, 1, 'long');
fseek(fid, 1952, 'bof');  m.h_G19_Acq1_CM_ACScaleVector = fread(fid, 1, 'float');
fseek(fid, 1984, 'bof');  t = fread(fid, 33, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G19_Acq1_CM_ACElementConnected = char(t);
fseek(fid, 2048, 'bof');  m.h_G19_Acq2_Mr_TotalMeasurementTime = fread(fid, 1, 'double');
fseek(fid, 2056, 'bof');  m.h_G19_Acq2_Mr_StartDelayTime = fread(fid, 1, 'double');
fseek(fid, 2064, 'bof');  m.h_G19_Acq2_Mr_NumberOfPhases = fread(fid, 1, 'long');
fseek(fid, 2068, 'bof');  m.h_G19_Acq2_Mr_SequenceControlMask = fread(fid, 1, 'long');
fseek(fid, 2076, 'bof');  m.h_G19_Acq2_Mr_Gap1216 = fread(fid, 1, 'long');
fseek(fid, 2084, 'bof');  m.h_G19_Acq2_Mr_NumberOfFourierLinesNominal = fread(fid, 1, 'long');
fseek(fid, 2088, 'bof');  m.h_G19_Acq2_Mr_NumberOfFourierLinesCurrent = fread(fid, 1, 'long');
fseek(fid, 2092, 'bof');  m.h_G19_Acq2_Mr_NumberOfFourierLinesAfterZero = fread(fid, 1, 'long');
fseek(fid, 2096, 'bof');  m.h_G19_Acq2_Mr_FirstMeasuredFourierLine = fread(fid, 1, 'long');
fseek(fid, 2100, 'bof');  m.h_G19_Acq2_Mr_AcquisitionColumns = fread(fid, 1, 'long');
fseek(fid, 2104, 'bof');  m.h_G19_Acq2_Mr_ReconstructionColumns = fread(fid, 1, 'long');
fseek(fid, 2108, 'bof');  m.h_G19_Acq2_Mr_NumberOfAverages = fread(fid, 1, 'long');
fseek(fid, 2112, 'bof');  m.h_G19_Acq2_Mr_FlipAngle = fread(fid, 1, 'double');
fseek(fid, 2120, 'bof');  m.h_G19_Acq2_Mr_NumberOfPrescans = fread(fid, 1, 'long');
fseek(fid, 2240, 'bof');  m.h_G19_Acq2_Mr_NumberOfSaturationRegions = fread(fid, 1, 'long');
fseek(fid, 2248, 'bof');  m.h_G19_Acq2_Mr_ImageRotationAngle = fread(fid, 1, 'double');
fseek(fid, 2256, 'bof');  m.h_G19_Acq2_Mr_DwellTime = fread(fid, 1, 'double');
fseek(fid, 2264, 'bof');  m.h_G19_Acq2_Mr_CoilIdMask = fread(fid, 1, 'long');
fseek(fid, 2276, 'bof');  m.h_G19_Acq2_Mr_Gap1296 = fread(fid, 1, 'long');
fseek(fid, 2288, 'bof');  m.h_G19_Acq2_Mr_CoilPosition_Sag = fread(fid, 1, 'double');
fseek(fid, 2296, 'bof');  m.h_G19_Acq2_Mr_CoilPosition_Cor = fread(fid, 1, 'double');
fseek(fid, 2304, 'bof');  m.h_G19_Acq2_Mr_CoilPosition_Tra = fread(fid, 1, 'double');
fseek(fid, 2312, 'bof');  m.h_G19_Acq2_Mr_TotalMeasurementTimeCur = fread(fid, 1, 'double');
fseek(fid, 2320, 'bof');  m.h_G19_Acq2_Mr_MeasurementStatusMask = fread(fid, 1, 'long');
fseek(fid, 2344, 'bof');  m.h_G19_Acq2_Mr_ACADCPairNumber = fread(fid, 1, 'long');
fseek(fid, 2348, 'bof');  m.h_G19_Acq2_Mr_CoilClassMask = fread(fid, 1, 'long');
fseek(fid, 2356, 'bof');  m.h_G19_Acq2_Mr_ACCombinationMask = fread(fid, 1, 'long');
fseek(fid, 2360, 'bof');  m.h_G19_Acq2_Mr_ACElementNumber = fread(fid, 1, 'long');
fseek(fid, 2364, 'bof');  m.h_G19_Acq2_Mr_ACElementSelectMask = fread(fid, 1, 'long');
fseek(fid, 2368, 'bof');  m.h_G19_Acq2_Mr_ACElementDataMask = fread(fid, 1, 'long');
fseek(fid, 2372, 'bof');  m.h_G19_Acq2_Mr_ACElementToADCConnect = fread(fid, 1, 'float');
fseek(fid, 2436, 'bof');  m.h_G19_Acq2_Mr_ACElementNoiseLevel = fread(fid, 1, 'float');
fseek(fid, 2500, 'bof');  m.h_G19_Acq2_Mr_PhaseEncodingVectorSag = fread(fid, 1, 'float');
fseek(fid, 2524, 'bof');  m.h_G19_Acq2_Mr_ReadoutVectorSag = fread(fid, 1, 'float');
fseek(fid, 2548, 'bof');  m.h_G19_Acq2_Mr_EpiStimuMoniMode = fread(fid, 1, 'long');
fseek(fid, 2560, 'bof');  m.h_G19_Acq3_Mr_MagneticFieldStrength = fread(fid, 1, 'double');
fseek(fid, 2568, 'bof');  m.h_G19_Acq3_Mr_ADCVoltage = fread(fid, 1, 'double');
fseek(fid, 2576, 'bof');  m.h_G19_Acq3_Mr_TransmitterAmplitude = fread(fid, 1, 'double');
fseek(fid, 2584, 'bof');  m.h_G19_Acq3_Mr_NumberOfTransmitterAmplitudes = fread(fid, 1, 'long');
fseek(fid, 2592, 'bof');  m.h_G19_Acq3_Mr_TransmitterCalibration = fread(fid, 1, 'double');
fseek(fid, 2600, 'bof');  m.h_G19_Acq3_Mr_ReceiverTotalGain = fread(fid, 1, 'double');
fseek(fid, 2608, 'bof');  m.h_G19_Acq3_Mr_ReceiverAmplifierGain = fread(fid, 1, 'double');
fseek(fid, 2616, 'bof');  m.h_G19_Acq3_Mr_ReceiverPreamplifierGain = fread(fid, 1, 'double');
fseek(fid, 2624, 'bof');  m.h_G19_Acq3_Mr_ReceiverCableAttenuation = fread(fid, 1, 'double');
fseek(fid, 2632, 'bof');  m.h_G19_Acq3_Mr_ReconstructionScaleFactor = fread(fid, 1, 'double');
fseek(fid, 2640, 'bof');  m.h_G19_Acq3_Mr_PhaseGradientAmplitude = fread(fid, 1, 'double');
fseek(fid, 2648, 'bof');  m.h_G19_Acq3_Mr_ReadoutGradientAmplitude = fread(fid, 1, 'double');
fseek(fid, 2656, 'bof');  m.h_G19_Acq3_Mr_SelectionGradientAmplitude = fread(fid, 1, 'double');
fseek(fid, 2688, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G19_Acq3_Mr_SensitivityCorrectionLabel = char(t);
fseek(fid, 2720, 'bof');  m.h_G19_Acq3_Mr_ADCOffset = fread(fid, 1, 'double');
fseek(fid, 2736, 'bof');  m.h_G19_Acq3_Mr_TransmitterAttenuator = fread(fid, 1, 'double');
fseek(fid, 2744, 'bof');  m.h_G19_Acq3_Mr_TransmitterReference = fread(fid, 1, 'double');
fseek(fid, 2752, 'bof');  m.h_G19_Acq3_Mr_ReceiverReferenceGain = fread(fid, 1, 'double');
fseek(fid, 2760, 'bof');  m.h_G19_Acq3_Mr_ReceiverFilterFrequency = fread(fid, 1, 'long');
fseek(fid, 2768, 'bof');  m.h_G19_Acq3_Mr_ReferenceScaleFactor = fread(fid, 1, 'double');
fseek(fid, 2776, 'bof');  m.h_G19_Acq3_Mr_TotalGradientDelayTime = fread(fid, 1, 'double');
fseek(fid, 2784, 'bof');  m.h_G19_Acq3_Mr_RfWatchdogMask = fread(fid, 1, 'long');
fseek(fid, 2792, 'bof');  m.h_G19_Acq3_Mr_RfPowerErrorIndicator = fread(fid, 1, 'double');
fseek(fid, 2848, 'bof');  m.h_G19_Acq3_Mr_AdjustmentStatusMask = fread(fid, 1, 'long');
fseek(fid, 2852, 'bof');  m.h_G19_Acq3_Mr_FlowSensitivity = fread(fid, 1, 'float');
fseek(fid, 2860, 'bof');  m.h_G19_Acq3_Mr_FoVRatio = fread(fid, 1, 'float');
fseek(fid, 2864, 'bof');  m.h_G19_Acq3_Mr_BaseRawMatrixSize = fread(fid, 1, 'long');
fseek(fid, 2868, 'bof');  m.h_G19_Acq3_Mr_NumberOf2DPhaseOversamplingLines = fread(fid, 1, 'long');
fseek(fid, 2872, 'bof');  m.h_G19_Acq3_Mr_NumberOf3DPhaseOversamplingPart = fread(fid, 1, 'long');
fseek(fid, 2876, 'bof');  m.h_G19_Acq3_Mr_EchoLinePosition = fread(fid, 1, 'long');
fseek(fid, 2880, 'bof');  m.h_G19_Acq3_Mr_EchoColumnPosition = fread(fid, 1, 'long');
fseek(fid, 2884, 'bof');  m.h_G19_Acq3_Mr_LinesPerSegment = fread(fid, 1, 'long');
fseek(fid, 2892, 'bof');  m.h_G19_Acq3_Mr_PhaseEncodingVectorCor = fread(fid, 1, 'float');
fseek(fid, 2916, 'bof');  m.h_G19_Acq3_Mr_ReadoutVectorCor = fread(fid, 1, 'float');
fseek(fid, 2944, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G19_Acq4_CM_ParameterFileName = char(t);
fseek(fid, 3009, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G19_Acq4_CM_SequenceFileName = char(t);
fseek(fid, 3074, 'bof');  t = fread(fid, 9, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G19_Acq4_CM_SequenceFileOwner = char(t);
fseek(fid, 3089, 'bof');  t = fread(fid, 9, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G19_Acq4_CM_SequenceDescription = char(t);
fseek(fid, 3200, 'bof');  m.h_G20_Rel_Study = fread(fid, 1, 'long');
fseek(fid, 3204, 'bof');  m.h_G20_Rel_Series = fread(fid, 1, 'long');
fseek(fid, 3208, 'bof');  m.h_G20_Rel_Acquisition = fread(fid, 1, 'long');
fseek(fid, 3212, 'bof');  m.h_G20_Rel_Image = fread(fid, 1, 'long');
fseek(fid, 3216, 'bof');  m.h_G20_Rel_Gap0030 = fread(fid, 1, 'long');
fseek(fid, 3232, 'bof');  m.h_G20_Rel_Gap0035 = fread(fid, 1, 'double');
fseek(fid, 3280, 'bof');  m.h_G20_Rel_Location = fread(fid, 1, 'long');
fseek(fid, 3292, 'bof');  m.h_G20_Rel_AcquisitionsInSeries = fread(fid, 1, 'long');
fseek(fid, 3560, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G20_Rel_StudyInstanceUID = char(t);
fseek(fid, 3625, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G20_Rel_SeriesInstanceUID = char(t);
fseek(fid, 3690, 'bof');  t = fread(fid, 17, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G20_Rel_StudyId = char(t);
fseek(fid, 3712, 'bof');  m.h_G21_Rel1_CM_Gap1010 = fread(fid, 1, 'double');
fseek(fid, 3720, 'bof');  m.h_G21_Rel1_CM_Target_X = fread(fid, 1, 'double');
fseek(fid, 3728, 'bof');  m.h_G21_Rel1_CM_Target_Y = fread(fid, 1, 'double');
fseek(fid, 3736, 'bof');  m.h_G21_Rel1_CM_RoiMask = fread(fid, 1, 'short');
fseek(fid, 3744, 'bof');  m.h_G21_Rel1_CM_FoV_Height = fread(fid, 1, 'double');
fseek(fid, 3752, 'bof');  m.h_G21_Rel1_CM_FoV_Width = fread(fid, 1, 'double');
fseek(fid, 3768, 'bof');  m.h_G21_Rel1_CM_ImagePosition_Sag = fread(fid, 1, 'double');
fseek(fid, 3776, 'bof');  m.h_G21_Rel1_CM_ImagePosition_Cor = fread(fid, 1, 'double');
fseek(fid, 3784, 'bof');  m.h_G21_Rel1_CM_ImagePosition_Tra = fread(fid, 1, 'double');
fseek(fid, 3792, 'bof');  m.h_G21_Rel1_CM_ImageNormal_Sag = fread(fid, 1, 'double');
fseek(fid, 3800, 'bof');  m.h_G21_Rel1_CM_ImageNormal_Cor = fread(fid, 1, 'double');
fseek(fid, 3808, 'bof');  m.h_G21_Rel1_CM_ImageNormal_Tra = fread(fid, 1, 'double');
fseek(fid, 3816, 'bof');  m.h_G21_Rel1_CM_ImageDistance = fread(fid, 1, 'double');
fseek(fid, 3824, 'bof');  m.h_G21_Rel1_CM_ImagePositioningHistoryMask = fread(fid, 1, 'short');
fseek(fid, 3832, 'bof');  m.h_G21_Rel1_CM_ImageRow_Sag = fread(fid, 1, 'double');
fseek(fid, 3840, 'bof');  m.h_G21_Rel1_CM_ImageRow_Cor = fread(fid, 1, 'double');
fseek(fid, 3848, 'bof');  m.h_G21_Rel1_CM_ImageRow_Tra = fread(fid, 1, 'double');
fseek(fid, 3856, 'bof');  m.h_G21_Rel1_CM_ImageColumn_Sag = fread(fid, 1, 'double');
fseek(fid, 3864, 'bof');  m.h_G21_Rel1_CM_ImageColumn_Cor = fread(fid, 1, 'double');
fseek(fid, 3872, 'bof');  m.h_G21_Rel1_CM_ImageColumn_Tra = fread(fid, 1, 'double');
fseek(fid, 3904, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G21_Rel1_CM_StudyName = char(t);
fseek(fid, 3936, 'bof');  m.h_G21_Rel1_CM_ImageMagnificationFactor = fread(fid, 1, 'double');
fseek(fid, 3944, 'bof');  m.h_G21_Rel1_CM_ImageScrollOffset_Row = fread(fid, 1, 'double');
fseek(fid, 3952, 'bof');  m.h_G21_Rel1_CM_ImageScrollOffset_Col = fread(fid, 1, 'double');
fseek(fid, 3960, 'bof');  m.h_G21_Rel1_CM_ImagePixelOffset = fread(fid, 1, 'long');
fseek(fid, 3968, 'bof');  m.h_G21_Rel2_Mr_PhaseCorRowSeq = fread(fid, 1, 'long');
fseek(fid, 3972, 'bof');  m.h_G21_Rel2_Mr_PhaseCorColSeq = fread(fid, 1, 'long');
fseek(fid, 3976, 'bof');  m.h_G21_Rel2_Mr_PhaseCorRowRec = fread(fid, 1, 'long');
fseek(fid, 3980, 'bof');  m.h_G21_Rel2_Mr_PhaseCorColRec = fread(fid, 1, 'long');
fseek(fid, 3984, 'bof');  m.h_G21_Rel2_Mr_NumberOf3DRawPartNom = fread(fid, 1, 'long');
fseek(fid, 3988, 'bof');  m.h_G21_Rel2_Mr_NumberOf3DRawPartCur = fread(fid, 1, 'long');
fseek(fid, 3992, 'bof');  m.h_G21_Rel2_Mr_NumberOf3DImaPart = fread(fid, 1, 'long');
fseek(fid, 3996, 'bof');  m.h_G21_Rel2_Mr_Actual3DImaPartNumber = fread(fid, 1, 'long');
fseek(fid, 4000, 'bof');  m.h_G21_Rel2_Mr_Gap1338 = fread(fid, 1, 'long');
fseek(fid, 4004, 'bof');  m.h_G21_Rel2_Mr_NumberOfSlicesNom = fread(fid, 1, 'long');
fseek(fid, 4008, 'bof');  m.h_G21_Rel2_Mr_NumberOfSlicesCur = fread(fid, 1, 'long');
fseek(fid, 4012, 'bof');  m.h_G21_Rel2_Mr_CurrentSliceNumber = fread(fid, 1, 'long');
fseek(fid, 4016, 'bof');  m.h_G21_Rel2_Mr_CurrentGroupNumber = fread(fid, 1, 'long');
fseek(fid, 4020, 'bof');  m.h_G21_Rel2_Mr_MipStartRow = fread(fid, 1, 'long');
fseek(fid, 4024, 'bof');  m.h_G21_Rel2_Mr_MipStopRow = fread(fid, 1, 'long');
fseek(fid, 4028, 'bof');  m.h_G21_Rel2_Mr_MipStartColumn = fread(fid, 1, 'long');
fseek(fid, 4032, 'bof');  m.h_G21_Rel2_Mr_MipStopColumn = fread(fid, 1, 'long');
fseek(fid, 4036, 'bof');  m.h_G21_Rel2_Mr_MipStartSlice = fread(fid, 1, 'long');
fseek(fid, 4040, 'bof');  m.h_G21_Rel2_Mr_MipStopSlice = fread(fid, 1, 'long');
fseek(fid, 4044, 'bof');  m.h_G21_Rel2_Mr_SignalMask = fread(fid, 1, 'long');
fseek(fid, 4048, 'bof');  m.h_G21_Rel2_Mr_Gap1350 = fread(fid, 1, 'long');
fseek(fid, 4052, 'bof');  m.h_G21_Rel2_Mr_DelayAfterTrigger = fread(fid, 1, 'long');
fseek(fid, 4056, 'bof');  m.h_G21_Rel2_Mr_RRInterval = fread(fid, 1, 'long');
fseek(fid, 4064, 'bof');  m.h_G21_Rel2_Mr_NumberOfTriggerPulses = fread(fid, 1, 'double');
fseek(fid, 4072, 'bof');  m.h_G21_Rel2_Mr_RepetitionTime = fread(fid, 1, 'double');
fseek(fid, 4088, 'bof');  m.h_G21_Rel2_Mr_GateThreshold = fread(fid, 1, 'double');
fseek(fid, 4096, 'bof');  m.h_G21_Rel2_Mr_GateRatio = fread(fid, 1, 'double');
fseek(fid, 4104, 'bof');  m.h_G21_Rel2_Mr_NumberOfInterpolatedImages = fread(fid, 1, 'long');
fseek(fid, 4108, 'bof');  m.h_G21_Rel2_Mr_NumberOfEchoes = fread(fid, 1, 'long');
fseek(fid, 4112, 'bof');  m.h_G21_Rel2_Mr_SecondEchoTime = fread(fid, 1, 'double');
fseek(fid, 4120, 'bof');  m.h_G21_Rel2_Mr_SecondRepetitionTime = fread(fid, 1, 'double');
fseek(fid, 4128, 'bof');  m.h_G21_Rel2_Mr_CardiacCode = fread(fid, 1, 'long');
fseek(fid, 4136, 'bof');  m.h_G21_Rel2_Mr_CurrentSliceDistanceFactor = fread(fid, 1, 'double');
fseek(fid, 4152, 'bof');  m.h_G21_Rel2_Mr_SlabThickness = fread(fid, 1, 'double');
fseek(fid, 4160, 'bof');  m.h_G21_Rel2_Mr_PhaseEncodingVectorTra = fread(fid, 1, 'float');
fseek(fid, 4184, 'bof');  m.h_G21_Rel2_Mr_ReadoutVectorTra = fread(fid, 1, 'float');
fseek(fid, 4208, 'bof');  m.h_G21_Rel2_Mr_EpiChangeValueMagnitude = fread(fid, 1, 'float');
fseek(fid, 4212, 'bof');  m.h_G21_Rel2_Mr_EpiChangeValue_x_comp = fread(fid, 1, 'float');
fseek(fid, 4216, 'bof');  m.h_G21_Rel2_Mr_EpiChangeValue_y_comp = fread(fid, 1, 'float');
fseek(fid, 4220, 'bof');  m.h_G21_Rel2_Mr_EpiChangeValue_z_comp = fread(fid, 1, 'float');
fseek(fid, 4224, 'bof');  m.h_G21_Rel3_Mr_EpiReconstructionPhase = fread(fid, 1, 'double');
fseek(fid, 4232, 'bof');  m.h_G21_Rel3_Mr_EpiReconstructionSlope = fread(fid, 1, 'double');
fseek(fid, 4240, 'bof');  m.h_G21_Rel3_Mr_EpiCapacity = fread(fid, 1, 'double');
fseek(fid, 4288, 'bof');  m.h_G21_Rel3_Mr_EpiInductance = fread(fid, 1, 'double');
fseek(fid, 4312, 'bof');  m.h_G21_Rel3_Mr_EpiSwitchConfigurationCode = fread(fid, 1, 'long');
fseek(fid, 4324, 'bof');  m.h_G21_Rel3_Mr_EpiSwitchHardwareCode = fread(fid, 1, 'long');
fseek(fid, 4336, 'bof');  m.h_G21_Rel3_Mr_EpiSwitchDelayTime = fread(fid, 1, 'long');
fseek(fid, 4360, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G21_Rel3_Mr_EpiFileName = char(t);
fseek(fid, 4452, 'bof');  m.h_G21_Rel3_Mr_NormalVector_Sag = fread(fid, 1, 'float');
fseek(fid, 4456, 'bof');  m.h_G21_Rel3_Mr_NormalVector_Cor = fread(fid, 1, 'float');
fseek(fid, 4460, 'bof');  m.h_G21_Rel3_Mr_NormalVector_Tra = fread(fid, 1, 'float');
fseek(fid, 4524, 'bof');  m.h_G21_Rel3_Mr_PositionVector_Sag = fread(fid, 1, 'float');
fseek(fid, 4528, 'bof');  m.h_G21_Rel3_Mr_PositionVector_Cor = fread(fid, 1, 'float');
fseek(fid, 4532, 'bof');  m.h_G21_Rel3_Mr_PositionVector_Tra = fread(fid, 1, 'float');
fseek(fid, 4596, 'bof');  m.h_G21_Rel3_Mr_SaturationThickness = fread(fid, 1, 'float');
fseek(fid, 4620, 'bof');  m.h_G21_Rel3_Mr_SaturationWidth = fread(fid, 1, 'float');
fseek(fid, 4644, 'bof');  m.h_G21_Rel3_Mr_SaturationDistance = fread(fid, 1, 'float');
fseek(fid, 4668, 'bof');  m.h_G21_Rel3_Mr_ACADCOffset = fread(fid, 1, 'float');
fseek(fid, 4796, 'bof');  m.h_G21_Rel3_Mr_ACPreamplifierGain = fread(fid, 1, 'float');
fseek(fid, 4864, 'bof');  m.h_G21_Rel3_Mr_VectorSizeOriginal = fread(fid, 1, 'long');
fseek(fid, 4868, 'bof');  m.h_G21_Rel3_Mr_VectorSizeExtended = fread(fid, 1, 'long');
fseek(fid, 4872, 'bof');  m.h_G21_Rel3_Mr_AcquiredSpectralRange = fread(fid, 1, 'float');
fseek(fid, 4876, 'bof');  m.h_G21_Rel3_Mr_VOIPosition_Sag = fread(fid, 1, 'float');
fseek(fid, 4880, 'bof');  m.h_G21_Rel3_Mr_VOIPosition_Cor = fread(fid, 1, 'float');
fseek(fid, 4884, 'bof');  m.h_G21_Rel3_Mr_VOIPosition_Tra = fread(fid, 1, 'float');
fseek(fid, 4888, 'bof');  m.h_G21_Rel3_Mr_VOISize_Sag = fread(fid, 1, 'float');
fseek(fid, 4892, 'bof');  m.h_G21_Rel3_Mr_VOISize_Cor = fread(fid, 1, 'float');
fseek(fid, 4896, 'bof');  m.h_G21_Rel3_Mr_VOISize_Tra = fread(fid, 1, 'float');
fseek(fid, 4900, 'bof');  m.h_G21_Rel3_Mr_CSIMatrixSizeOriginal_Lines = fread(fid, 1, 'long');
fseek(fid, 4904, 'bof');  m.h_G21_Rel3_Mr_CSIMatrixSizeOriginal_Columns = fread(fid, 1, 'long');
fseek(fid, 4908, 'bof');  m.h_G21_Rel3_Mr_CSIMatrixSizeOriginal_Partitions = fread(fid, 1, 'long');
fseek(fid, 4912, 'bof');  m.h_G21_Rel3_Mr_CSIMatrixSizeExtended_Lines = fread(fid, 1, 'long');
fseek(fid, 4916, 'bof');  m.h_G21_Rel3_Mr_CSIMatrixSizeExtended_Columns = fread(fid, 1, 'long');
fseek(fid, 4920, 'bof');  m.h_G21_Rel3_Mr_CSIMatrixSizeExtended_Partitions = fread(fid, 1, 'long');
fseek(fid, 4924, 'bof');  m.h_G21_Rel3_Mr_SpatialGridShift_Lines = fread(fid, 1, 'float');
fseek(fid, 4928, 'bof');  m.h_G21_Rel3_Mr_SpatialGridShift_Columns = fread(fid, 1, 'float');
fseek(fid, 4932, 'bof');  m.h_G21_Rel3_Mr_SpatialGridShift_Partitions = fread(fid, 1, 'float');
fseek(fid, 4936, 'bof');  m.h_G21_Rel3_Mr_SignalMinimum = fread(fid, 1, 'float');
fseek(fid, 4940, 'bof');  m.h_G21_Rel3_Mr_SignalMaximum = fread(fid, 1, 'float');
fseek(fid, 4944, 'bof');  m.h_G21_Rel3_Mr_SpecInfoMask = fread(fid, 1, 'long');
fseek(fid, 4948, 'bof');  m.h_G21_Rel3_Mr_EpiTimerateChangeMagnitude = fread(fid, 1, 'float');
fseek(fid, 4952, 'bof');  m.h_G21_Rel3_Mr_EpiTimerate_x_comp = fread(fid, 1, 'float');
fseek(fid, 4956, 'bof');  m.h_G21_Rel3_Mr_EpiTimerate_y_comp = fread(fid, 1, 'float');
fseek(fid, 4960, 'bof');  m.h_G21_Rel3_Mr_EpiTimerate_z_comp = fread(fid, 1, 'float');
fseek(fid, 4964, 'bof');  m.h_G21_Rel3_Mr_EpiTimerateLegalLimit1 = fread(fid, 1, 'float');
fseek(fid, 4968, 'bof');  m.h_G21_Rel3_Mr_EpiOperationModeFlag = fread(fid, 1, 'long');
fseek(fid, 4972, 'bof');  m.h_G21_Rel3_Mr_EpiFieldCalculationSafetyFactor = fread(fid, 1, 'float');
fseek(fid, 4976, 'bof');  m.h_G21_Rel3_Mr_EpiLegalLimit1ChangeValue = fread(fid, 1, 'float');
fseek(fid, 4980, 'bof');  m.h_G21_Rel3_Mr_EpiLegalLimit2ChangeValue = fread(fid, 1, 'float');
fseek(fid, 4984, 'bof');  m.h_G21_Rel3_Mr_EpiRiseTime = fread(fid, 1, 'float');
fseek(fid, 4992, 'bof');  m.h_G28_Pre_ImageDimension = fread(fid, 1, 'short');
fseek(fid, 4994, 'bof');  m.h_G28_Pre_Rows = fread(fid, 1, 'short');
fseek(fid, 4996, 'bof');  m.h_G28_Pre_Columns = fread(fid, 1, 'short');
fseek(fid, 5000, 'bof');  m.h_G28_Pre_PixelSize_Row = fread(fid, 1, 'double');
fseek(fid, 5008, 'bof');  m.h_G28_Pre_PixelSize_Column = fread(fid, 1, 'double');
fseek(fid, 5024, 'bof');  m.h_G28_Pre_BitsAllocated = fread(fid, 1, 'short');
fseek(fid, 5026, 'bof');  m.h_G28_Pre_BitsStored = fread(fid, 1, 'short');
fseek(fid, 5028, 'bof');  m.h_G28_Pre_HighBit = fread(fid, 1, 'short');
fseek(fid, 5030, 'bof');  m.h_G28_Pre_PixelRepresentation = fread(fid, 1, 'short');
fseek(fid, 5048, 'bof');  m.h_G28_Pre_RescaleIntercept = fread(fid, 1, 'long');
fseek(fid, 5052, 'bof');  m.h_G28_Pre_RescaleSlope = fread(fid, 1, 'long');
fseek(fid, 5056, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G28_Pre_PatientName_DICOM = char(t);
fseek(fid, 5121, 'bof');  t = fread(fid, 65, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G28_Pre_PatientId_DICOM = char(t);
fseek(fid, 5188, 'bof');  m.h_G28_Pre_StudyDate_DICOM_Year = fread(fid, 1, 'long');
fseek(fid, 5192, 'bof');  m.h_G28_Pre_StudyDate_DICOM_Month = fread(fid, 1, 'long');
fseek(fid, 5196, 'bof');  m.h_G28_Pre_StudyDate_DICOM_Day = fread(fid, 1, 'long');
fseek(fid, 5200, 'bof');  m.h_G28_Pre_StudyTime_DICOM_Hour = fread(fid, 1, 'long');
fseek(fid, 5204, 'bof');  m.h_G28_Pre_StudyTime_DICOM_Minute = fread(fid, 1, 'long');
fseek(fid, 5208, 'bof');  m.h_G28_Pre_StudyTime_DICOM_Second = fread(fid, 1, 'long');
fseek(fid, 5212, 'bof');  m.h_G28_Pre_StudyTime_DICOM_Fraction = fread(fid, 1, 'long');
fseek(fid, 5284, 'bof');  m.h_G29_Pre_SortCode = fread(fid, 1, 'long');
fseek(fid, 5288, 'bof');  m.h_G29_Pre_Splash = fread(fid, 1, 'long');
fseek(fid, 5504, 'bof');  t = fread(fid, 13, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_PatientNumber = char(t);
fseek(fid, 5517, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_PatientSexAndAge = char(t);
fseek(fid, 5529, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_PatientPosition = char(t);
fseek(fid, 5541, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_ImageNumber = char(t);
fseek(fid, 5553, 'bof');  t = fread(fid, 6, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_Label = char(t);
fseek(fid, 5559, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_DateOfMeasurement = char(t);
fseek(fid, 5571, 'bof');  t = fread(fid, 6, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_TimeOfMeasurement = char(t);
fseek(fid, 5577, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_TimeOfAcquisition = char(t);
fseek(fid, 5589, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_NumberOfAcquisitions = char(t);
fseek(fid, 5601, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_CommentNo1 = char(t);
fseek(fid, 5628, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_CommentNo2 = char(t);
fseek(fid, 5655, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_InstallationName = char(t);
fseek(fid, 5682, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_SoftwareVersion = char(t);
fseek(fid, 5694, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_Matrix = char(t);
fseek(fid, 2706, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_TypeOfMeasurement = char(t);
fseek(fid, 5718, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_ScanNumber = char(t);
fseek(fid, 5730, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_RepetitionTime = char(t);
fseek(fid, 5742, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_EchoTime = char(t);
fseek(fid, 5754, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_GatingAndTrigger = char(t);
fseek(fid, 5766, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_TubeCurrent = char(t);
fseek(fid, 5778, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_TubeVoltage = char(t);
fseek(fid, 5790, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_SliceThickness = char(t);
fseek(fid, 5802, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_SlicePosition = char(t);
fseek(fid, 5814, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_SliceOrientationNo1 = char(t);
fseek(fid, 5826, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_SliceOrientationNo2 = char(t);
fseek(fid, 5838, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_FieldOfView = char(t);
fseek(fid, 5850, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_ZoomCenter = char(t);
fseek(fid, 5862, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_GantryTilt = char(t);
fseek(fid, 5874, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_TablePosition = char(t);
fseek(fid, 5886, 'bof');  t = fread(fid, 4, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_MipHeadLine = char(t);
fseek(fid, 5890, 'bof');  t = fread(fid, 16, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_MipLine = char(t);
fseek(fid, 5906, 'bof');  t = fread(fid, 16, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_MipColumn = char(t);
fseek(fid, 5922, 'bof');  t = fread(fid, 16, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_MipSlice = char(t);
fseek(fid, 5938, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_StudyNumber = char(t);
fseek(fid, 5950, 'bof');  t = fread(fid, 6, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_Contrast = char(t);
fseek(fid, 5956, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_PatientBirthdate = char(t);
fseek(fid, 5968, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_SequenceInformation = char(t);
fseek(fid, 5980, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_SaturationRegions = char(t);
fseek(fid, 5992, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_DataSetId = char(t);
fseek(fid, 6019, 'bof');  t = fread(fid, 12, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_MagnificationFactor = char(t);
fseek(fid, 6031, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_ManufacturerModel = char(t);
fseek(fid, 6058, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_PatientName = char(t);
fseek(fid, 6085, 'bof');  t = fread(fid, 9, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_TimeOfScanning = char(t);

fclose(fid);

% eof
