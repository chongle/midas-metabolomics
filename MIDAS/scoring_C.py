#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

import weightedscore


def BinarySearch_Upper(target_list, bounder_value):
    # return the upper bound of the index range
    ele_num = len(target_list)
    upper_index = ele_num-1
    if (target_list[upper_index] < bounder_value):
        return upper_index
    lower_index = 0
    if (target_list[lower_index] > bounder_value) :
        return -1
    while ((upper_index - lower_index)>1):
        middle_index = (lower_index + upper_index)/2
        if (target_list[middle_index]> bounder_value ) :
            upper_index = middle_index
        else :
            lower_index = middle_index
    # upper_index is the index of element with value just higher than the upper bound
    # lower_index is the index of element with value just lower than the upper bound
    return lower_index
    
def BinarySearch_Lower(target_list, bounder_value):
    ele_num = len(target_list)
    upper_index = ele_num-1
    if (target_list[upper_index] < bounder_value):
        return -1
    lower_index = 0
    if (target_list[lower_index] > bounder_value) :
        return lower_index
    while ((upper_index - lower_index)>1):
        middle_index = (lower_index + upper_index)/2
        if (target_list[middle_index]>= bounder_value ) : # !!!!
            upper_index = middle_index
        else :
            lower_index = middle_index
    #upper_index is the index of element with value just higher than the lower bound
    #lower_index is the index of element with value just lower than the lower bound
    return upper_index

def GetRelatedCompound(Compound_list, precursor_mz, precursor_accuracy):
    compound_mass_list = [each_compound[2] for each_compound in Compound_list]
#    print compound_mass_list
    upper_compound_mass = precursor_mz + precursor_accuracy
    lower_compound_mass = precursor_mz - precursor_accuracy
#    print upper_compound_mass, lower_compound_mass
    upper_compound_index=BinarySearch_Upper(compound_mass_list, upper_compound_mass)
    #print upper_compound_index, compound_mass_list[upper_compound_index], upper_compound_mass, compound_mass_list[upper_compound_index+1]
    lower_compound_index=BinarySearch_Lower(compound_mass_list, lower_compound_mass)
    #print lower_compound_index, compound_mass_list[lower_compound_index-1], lower_compound_mass, compound_mass_list[lower_compound_index] 
    QueryCompound_list = []
    if ((lower_compound_index != -1) and (upper_compound_index != -1)) :
        QueryCompound_list = Compound_list[lower_compound_index : upper_compound_index+1]
    return QueryCompound_list


def NormalizeIntensity(allPeaks_list) :
    max_intensity = max([each_peak[1] for each_peak in allPeaks_list])
    for i in range(len(allPeaks_list)) :
        allPeaks_list[i][1] = allPeaks_list[i][1]/max_intensity
#        allPeaks_list[i][1] = 1
    return allPeaks_list


def score_main(Compound_list, sOutput_Filename, bSpectrumDetails, bBreakRing, dCurrentPrecursor_type, bRankSum, dCurrentParentMass, current_peaks_list, sCurrentScanNumber, mylock, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, sFT2_basename, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, dPrecursorMofZ, sRetentionTime, sAnnotation_Filename):
    #print sCurrentScanNumber
    precursor_accuracy = dMass_Tolerance_Parent_Ion
    Compound_Scores_list = []
    allPeaks_list = NormalizeIntensity(current_peaks_list)
    QueryCompound_list = GetRelatedCompound(Compound_list, dCurrentParentMass, precursor_accuracy)
    # Compounds satisfied with precursor mass range are saved in QueryCompound_list 
    #print len(QueryCompound_list)
    for each_compound in QueryCompound_list :
        current_mol = Chem.MolFromInchi(each_compound[1])
        current_fragments_list = Chem.GetMolFrags(current_mol, asMols=True, sanitizeFrags=False)
        if (len(current_fragments_list) != 1) :
            continue
        # get scores. dCurrentEnergy and sOtherInfo don't have valid value in the current version
        dCurrentScore, dCurrentEnergy, iIdentifiedPeak, sAnnotation_list, sOtherInfo, sAFT2_Info_list= weightedscore.OwnScore( allPeaks_list, current_mol, bBreakRing, dCurrentPrecursor_type, bRankSum, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, iParentMassWindow_list, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list)
        if (dCurrentPrecursor_type == 1) :
            sScanType = "PositiveIon"
        else :
            sScanType = "NegativeIon"
        Compound_Scores_list.append([dCurrentScore, dCurrentEnergy, each_compound[0], each_compound[1], each_compound[3], iIdentifiedPeak, sOtherInfo, sAnnotation_list, sScanType, each_compound[2], each_compound[4], sAFT2_Info_list])
    iPeakNum = len(allPeaks_list)
    Compound_Scores_list.sort(key=lambda e:e[0], reverse=True)
    #print "aa"
    output_str = ""
    annotation_str = ""
    dPreviousScore = 1000 # initial score is bigger than any real scores
    dCurrentRank   = 1
    for i in range(len(Compound_Scores_list)) :
        # here we output top 5 candidates
        if (i>=5):
            break
    #    print sCurrentScanNumber, i 
        each_compound_info = Compound_Scores_list[i]
        sScanType  = each_compound_info[8]
        dParentMassError = each_compound_info[9] - dCurrentParentMass
        sIdentifier = each_compound_info[2]
        sName = each_compound_info[10]
        sOriginalInChI = each_compound_info[3]
        sIdentifiedInChI = sOriginalInChI
        sModifications = "NA"
        sLinks = each_compound_info[4]
        iIdentifiedPeak = each_compound_info[5]
        dCurrentScore = each_compound_info[0]

        if (dCurrentScore < dPreviousScore) or (i==0):
            dCurrentRank  = i+1
            dPreviousScore = dCurrentScore

        output_str += sFT2_basename + "\t" + sCurrentScanNumber + "\t"+str(dPrecursorMofZ)+"\t"+sRetentionTime+"\t" + sScanType + "\t" + str(dCurrentRank) + "\t" + str(dParentMassError) + "\t"
        annotation_str += "M\t" + sFT2_basename + "\t" + sCurrentScanNumber + "\t"+str(dPrecursorMofZ)+"\t"+sRetentionTime+"\t" + sScanType + "\t" + str(dCurrentRank) + "\t" + str(dParentMassError) + "\t"

        output_str += str(iIdentifiedPeak)+">"+str(iPeakNum)+"\t"+str(dCurrentScore)+"\t"+sIdentifier+"\t"+sName+"\t"+sOriginalInChI+"\t"
        annotation_str += str(iIdentifiedPeak)+">"+str(iPeakNum)+"\t"+str(dCurrentScore)+"\t"+sIdentifier+"\t"+sName+"\t"+sOriginalInChI+"\t"
        
        output_str += sIdentifiedInChI+"\t"+sModifications+"\t"+sLinks+"\n"
        annotation_str += sIdentifiedInChI+"\t"+sModifications+"\t"+sLinks+"\n"
        current_sAFT2_Info_list = Compound_Scores_list[i][11]
        for each_AFT2_line in current_sAFT2_Info_list :
            annotation_str += each_AFT2_line+ "\n"
   #     print dCurrentScore

    mylock.acquire()
    #sequential part to avoid writing conflict
    #meb file
    Output_File = open(sOutput_Filename, "a")
    Output_File.write(output_str)
    Output_File.close()
    #AFT2 file
    Annotation_File = open(sAnnotation_Filename, "a")
    Annotation_File.write(annotation_str)
    Annotation_File.close()
    #print sCurrentScanNumber
    mylock.release()
    #end of sequential part





