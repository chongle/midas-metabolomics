#!/usr/bin/python


import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem




def OwnScore( allPeaks_list, current_mol, bBreakRing, precursor_type, bRankSum, dMass_Tolerance_Fragment_Ions,iFragmentation_Depth, iParentMassWindow_list, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list):
#    global dMass_Tolerance_Fragment_Ions
#    dMass_Tolerance_Fragment_Ions = 0.01
#    print iParentMassWindow_list, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list
    sOtherInfo = ""
    dCurrentScore = 0
    dCurrentEnergy = 0
    iIdentifiedPeak= 0
    peakmatch_list = [[] for each_peak in allPeaks_list]
    total_fragment_list, observed_fragment_list, all_compound_match  = ExhaustBonds(current_mol, allPeaks_list, peakmatch_list, bBreakRing,precursor_type, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list)
   # print "Exhaust done!"
    dInnerProduct_Intensity = 0
    dSum_Intensity = 0
    for each_peak in allPeaks_list :
        dInnerProduct_Intensity += each_peak[1] * each_peak[1]
        dSum_Intensity += each_peak[1]
    sAnnotation_list = []
    sAFT2_Info_list  = []
    for i in range(len(allPeaks_list)) :
        sAFT2_Info  = ""
        for each_measure_info in allPeaks_list[i] :
            sAFT2_Info += str(each_measure_info) + "\t"
        for miss_info_id in range(6-len(allPeaks_list[i])):
            sAFT2_Info += "0\t"
        if peakmatch_list[i] :
            iIdentifiedPeak += 1 
            dCurrentEnergy  += peakmatch_list[i][0][5]
            dErrorDa         = peakmatch_list[i][0][4]
            dSubScore        = peakmatch_list[i][0][6]
            dCurrentScore   += dSubScore
            sAFT2_Info += peakmatch_list[i][0][2] + "\t" + peakmatch_list[i][0][0] + "\t" + str(dErrorDa) + "\t" + str(peakmatch_list[i][0][8])
        else :
            sAFT2_Info += "NA\tNA\tNA\tNA"
        sAFT2_Info_list.append(sAFT2_Info)
    if (bRankSum) :
        iUnIdentifiedPeak = len(allPeaks_list) - iIdentifiedPeak
        dCurrentScore = CalculateRankSumScore(all_compound_match, iUnIdentifiedPeak)
    return dCurrentScore, dCurrentEnergy, iIdentifiedPeak, sAnnotation_list, sOtherInfo, sAFT2_Info_list

def CalculateRankSumScore(all_compound_match, iUnIdentifiedPeak) :
    iU = 0
    iUnobservedCompound = 0
    all_compound_match.sort(key=lambda e:e[1])
    for each_compound in all_compound_match :
        bObserved = each_compound[0]
        if (bObserved) :
            iU += iUnobservedCompound
        else :
            iUnobservedCompound += 1
    iU += iUnobservedCompound * iUnIdentifiedPeak
    n_1 = len(all_compound_match) - iUnobservedCompound + iUnIdentifiedPeak
    n_2 = iUnobservedCompound
    RankSumMean = n_1 * n_2 /2.0
    RankSumD    = math.sqrt(n_1*n_2*(n_1+n_2+1)/12.0)
    if (RankSumD == 0) :
        RankSumD = 1
    dRankSumScore = math.fabs((iU - RankSumMean)/RankSumD)
    return (-1)*dRankSumScore

def TreeLikeBreakBondsDepthFirst(current_mol, iBondsNum, allPeaks_list, iDepth, peakmatch_list,  bBreakRing, precursor_type, total_fragment_list, observed_fragment_list, all_compound_match, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, root_observed_status, dF_root) :
    root_node = [Chem.EditableMol(current_mol),[],0] # editable_mol,list of list of removed bonds, depth
    current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(root_node[0].GetMol(), bBreakRing)
    current_ringbonds_iter = itertools.combinations(current_ring_bonds_list, 2)
    current_ringbonds_combination_list = list(current_ringbonds_iter)
    unprocessedKid = []
    dF_self = dF_root
    #dF_self = 1.0
    C_value = 0
    #C_self  = 0
    #root_observed_status = 1 # 1 no 0 yes
    # print "ee"
    storedNodes = [[root_node, current_linear_bonds_list, current_ringbonds_combination_list, unprocessedKid, C_value, root_observed_status, dF_self]]
    while (len(storedNodes) > 0) :
        if (len(storedNodes[-1][3]) > 0) : # unprocessed kid
            new_item = processKid(storedNodes[-1][3][0][0], storedNodes[-1][3][0][1], storedNodes[-1][0][2]+1, allPeaks_list, peakmatch_list,  bBreakRing, precursor_type, total_fragment_list, observed_fragment_list, all_compound_match, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, storedNodes[-1][3][0][2], storedNodes[-1][4], storedNodes[-1][5], storedNodes[-1][6]) 
            del storedNodes[-1][3][0]
            if (new_item[0][2] < iDepth) :
                storedNodes.append(new_item)
        elif (len(storedNodes[-1][1]) > 0) : # linear bond
            remove_bond = storedNodes[-1][1][0]
            current_fragments_list, bValidOperation = RemoveBonds(storedNodes[-1][0][0].GetMol(), [remove_bond]) 
            if (bValidOperation) :
                FragmentBonds_list  =  list( storedNodes[-1][0][1] )
                FragmentBonds_list.append(remove_bond)
                C_self = 1
                for i in range(2) :
                    storedNodes[-1][3].append([Chem.EditableMol(current_fragments_list[i]), FragmentBonds_list, C_self])
            else :
                print "wrong!"
                sys.exit(1)
            del storedNodes[-1][1][0]
        elif (len(storedNodes[-1][2]) > 0) : # ring bonds
            remove_first_bond  = storedNodes[-1][2][0][0]
            remove_second_bond = storedNodes[-1][2][0][1]
            current_fragments_list, bValidOperation = RemoveBonds(storedNodes[-1][0][0].GetMol(), [remove_first_bond, remove_second_bond])
            if (bValidOperation) :
                FragmentBonds_list  =  list( storedNodes[-1][0][1] )
                FragmentBonds_list.append(remove_first_bond)
                FragmentBonds_list.append(remove_second_bond)
                C_self = 2
                for i in range(2) :
                    storedNodes[-1][3].append([Chem.EditableMol(current_fragments_list[i]), FragmentBonds_list, C_self])
            del storedNodes[-1][2][0]
        else :
            del storedNodes[-1]

def ExactBondsInfo( dF_parent, C_self, parent_observed_status) :
    dF       = 1.0
    
    if (parent_observed_status == 0) : # yes
        dF = math.pow(0.5, C_self) * dF_parent
    else :
        dF = math.pow(0.1, C_self) * dF_parent
    return  dF

def mypnorm (dMean, dStandardDeviation, dRandomVariable) :
    dZScore = ( dRandomVariable - dMean ) / dStandardDeviation
    dProbability = 0.5 * math.erfc( -1 * dZScore / math.sqrt( 2.0 ) )
    return dProbability

def SubScore(dIntensity, dErrorDa, dF, dMass_Tolerance_Fragment_Ions, current_mz_offset):
    dSubScore = 0
    dErrorScore = ( 1.0 - mypnorm( 0, ( dMass_Tolerance_Fragment_Ions / 2.0), math.fabs( dErrorDa  ) ) ) * 2.0
    dSubScore = dIntensity * dErrorScore * dF
    return dSubScore

def MapMass(current_dMass, allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, FragmentBonds_list,  precursor_type, total_fragment_list, observed_fragment_list, current_depth, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, C_value, dF_parent, C_self, parent_observed_status) :
#    print iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list
    total_fragment_list[current_depth] += 1
    z_list = [1] 
    if  (precursor_type == 1) :
        mz_windows_list = iPositive_Ion_Fragment_Mass_Windows_list
    else :
        mz_windows_list = iNegative_Ion_Fragment_Mass_Windows_list
    #print mz_windows_list
    dHMass = 1.007825
    bFindPeak = False
    observed_status = 1
    dBDE = 0
    sBondInfo = ""
    dF = ExactBondsInfo( dF_parent, C_self, parent_observed_status)
    for i in range(len(allPeaks_list)) :
        each_peak = allPeaks_list[i]
        dMeasuredMZ = each_peak[0]
        for current_z in z_list :
            for current_mz_offset in mz_windows_list :
                mzdiff =math.fabs ((current_dMass + current_mz_offset*dHMass)/current_z  - dMeasuredMZ)
                if mzdiff <= dMass_Tolerance_Fragment_Ions :
                    observed_status = 0
                    dErrorDa = dMeasuredMZ - ((current_dMass+current_mz_offset*dHMass)/current_z)
                    dSubScore= SubScore(allPeaks_list[i][1], dErrorDa, dF, dMass_Tolerance_Fragment_Ions, current_mz_offset)
                    current_newmatch = [str(current_mz_offset),str(current_sFragmentFormula),str(current_smiles),sBondInfo,dErrorDa,dBDE, dSubScore, current_mz_offset, current_depth]
                    observed_fragment_list[current_depth] += 1
                    if (not peakmatch_list[i]) :
                        peakmatch_list[i].append(current_newmatch)
                    else :
                        old_dSubScore = peakmatch_list[i][0][6]
                        if (dSubScore > old_dSubScore) :
                           peakmatch_list[i][0]=current_newmatch
                    bFindPeak = True
                    break
    return bFindPeak, dBDE, observed_status, dF


def DumpOneFragment(current_fragment_mol, FragmentBonds_list) :
    current_dMass    = Descriptors.ExactMolWt(current_fragment_mol)
    current_sFragmentFormula = AllChem.CalcMolFormula(current_fragment_mol)
    sBondsTypes = "{"
    for current_bond in FragmentBonds_list :
        if (sBondsTypes != "{") :
            sBondsTypes += ","
        sBondsTypes += str(current_bond.GetBondType())
    sBondsTypes += "}"
    #current_SanitizedMol = Chem.SanitizeMol(current_fragment_mol)
    #current_inchi =  Chem.MolToInchi(current_fragment_mol)
    current_smiles=  Chem.MolToSmiles (current_fragment_mol)
    #print current_dMass, current_sFragmentFormula, sBondsTypes, current_smiles
    return current_dMass, current_sFragmentFormula, current_smiles

def ClassifyBonds(current_mol, bBreakRing) :
    ring_bonds_list   = []
    linear_bonds_list = []
    iBondsNum = current_mol.GetNumBonds() 
    for i in range(iBondsNum) :
        current_bond = current_mol.GetBondWithIdx(i)
        if current_bond.IsInRing()  :
            if (bBreakRing) :
                ring_bonds_list.append(current_bond)
        else :
            linear_bonds_list.append(current_bond)
    return ring_bonds_list, linear_bonds_list

def RemoveBonds(current_mol, bonds_list) :
   # print len( Chem.GetMolFrags(current_mol, asMols=True)  )
    em = Chem.EditableMol(current_mol)
    for each_removable_bond in bonds_list :
        idx_beginAtom = each_removable_bond.GetBeginAtomIdx()
        idx_endAtom   = each_removable_bond.GetEndAtomIdx()
        em.RemoveBond(idx_beginAtom, idx_endAtom)
    current_modified_mol = em.GetMol()
    current_fragments_list = Chem.GetMolFrags(current_modified_mol, asMols=True, sanitizeFrags=False)
    #print len(current_fragments_list)
    if len(current_fragments_list) == 2 :
        bValidOperation = True
    else :
        bValidOperation = False
    #print len( Chem.GetMolFrags(current_mol, asMols=True)  )
    return current_fragments_list, bValidOperation


def processKid(current_editable_mol, current_removebond_list, iCurrent_depth, allPeaks_list, peakmatch_list,   bBreakRing, precursor_type, total_fragment_list, observed_fragment_list, all_compound_match, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, C_self, C_parent, parent_observed_status, dF_parent) :
    #print current_editable_mol
    C_value = C_self * parent_observed_status + C_parent
    current_mol = current_editable_mol.GetMol()
    current_dMass, current_sFragmentFormula, current_smiles=DumpOneFragment(current_mol, current_removebond_list)
    bFindPeak, dCurrentFragmentBDE, observed_status, dF_self  = MapMass(current_dMass, allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, current_removebond_list,  precursor_type, total_fragment_list, observed_fragment_list, iCurrent_depth, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, C_value, dF_parent, C_self, parent_observed_status)
    all_compound_match.append([bFindPeak, dCurrentFragmentBDE])
    current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(current_mol, bBreakRing)
    current_ringbonds_iter = itertools.combinations(current_ring_bonds_list, 2)
    current_ringbonds_combination_list = list(current_ringbonds_iter)
    unprocessedKid = []

    new_item = [[current_editable_mol, current_removebond_list, iCurrent_depth], current_linear_bonds_list, current_ringbonds_combination_list, unprocessedKid, C_value, observed_status, dF_self]
    return new_item


def ExhaustBonds(current_mol, allPeaks_list, peakmatch_list, bBreakRing,precursor_type, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list) :
    total_fragment_list = [0 for i in range(4)]
    observed_fragment_list = [0 for i in range(4) ]
    all_compound_match = []
    
    current_sFragmentFormula = AllChem.CalcMolFormula(current_mol)
    current_smiles = Chem.MolToSmiles(current_mol)
#    print Descriptors.ExactMolWt(current_mol), "NULL", AllChem.CalcMolFormula(current_mol)
    dF_self = 1.0
    root_observed_status = 0
    bRootFragment = True
    bFindPeak, dCurrentFragmentBDE, root_observed_status, dF_self = MapMass(Descriptors.ExactMolWt(current_mol), allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, [],   precursor_type, total_fragment_list,  observed_fragment_list, 0, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, 0, 1, 0, 0)
    all_compound_match.append([bFindPeak, dCurrentFragmentBDE])
    iBondsNum = current_mol.GetNumBonds()
    TreeLikeBreakBondsDepthFirst(current_mol, iBondsNum, allPeaks_list, iFragmentation_Depth, peakmatch_list,   bBreakRing, precursor_type, total_fragment_list, observed_fragment_list, all_compound_match, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, root_observed_status, dF_self)
    return total_fragment_list, observed_fragment_list, all_compound_match
