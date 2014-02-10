#!/usr/bin/python


import sys, getopt, warnings, os, re

import itertools, copy, math

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem




def OwnScore( allPeaks_list, current_mol, bBreakRing, precursor_type, bRankSum, dMass_Tolerance_Fragment_Ions,iFragmentation_Depth, iParentMassWindow_list, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list):
# This function manage to score a a given compound against a given scan
# allPeaks_list: list of peaks; current_mol: mol of the original given compound; bBreakRing: whether break ring bonds; precursor_type: precursor type 1 is positive, otherwise, negative; bRankSum: always false; dMass_Tolerance_Fragment_Ions: max mass error allowed for fragment mass ; iFragmentation_Depth : max depth in the depth-first search of the fragmentation tree; iParentMassWindow_list: parent mass window; iPositive_Ion_Fragment_Mass_Windows_list: fragment mass window under positive mode; iNegative_Ion_Fragment_Mass_Windows_list : fragment mass window under negative mode
# This function manage to score a a given compound against a given scan
# return value
# dCurrentScore: score of the compound; dCurrentEnergy: not used in the current version; iIdentifiedPeak: number of peaks identified; sAnnotation_list: not used in the current version ; sOtherInfo :not used in the current version; sAFT2_Info_list: list of strings for annotation file (.AFT2).
    sOtherInfo = ""  # won't get real value in the current version
    dCurrentScore = 0
    dCurrentEnergy = 0 # won't get real value in the current version
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
    if (bRankSum) : # RankSum mode won't be picked in the current version
        iUnIdentifiedPeak = len(allPeaks_list) - iIdentifiedPeak
        dCurrentScore = CalculateRankSumScore(all_compound_match, iUnIdentifiedPeak)
    return dCurrentScore, dCurrentEnergy, iIdentifiedPeak, sAnnotation_list, sOtherInfo, sAFT2_Info_list

def CalculateRankSumScore(all_compound_match, iUnIdentifiedPeak) :
    # this function calculate the ranksum score, which won't be used in the current version
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

# This function conducts a Depth-first search. Before reaching the max depth, a fragment will be tried to generate new fragments by breaking linear bond or ring bonds. Each time, two new fragments are generated, one will be processed soon, another will be put in the unprocesedkid filed. If a fragment can't generate new fragments, it will be deleted from the data structure.
#current_mol: original compound mol (mol is a data structure used in RDKit for representing a compound or fragment); iBondsNum: total number of bonds; allPeaks_list: list of peak information; iDepth: max depth in the depth-first search of the fragmentation tree; peakmatch_list: list of peak match information ; bBreakRing: whether break ring bonds; precursor_type: 1 positive precursor, -1 negative; total_fragment_list: not used in the current version; observed_fragment_list: not used in the current version; all_compound_match: not used in the current version; dMass_Tolerance_Fragment_Ions: max mass error for fragment when it is against a peak; iPositive_Ion_Fragment_Mass_Windows_list: fragment mass window under positive mode; iNegative_Ion_Fragment_Mass_Windows_list: fragment mass window under negative mode; root_observed_status: 0 original compound is obersved, 1 not; dF_root: plausibility score of the original compound 

    root_node = [Chem.EditableMol(current_mol),[],0] # editable_mol,list of list of removed bonds, depth
    current_ring_bonds_list, current_linear_bonds_list = ClassifyBonds(root_node[0].GetMol(), bBreakRing)
    current_ringbonds_iter = itertools.combinations(current_ring_bonds_list, 2)
    current_ringbonds_combination_list = list(current_ringbonds_iter)
    unprocessedKid = []
# dF_self is S(F) of the current fragment
    dF_self = dF_root
#C_vlaue is not really used in the current version
    C_value = 0
    
    #storedNodes is a list containing unprocessed fragments
    storedNodes = [[root_node, current_linear_bonds_list, current_ringbonds_combination_list, unprocessedKid, C_value, root_observed_status, dF_self]]
    while (len(storedNodes) > 0) :
        #check the last fragment, if it has unproecessed kid, deal with it, otherwise generate new kid. If even no new kids, delete this fragment
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
    #dF is plausibility score; dF_parent is the plausibility score of the fragment;parent_observed_status: 1 parent fragment is observed, 0 not;C_self is the number of cleaved bonds from the parent fragment to the current fragment
    #C_self = 0, if the current fragment is the original compound; C_self = 1 for breaking linear bond case; C_self = 2 for breaking ring bonds case

    dF       = 1.0
    
    if (parent_observed_status == 0) : # yes
        dF = math.pow(0.5, C_self) * dF_parent
    else :
        dF = math.pow(0.1, C_self) * dF_parent
    return  dF

def mypnorm (dMean, dStandardDeviation, dRandomVariable) :
# calculate the probability of a value given mean and standard deviation under normal distribution
#dMean: mean value;  dStandardDeviation: standard deviation; dRandomVariable: value of the random variable
    dZScore = ( dRandomVariable - dMean ) / dStandardDeviation
    dProbability = 0.5 * math.erfc( -1 * dZScore / math.sqrt( 2.0 ) )
    return dProbability

def SubScore(dIntensity, dErrorDa, dF, dMass_Tolerance_Fragment_Ions, current_mz_offset):
#This function calculate the subscore between a fragment against a peak, returns the subscore
#dIntensity: peak intensity; dErrorDa: mass difference between peak and fragment; dF: plausibility score; dMass_Tolerance_Fragment_Ions: max mass error allowed for the difference between a fragment and a peak; current_mz_offset: offset value;
    dSubScore = 0
    dErrorScore = ( 1.0 - mypnorm( 0, ( dMass_Tolerance_Fragment_Ions / 2.0), math.fabs( dErrorDa  ) ) ) * 2.0
    dSubScore = dIntensity * dErrorScore * dF
    return dSubScore

def MapMass(current_dMass, allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, FragmentBonds_list,  precursor_type, total_fragment_list, observed_fragment_list, current_depth, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, C_value, dF_parent, C_self, parent_observed_status) :
# claculate score in this function. This function tries to match a fragment against peaks. If a peak has more than one matches, choose the fragment with highest score
#current_dMass: mass of the current fragment; allPeaks_list: all peaks information; peakmatch_list: list of peak match information; current_sFragmentFormula: formula of the current fragment ; current_smiles: smiles of the current fragment; FragmentBonds_list: bonds broken for generating the current fragment from the original compound;  precursor_type: 1 positive precursor, -1 negative; total_fragment_list: not used in the current version; observed_fragment_list: not used in the current version; current_depth: current depth of the fragmentation in the fragmentation tree; dMass_Tolerance_Fragment_Ions: max mass error allowed for the fragment; iPositive_Ion_Fragment_Mass_Windows_list: fragment mass window under positive mode; iNegative_Ion_Fragment_Mass_Windows_list: fragment mass window under negative mode; C_value: not used in the current version; dF_parent: plausibility score of the parent fragment; C_self : the number of cleaved bonds from the parent fragment to the current fragment; parent_observed_status: 0 parent fragment is observed, 1 not
# claculate score in this function. This function tries to match a fragment against peaks. If a peak has more than one matches, choose the fragment with highest score
#return value
#bFindPeak: whether at least one peak matches this fragment; dBDE: not used in the current version; observed_status: 0 this fragment is observed, 1 not ; dF: plausibility score
    total_fragment_list[current_depth] += 1
    z_list = [1] 
    if  (precursor_type == 1) :
        mz_windows_list = iPositive_Ion_Fragment_Mass_Windows_list
    else :
        mz_windows_list = iNegative_Ion_Fragment_Mass_Windows_list
    #print mz_windows_list
    dHMass = 1.007825
    bFindPeak = False
    observed_status = 1 # 0 for observed; 1 for not
    dBDE = 0 # won't have real value in the current version
    sBondInfo = "" # won't have real value in the current version
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
# current_fragment_mol: mol of current fragment; FragmentBonds_list: bonds broken for generating the current fragment from the original compound, it is not used in this function 
    # returns  mass (current_dMass), formula (current_smiles), and smile of a fragment (current_sFragmentFormula)
    current_dMass    = Descriptors.ExactMolWt(current_fragment_mol)
    current_sFragmentFormula = AllChem.CalcMolFormula(current_fragment_mol)
    current_smiles=  Chem.MolToSmiles (current_fragment_mol)
    return current_dMass, current_sFragmentFormula, current_smiles

def ClassifyBonds(current_mol, bBreakRing) :
    # classify ring bonds and linear bonds
# current_mol: the mol of current fragment; bBreakRing: whether consider ring bonds
#ring_bonds_list: list of ring bonds; linear_bonds_list: list of linear bonds.
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
    # split the current fragment by breaking one linear bond or two ring bonds
# current_mol: current fragment mol;  bonds_list: list of bonds, which should be broken (1 linear bond or 2 ring bonds)
# return values
# current_fragments_list: new generated fragments; bValidOperation : whether fragments generated this time are valid, if generated exactly 2, yes, otherwise, invalid.
    em = Chem.EditableMol(current_mol)
    for each_removable_bond in bonds_list :
        idx_beginAtom = each_removable_bond.GetBeginAtomIdx()
        idx_endAtom   = each_removable_bond.GetEndAtomIdx()
        em.RemoveBond(idx_beginAtom, idx_endAtom)
    current_modified_mol = em.GetMol()
    current_fragments_list = Chem.GetMolFrags(current_modified_mol, asMols=True, sanitizeFrags=False)
    # After removing one linear bond, or two ring bonds, we expect to see two new fragments, and we don't accept other cases 
# e.g., if we have more than two new fragments, or only one fragment (it could happen after two ring bonds are removed), we don't accept these bad cases. 
    if len(current_fragments_list) == 2 :
        bValidOperation = True
    else :
        bValidOperation = False
    #print len( Chem.GetMolFrags(current_mol, asMols=True)  )
    return current_fragments_list, bValidOperation


def processKid(current_editable_mol, current_removebond_list, iCurrent_depth, allPeaks_list, peakmatch_list,  bBreakRing, precursor_type, total_fragment_list, observed_fragment_list, all_compound_match, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, C_self, C_parent, parent_observed_status, dF_parent) :
    # This function score a fragment by calling MapMass. Also this function preprocess a fragment for generating kids
    #C_value is not used in the current version
#current_editable_mol: the fragment mol; current_removebond_list: list of bonds broken for generating the current fragment; iCurrent_depth: depth of the current fragment in the fragmentation tree; allPeaks_list: list of peaks in the scan; peakmatch_list: peak match info, store subscore information;   bBreakRing: whether break ring bonds; precursor_type: 1 positive precursor, otherwise negative; total_fragment_list: not used in the current version; observed_fragment_list: not used in the current version; all_compound_match: not used in the current version; dMass_Tolerance_Fragment_Ions: max mass error allowed for a fragment; iPositive_Ion_Fragment_Mass_Windows_list: fragment mass window under positive mode iNegative_Ion_Fragment_Mass_Windows_list: fragment mass window under negative mode; C_self: the number of cleaved bonds from the parent fragment to the current fragment; C_parent : C_self of the parent fragment;  parent_observed_status: 0 parent fragment observed, 1 not; dF_parent: the plausibility score of parent fragment
#This function turns a processed fragment
# current_editable_mol: the fragment mol; current_removebond_list: bonds broken for generating the current fragment; iCurrent_depth : the depth of the current fragment in the fragmentation tree; current_linear_bonds_list: linear bonds of the current fragment; current_ringbonds_combination_list: put 2 ring bonds together as a combination, the list of these combinations; unprocessedKid: unprocessed new fragment; C_value: not used in the current version; observed_status: 0, the current fragment is observed, 1 not; dF_self: the plausibility score of the current fragment.
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
# This function manages to calculate subscores of each peak. subscore information are stored in peakmatch_list
# current_mol: mol of the compound; allPeaks_list: list of peaks ; peakmatch_list: store subscore and related information of each peak; bBreakRing : whether break ring bonds; precursor_type: 1 positive precursor, otherwise negative; dMass_Tolerance_Fragment_Ions: max mass error allowed for a fragment; iFragmentation_Depth: max depth in the depth-first search of the fragmentation tree; iPositive_Ion_Fragment_Mass_Windows_list: fragment mass window under positive mode; iNegative_Ion_Fragment_Mass_Windows_list: fragment mass window under negative mode.
# This function manages to calculate subscores of each peak. subscore information are stored in peakmatch_list
# no return values will not accutally used in the current version
    total_fragment_list = [0 for i in range(4)]
    observed_fragment_list = [0 for i in range(4) ]
    all_compound_match = []
    
    current_sFragmentFormula = AllChem.CalcMolFormula(current_mol)
    current_smiles = Chem.MolToSmiles(current_mol)
    # dF_self S(F) of current fragment
    dF_self = 1.0
    root_observed_status = 0 # initialization, 0 observed, 1 not
    bRootFragment = True
    #calculate score of the original compound (root fragment)
    bFindPeak, dCurrentFragmentBDE, root_observed_status, dF_self = MapMass(Descriptors.ExactMolWt(current_mol), allPeaks_list, peakmatch_list, current_sFragmentFormula, current_smiles, [],   precursor_type, total_fragment_list,  observed_fragment_list, 0, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, 0, 1, 0, 0)
    all_compound_match.append([bFindPeak, dCurrentFragmentBDE])
    iBondsNum = current_mol.GetNumBonds()
    TreeLikeBreakBondsDepthFirst(current_mol, iBondsNum, allPeaks_list, iFragmentation_Depth, peakmatch_list,   bBreakRing, precursor_type, total_fragment_list, observed_fragment_list, all_compound_match, dMass_Tolerance_Fragment_Ions, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, root_observed_status, dF_self)
    return total_fragment_list, observed_fragment_list, all_compound_match
