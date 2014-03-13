#!/usr/bin/python

import sys, getopt, warnings, os, re

import itertools, copy, math
import multiprocessing
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

import scoring_C
import parseconfig

def parse_options (argv) :
    opts, args = getopt.getopt(argv[1:], "hf:c:o:a:",
                                ["help", "ft2-filename", "configure-filename", "output-result-filename", "output-annotation-filename"])
    sFT2_Filename      = ""
    sOutput_Filename    = ""
    sConfigure_Filename = ""
    sAnnotation_Filename = ""

    for option, value in opts:
        if option in ("-h", "--help"):
            print "-f ft2-filename -c configure-filename -o output-result-filename -a output-annotation-filename"
            sys.exit(0)
        if option in ("-f", "--ft2-filename"):
            sFT2_Filename = value
        if option in ("-c", "--configure-filename"):
            sConfigure_Filename = value
        if option in ("-o", "--output-filename"):
            sOutput_Filename = value
        if option in ("-a", "--output-annotation-filename"):
            sAnnotation_Filename = value
    
    if (sFT2_Filename == "") :
        print "Please specify FT2 file name"
        sys.exit(1)
    if (sOutput_Filename == "") :
        sOutput_Filename =  os.path.splitext(sFT2_Filename)[0]+".meb"
    if (sAnnotation_Filename == "") :
        sAnnotation_Filename = os.path.splitext(sFT2_Filename)[0]+".AFT2"
    
    if (sConfigure_Filename == "") :
        print "Please specify configure file name"
        sys.exit(1)

    return [sFT2_Filename,  sOutput_Filename, sConfigure_Filename, sAnnotation_Filename]


## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

       # if len(file_list) == 0:
        #    print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            # die("Program exit!")
	 #   sys.exit(0)
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        sys.exit(1)

    return file_list


def HandleOneRealHit(compound_list, sOutput_Filename, bSpectrumDetails, bBreakRing, dCurrentPrecursor_type, bRankSum, dCurrentParentMass, current_peaks_list, sCurrentScanNumber, mylock, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, sFT2_basename, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, dPrecursorMofZ, sRetentionTime, sAnnotation_Filename) :
 
#Call scoring_C.score_main
#compound_list :list of compounds of database; sOutput_Filename : outputfile (.meb) file name;  bSpectrumDetails: not used in the current version; bBreakRing : whether break ring bonds; dCurrentPrecursor_type : Precusor type, 1 indicates positive, otherwise negative;  bRankSum: alway false, not used in the current version; dCurrentParentMass: precursor mass; current_peaks_list: peak information; sCurrentScanNumber: scan id; mylock: lock for synchronization; iParentMassWindow_list: parent mass window; dMass_Tolerance_Parent_Ion : parent mass error range;  dMass_Tolerance_Fragment_Ions: peak mass error range; iFragmentation_Depth : search depth in the fragmentation tree;  sFT2_basename: FT2 file base name; iPositive_Ion_Fragment_Mass_Windows_list: fragment mass window if positive; iNegative_Ion_Fragment_Mass_Windows_list: fragment mass window if negative; dPrecursorMofZ: precursor M/Z;  sRetentionTime: retention time;  sAnnotation_Filename: annotation file (.AFT2) file name.
 

    scoring_C.score_main(compound_list, sOutput_Filename, bSpectrumDetails, bBreakRing, dCurrentPrecursor_type, bRankSum, dCurrentParentMass, current_peaks_list, sCurrentScanNumber, mylock, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, sFT2_basename, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, dPrecursorMofZ, sRetentionTime, sAnnotation_Filename)



def CalCompoundMass(sInchi) :
#sInchi : string of Inchi
# This function calculates the mass of the compound according to its Inchi data, returns the mass
    current_mol  = Chem.MolFromInchi(sInchi)
    current_mass = Descriptors.ExactMolWt(current_mol)
    return current_mass

def ReadCompoundFile(compound_filename) :
# Read compound file 
#compound_filename : name of the database file
# return all information of all compounds in the database
    #print "Start Reading Database..."
    compound_list  = [] # each entry: ID, Inchi, Mass, DBLink, Name
    compound_file = open(compound_filename)
    bTitleLine = True
    for each_line in compound_file :
        if (bTitleLine) :
            bTitleLine = False
            continue
        each_line = each_line.strip()
    #    print each_line
        if ((each_line == "") or (each_line.startswith("#"))) :
            continue
        current_compound_info = each_line.split("\t")
    # format check
        if (len(current_compound_info) != 4) :
            print "illegal compound", each_line
            sys.exit(1)
        else :
            dMass     = CalCompoundMass(current_compound_info[2])
            #compound_list.append([current_compound_info[0], updated_inchi, dMass, current_compound_info[3], current_compound_info[1]])
            compound_list.append([current_compound_info[0], current_compound_info[2], dMass, current_compound_info[3], current_compound_info[1]])

    compound_file.close()
    compound_list.sort(key=lambda e:e[2]) # sort based on mass
    #print "Finish reading Database."
    return compound_list


def HandleAllRealHits(all_scans_list, compound_filename, sOutput_Filename, process_number, bSpectrumDetails, bBreakRing, bRankSum, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, sFT2_basename, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, sAnnotation_Filename) :
   
    #This function has the multi-processing part

#all_scans_list: list of all scans; compound_filename : database file name; sOutput_Filename: output file (.meb) name;  process_number: the max number of processes allowed to parallel handle compounds; bSpectrumDetails: not used in the current version; bBreakRing: whether break ring bonds; bRankSum: not used in the current version, always false; iParentMassWindow_list: parent mass window; dMass_Tolerance_Parent_Ion: max mass error allowed for precursor mass; dMass_Tolerance_Fragment_Ions: max mass error allowed for peaks; iFragmentation_Depth: max depth in the depth-first search of the fragmentation tree; sFT2_basename: FT2 file name; iPositive_Ion_Fragment_Mass_Windows_list: fragmenet mass window under positive mode; iNegative_Ion_Fragment_Mass_Windows_list: fragment mass window under negative mode; sAnnotation_Filename: annotation file (.AFT2) name.

    #dProtonMass = 1.007825 # proton mass
    mypool = Pool(processes=process_number)
    compound_list = ReadCompoundFile(compound_filename)
    manager = multiprocessing.Manager() 
    mylock = manager.Lock() 
    for each_scan in all_scans_list :
        current_peaks_list     = each_scan[0]
        dCurrentParentMass     = each_scan[1]
        sCurrentScanNumber     = each_scan[2]
        dCurrentPrecursor_type = each_scan[3]
        dPrecursorMofZ         = each_scan[4]
        sRetentionTime         = each_scan[5]
        #print dCurrentParenMass, sCurrentScanNumber, dCurrentPrecursor_type
        result = mypool.apply_async(HandleOneRealHit,(compound_list, sOutput_Filename, bSpectrumDetails, bBreakRing, dCurrentPrecursor_type, bRankSum, dCurrentParentMass, current_peaks_list, sCurrentScanNumber, mylock, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, sFT2_basename, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, dPrecursorMofZ, sRetentionTime, sAnnotation_Filename))
    mypool.close()
    mypool.join()


def ReadWindowInfo(sWindow) :
#This function parse the string of window setting and return it
# sWindow is the string containing the window setting
    sWindow = sWindow.strip(",")
    sWindow_list = sWindow.split(",")
    iWindow_list = []
    for each_mass in sWindow_list :
        iWindow_list.append(int(each_mass))

    return iWindow_list

def ReadConfigureFile(sConfigure_Filename) :
#Read configure file
#sConfigure_Filename : configure file name
# return value 
# compound_filename : path of database file; iDefault_Polarity: default polarity, 1 positive, -1  negative; iDefault_Charge_State: default charge state; iParentMassWindow_list : parent mass window; dMass_Tolerance_Parent_Ion: max error allowed for the precursor mass; dMass_Tolerance_Fragment_Ions: max error allowed for the peaks; bBreakRing : whether break ring bonds; iFragmentation_Depth: max depth of depth-first search will reach in the fragmentation tree; process_number: number of processes will created for parallel handling compounds; bRankSum: not used in this version; iPositive_Ion_Fragment_Mass_Windows_list: peak mass window for positive mode; iNegative_Ion_Fragment_Mass_Windows_list: peak mass window for negative mode.
    wholeDict = parseconfig.parseConfigKeyValues(sConfigure_Filename)
    compound_filename = wholeDict.get("[Metabolite_Identification]Metabolite_Database")
    #sIonization_Mode  = wholeDict.get("[Metabolite_Identification]Ionization_Mode")
    sDefault_Polarity = wholeDict.get("[Metabolite_Identification]Default_Polarity")
    iDefault_Charge_State = int(wholeDict.get("[Metabolite_Identification]Default_Charge_State"))
    sParent_Mass_Windows = wholeDict.get("[Metabolite_Identification]Parent_Mass_Windows")
    sPositive_Ion_Fragment_Mass_Windows = wholeDict.get("[Metabolite_Identification]Positive_Ion_Fragment_Mass_Windows")
    sNegative_Ion_Fragment_Mass_Windows = wholeDict.get("[Metabolite_Identification]Negative_Ion_Fragment_Mass_Windows")
    dMass_Tolerance_Parent_Ion = float( wholeDict.get("[Metabolite_Identification]Mass_Tolerance_Parent_Ion") )
    dMass_Tolerance_Fragment_Ions = float(wholeDict.get("[Metabolite_Identification]Mass_Tolerance_Fragment_Ions"))
    sBreak_rings = wholeDict.get("[Metabolite_Identification]Break_rings")
    iFragmentation_Depth = int(wholeDict.get("[Metabolite_Identification]Fragmentation_Depth"))
    process_number = int(wholeDict.get("[Metabolite_Identification]Number_of_Processes"))
    bRankSum = False

    if (sDefault_Polarity == "positive") :
        iDefault_Polarity = 1
    else:
        iDefault_Polarity = -1

    iParentMassWindow_list = ReadWindowInfo(sParent_Mass_Windows)
    iPositive_Ion_Fragment_Mass_Windows_list = ReadWindowInfo(sPositive_Ion_Fragment_Mass_Windows)
    iNegative_Ion_Fragment_Mass_Windows_list = ReadWindowInfo(sNegative_Ion_Fragment_Mass_Windows)

    #print iParentMassWindow_list, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list 

    if ((sBreak_rings == "true") or (sBreak_rings == "True") or (sBreak_rings == "TRUE")):
        bBreakRing = True
    else :
        bBreakRing = False
    
    return [compound_filename, iDefault_Polarity, iDefault_Charge_State, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, bBreakRing, iFragmentation_Depth, process_number, bRankSum, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list]


def OutputFileHead(sConfigure_Filename) :
# organize the string for the beginning part of the output file, return this string
# sConfigure_Filename: file name of configure file
    sOutputhead = ""
    ConfigureFile = open(sConfigure_Filename)
    for each_line in ConfigureFile :
        each_line = each_line.strip()
        sOutputhead += "#\t" + each_line + "\n"
    
    sOutputhead += "#\n"
    sOutputhead += "#\tFilename: Filename of the input MS2 file\n"
    sOutputhead += "#\tScanNumber: Scan number(s). Merged scans will have multiple numbers\n"
    sOutputhead += "#\tPrecursorMZ: precursor m/z\n"
    sOutputhead += "#\tRetentionTime: retention time\n"
    sOutputhead += "#\tScanType: Scan types, e.g. positive/negative ion mode\n"
    sOutputhead += "#\tRank: Rank of this identification\n"
    sOutputhead += "#\tParentMassError: mass error of the parent ion\n"
    sOutputhead += "#\tExplainedPeaks: e.g. 8>10 means, out of the 10 peaks, 8 are explained by this compound\n"
    sOutputhead += "#\tScore: Score of this identification\n"
    sOutputhead += "#\tIdentifier: Identifier from the database\n"
    sOutputhead += "#\tName: Name from the database\n"
    sOutputhead += "#\tOriginalInChI: Original InChI from the database\n"
    sOutputhead += "#\tIdentifiedInChI: Identified InChI with potential modifications\n"
    sOutputhead += "#\tModifications: Modifications added to the metabolite\n"
    sOutputhead += "#\tLinks: Links from the database\n"


    sOutputhead += "#\n"
    
    sOutputhead += "Filename\tScanNumber\tPrecursorMZ\tRetentionTime\tScanType\tRank\tParentMassError\tExplainedPeaks\tScore\tIdentifier\tName\tOriginalInChI\tIdentifiedInChI\tModifications\tLinks\n"

    ConfigureFile.close()
    return sOutputhead

def ReadFT2Files(sFT2_Filename, iDefault_Polarity, iDefault_Charge_State) :
# Read FT2 file, return all scan information
#sFT2_Filename: FT2 file name; iDefault_Polarity: 1 positive polarity, -1 negative polarity; iDefault_Charge_State: default charge state.

    FT2_File = open(sFT2_Filename)
    dProtonMass = 1.007825 # proton mass
    all_scans_list = []
    current_peaks_list = []
    dPrecursorMofZ = 0
    sScanNumber = "0"
    precursor_type = 0
    sScanHead = ""
    sRetentionTime = "NA"
    for each_line in FT2_File :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if not(each_line[0].isdigit()) :
            if (current_peaks_list) :
                #print current_peaks_list
                current_peaks_list.sort(key=lambda e:e[0])
                all_scans_list.append([current_peaks_list, dPrecursorMofZ - precursor_type*dProtonMass, sScanNumber, precursor_type, dPrecursorMofZ, sRetentionTime, sScanHead])
                current_peaks_list = []
                sRetentionTime = "NA"
                sScanHead = ""
            all_info_list = each_line.split("\t")
            if (all_info_list[0] == "S") :
                sScanHead += each_line+"\n"
                sScanNumber = all_info_list[1]
                dPrecursorMofZ = float(all_info_list[3])
            if (all_info_list[0] == "Z") :
                sScanHead += each_line+"\n"
                sZ = all_info_list[1]
                if (sZ == "+0") :
                    precursor_type = iDefault_Charge_State
                elif (sZ == "-0") :
                    precursor_type = (-1) * iDefault_Charge_State
                elif (sZ == "0"):
                    precursor_type = iDefault_Polarity * iDefault_Charge_State
                else :
                    precursor_type = int(sZ)
            if (all_info_list[0] == "I") and (len(all_info_list) >= 3):
                sScanHead += each_line+"\n"
                if (all_info_list[1] == "RetentionTime") :
                    sRetentionTime = all_info_list[2]
        else :
            all_peakinfo = each_line.split("\t")
            current_mofz = float(all_peakinfo[0])
            current_intensity = float(all_peakinfo[1])
            current_peaks_list.append([current_mofz, current_intensity, 0])
            #print current_peaks_list
    if current_peaks_list :
        current_peaks_list.sort(key=lambda e:e[0])
        all_scans_list.append([current_peaks_list, dPrecursorMofZ - precursor_type*dProtonMass, sScanNumber, precursor_type, dPrecursorMofZ, sRetentionTime])
    FT2_File.close()

    return all_scans_list

def main(argv=None):
#This is main function, i.e., the program starts here.
    if argv is None:
        argv = sys.argv
#Parse parameters
        [sFT2_Filename, sOutput_Filename, sConfigure_Filename, sAnnotation_Filename] = parse_options(argv)
        sFT2_basename = os.path.basename(sFT2_Filename)
#Read Configure file
        [compound_filename, iDefault_Polarity, iDefault_Charge_State, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, bBreakRing, iFragmentation_Depth, process_number, bRankSum, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list] = ReadConfigureFile(sConfigure_Filename)
# Read FT2 file
        all_scans_list = ReadFT2Files(sFT2_Filename, iDefault_Polarity, iDefault_Charge_State)
# Write the head part of .meb file
        sOutputhead = OutputFileHead(sConfigure_Filename)
        Output_File = open(sOutput_Filename, "w")
        Output_File.write(sOutputhead)
        Output_File.close()
# Write the head part of .AFT2 file
        Annotation_File = open(sAnnotation_Filename, "w")
        Annotation_File.write("Filename\tScanNumber\tPrecursorMZ\tRetentionTime\tScanType\tRank\tParentMassError\tExplainedPeaks\tScore\tIdentifier\tName\tOriginalInChI\tIdentifiedInChI\tModifications\tLinks\n")
        Annotation_File.write("m/z\tIntensity\tResolution\tBaseline\tNoise\tCharge\tSMILES\tProtonOffset\tm/zError\tFragmentationLevel\n")
        Annotation_File.close()
        bSpectrumDetails = False
#Analyze scans
        HandleAllRealHits(all_scans_list, compound_filename, sOutput_Filename, process_number, bSpectrumDetails, bBreakRing, bRankSum, iParentMassWindow_list, dMass_Tolerance_Parent_Ion, dMass_Tolerance_Fragment_Ions, iFragmentation_Depth, sFT2_basename, iPositive_Ion_Fragment_Mass_Windows_list, iNegative_Ion_Fragment_Mass_Windows_list, sAnnotation_Filename)



## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()


