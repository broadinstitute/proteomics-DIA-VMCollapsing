#!/usr/bin/env python
# coding: utf-8

# In[1]:


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[2]:


from alphamap.importing import import_data
from alphamap.preprocessing import format_input_data
from pyteomics import fasta

import re
import math
import os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

from statistics import mean,stdev,mode
import argparse
import sys

import yaml
from yaml.loader import SafeLoader

import Bio
from Bio import SeqIO

import copy


# In[3]:


wd = 'Z:\Helium_Tan\R02_PTMDIA\Exploris\Spectronaut\DDALibrary' + '\\'
with open(wd + 'input.yaml', 'r') as file:
    params = yaml.safe_load(file)


# In[5]:


#Unpack yaml file input parameters
file_name = params['main']['report']
fasta_path = params['main']['fasta_path']                  #File path that points directly to the file
engine = params['main']['engine']                          #Search engine
organism = params['main']['organism'] 

target_ptm = params[engine]['target_ptm']                  #User-indicated target PTM
experimental = params[engine]['experimental']              #Software-specifc list of annotations for experimental modifications that should be included in the search
heavies = params[engine]['heavies']                        #Software-specifc list of annotations for heavy modifications that should be included in the search
file_col = params[engine]['file']                          #Software-specific name of column that indicates file name
all_sequences = params[engine]['all_sequences']            #Column for modified sequences
prot_column = params[engine]['prot_column']                #Software-specific name of column containing protein IDs

ptm_phrases = params['phrases']                            
for k in ptm_phrases:
    if target_ptm in k:
        ptm_phrase = ptm_phrases[k]                        #PTM phrase is the assigned phrase that the target_ptm annotation is mapped to

threshold = params[engine]['ptm_confidence_threshold']     #Confidence localization threshold
keep = params[engine]['keep']                              #List of software-specific columns that are kept in the final collapsed report

print(threshold)


# In[6]:


if engine == 'SN':
    report = pd.read_csv(wd + file_name, delimiter ='\t', low_memory= False)    #Read in Spectronaut "Normal Report"
    if threshold != 0:
        report = report.loc[report['EG.PTMAssayProbability'] >= threshold]      #Apply site localization threshold to the report
    samples = (report[file_col]).unique()                                       #List of unique samples/ LC-MS runs that IDs may come from
    
elif engine == 'DIANN':                                                        
    report = pd.read_csv(wd+ file_name, sep = '\t')
    matrix = pd.read_csv(wd + 'report.pr_matrix.tsv', sep = '\t', low_memory = False)  
    samples = (report[file_col]).unique()


# # The following functions are used by both SN and DIA-NN.

# In[7]:


def find_prot_positions(sequence, prot_positions, engine, target_ptm, heavies):
    '''
    Find positions of target PTM and heavies in a sequence

    Args: 
        sequence: Modified sequence with experimental mods removed
        prot_positions: Positions of peptide on protein
        engine
        target_ptm
        heavies
        
    Returns: List of annotated target PTM and heavy PTM protein positions
   
    '''
    all = list(heavies)
    all.append(target_ptm)

    if engine == 'SN':
        all_regex = '[A-Z]('+ '|'.join(['\\[\\' + x + '\\]' for x in all]) + ')'        #Generates regex expression that searches for all target ptm and heavy annotations and the amino acid
        prot_positions = re.split(r';|,', prot_positions)                               #One odd case in SN output where IDs were split by comma instead of semicolon, this covers both cases
    
    if engine == 'DIANN':
        all_regex = '[A-Z]('+ '|'.join(['\\(' + x + '\\)' for x in all]) + ')'

    finder = re.finditer(all_regex, sequence)                                       #Find matches within sequence
    naked = re.sub(all_regex, 'X', sequence)                                        #Replace all matches with a single character to find their position in naked sequence

    annotations = [match.group() for match in finder]                               #All annotations with amino acid residue
    pep_positions = [match.start() for match in re.finditer('X', naked)]            #Peptide positions of matches (0-indexed)
        
    full = []
    
    for p in prot_positions:                                                        #PEP_positions relating to different protein IDs are separated by ';'
        
        protein_pos = [str(int(p) + x) for x in pep_positions]                      #Add the peptide position of the match to the protein position of peptide --> protein position of the match
        ann = [''.join(item) for item in zip(*[annotations, protein_pos])]          #Combine the annotation (regex match) with the updated protein position
        full.append(ann)
        
    return(full)


# In[8]:


def remove_experimental_mods(sequence):   
    """
    Removes annotations relating to N-term acetylation, oxidation of methionine, and carbamidomethylation. Also removes '_' flanking sequences to ease indexing in Spectronaut.

    Args:
        sequence: Sequence containing experimentally-introduced mod annotations listed above.

    Returns: Peptide sequence without experimentally-introduced mod annotations listed above.
    """
    for item in experimental:
        sequence = sequence.replace(item, '') #Remove carbamidomethylation, N-term acetylation, M-oxidation
        
    return(sequence)


# # DIA-NN specific functions.

# In[9]:


def convert_mods(sequence):  #DIA-NN
    '''
    Converts modifications from alphamap output to match DIA-NN modifications.

    Args: 
        sequence: Peptide sequence with alphamap modifications
        
    Returns: Peptide sequence with DIA-NN modifications
    '''
    
    dict = {

        '\[Acetyl (.*?)\]': '(UniMod:1)',
        '\[Amidated (.*?)\]': '(UniMod:2)',
        '\[Carbamidomethyl (.*?)\]': '(UniMod:4)',
        '\[Carbamyl (.*?)\]' : '(UniMod:5)',
        '\[Deamidation (.*?)\]' : '(UniMod:7)',
        '\[Phospho (.*?)\]': '(UniMod:21)',
        '\[Dehydrated (.*?)\]' : '(UniMod:23)',
        '\[Pyro-carbamidomethyl (.*?)\]':'(UniMod:26)',
        '\[Glu->pyro-Glu\]' : '(UniMod:27)',
        '\[Gln->pyro-Glu\]' : '(UniMod:28)',
        '\[Cation:Na (.*?)\]' : '(UniMod:30)',
        '\[Methyl (.*?)\]': '(UniMod:34)',
        '\[Oxidation (.*?)\]': '(UniMod:35)',
        '\[Dimethyl (.*?)\]': '(UniMod:36)',
        '\[Trimethyl (.*?)\]' : '(UniMod:37)',
        '\[Sulfo (.*?)\]': '(UniMod:40)',
        '\[Cys-Cys\]': '(UniMod:55)',
        '\[GlyGly (.*?)\]': '(UniMod:121)',
        '\[Delta:H(2)C(2) (.*?)\]': '(UniMod:254)',
        '\[Cysteinyl\]': '(UniMod:312)',
        '\[Trioxidation (.*?)\]': '(UniMod:345)',
        '\[Hydroxyproline\]': '(UniMod:408)',
        '\[Dioxidation (.*?)\]':'(UniMod:425)',
        '\[Dethiomethyl (.*?)\]': '(UniMod:526)',
        '\[QQTGG (.*?)\]' : '(UniMod:877)'

    }

    for x in dict.keys():
        sequence = re.sub(x, dict[x], sequence)

    return(sequence)


# In[10]:


def alphabetize(proteins):  #DIA-NN
    """
    Alphabetizes order of protein IDs in DIA-NN matrix 'Protein.Ids' column. This is needed for the merging of the alphamap output to the matrix.

    Args: 
        proteins: string of protein IDs separated by ';'

    Returns: String of alphabetized protein IDs separated by ';'
    """
    ls = sorted(proteins.split(';'))    #Split and sort
    str = ';'.join(ls)                  #Rejoin with semicolon into str

    return(str)


# In[11]:


def preprocess(report, threshold, file_name, fasta_path):  #DIA-NN
    """
    Localizes report for site confidence >= threshold, exports localized report. Performs alphamap processing of report to get peptide position and PTM information. Adjusts PTM annotations.

    Args: 
        report: Long-form DIA-NN report
        threshold
        file_name
        fasta_path

    Returns: {'localized_report': Localized dataframe report, 'formatted_data': alphamap output} for further processing of DIA-NN data.
    """
    #Localize report, export localized report
    if threshold != 0:
        localized_report = report[report['PTM.Site.Confidence'] >= threshold]
        localized_report.to_csv(wd + 'Localized_Report.tsv', sep = '\t', index = False)
        report_path = wd + 'Localized_Report.tsv'
        
    else: 
        localized_report = report
        report_path = wd + file_name

    print('Report with confidently localized entries exported.')

    #Prepare report and fasta for Alphamap processing
    diann_report = import_data(report_path)
    fasta_file = fasta.IndexedUniProt(fasta_path)
    print('Prepared for alphamap processing')

    #Remove contaminants and heavy annotations (does not remove full sequences containing heavies)
    diann_report = diann_report[diann_report["all_protein_ids"].str.contains("contaminant") == False].reset_index(drop=True)
     
    #Removes heavy annotations from sequence (but keeps sequence) because they cannot be processed by alphamap
    for h in heavies:
        diann_report = diann_report.replace('\('+ h + '\)','', regex = True)

    #Run report throught alphamap to get formatted data with peptide locations
    formatted_data = format_input_data(df=diann_report, fasta=fasta_file, modification_exp = r'\[.*?\]')                        #Run through alphamap to get peptide locations on protein
    formatted_data.rename(columns={'modified_sequence': 'Modified.Sequence', 'all_protein_ids': "Protein.Ids"}, inplace=True)   #Rename 'all_protein_ids' column to 'Protein IDs to match the matrix it will be merged with

    #Update peptide start and end positions to be 1-indexed instead of 0-indexed
    formatted_data['start'] = formatted_data['start'] + 1
    formatted_data['end'] = formatted_data['end'] + 1

    #Change annotations in the formatted data to match matrix it will be merged with
    formatted_data["Modified.Sequence"] = formatted_data["Modified.Sequence"].apply(convert_mods)
    
    formatted_data = formatted_data.sort_values('unique_protein_id')        #Sort unique protein IDs so when STY locations get appended, they match order of alphabetized protein IDs
    
    return({'localized_report': localized_report,'formatted_data': formatted_data})


# In[12]:


def condense_formatted_data(formatted_data):  #DIA-NN   
    """
    Adds PTM residue protein positions and condenses alphamap output such that PTM protein locations are given for all protein IDs a sequence maps to in a single row.
    
    Args: 
        formatted_data: Alphamap output where each entry corresponds to a unique protein ID a sequence maps to.
    
    Returns: Formatted output collapsed with PTM protein locations for all protein IDs a sequence maps to.
    
    """
    
    #The Alphamap report creates individual entry for every protein ID a sequence is mapped to, here we group by 'Protein.ids' (all protein IDs a sequence maps to) and 'Modified Sequence'
    condense = ['start','end']
    list_cols = {col: lambda x: list(x) for col in condense}                                              #Condenses entries from individual rows to a list of PTM protein locations
    
    keep = ['Protein.Ids','Modified.Sequence','naked_sequence','PTMtypes','PTMsites']
    other_cols = {col: 'first' for col in keep}                                                           #Keeps first of all other entries because they are identical across all rows of the same group

    col_map = {**other_cols, **list_cols}                                                                 #Determines how grouping should be done

    #Aggregates data according to aggregation map, such that each row is a unique combination of modified sequence and LIST of protein IDs the sequence maps to
    formatted_data = formatted_data.groupby(['Protein.Ids','Modified.Sequence']).agg(col_map).reset_index(drop=True) 
    
    return(formatted_data)


# In[13]:


def localize_matrix(localized_report, matrix): #DIA-NN
    """
    Filters DIA-NN matrix file-specific quant for sequences that are only in the localized report.
    
    Args: 
        localized_report: Long-form DIA-NN report with site localization threshold applied
        matrix: wide-form DIA-NN report without site localization filtering
    
    Returns: Matrix where intensities for each sample are only reported on file-specific localized sequences.
    
    """
    
    localized = {}                         #{Sample: List of unique localized sequences}
    for s in samples:           
        one = localized_report[localized_report['File.Name'] == s]   #Only entries from that run in report
        if s not in localized:                   #Key is sample name, list of unique modified sequences that are localized in that run
            localized[s] = None
            localized[s] = one['Modified.Sequence'].unique()
            
        else:
            localized[s] = one['Modified.Sequence'].unique()
            
            
    for index, row in matrix.iterrows():          #Only iterate through matrix once
        seq = row['Modified.Sequence']
        
        for s in samples:         
            if seq not in localized[s]:
                matrix.at[index,s] = None         #If the sequence in the matrix is not localized, replace intensity with None
                
    return(matrix)


# In[14]:


def merge(formatted_data, localized_matrix, heavies):  #DIA-NN
    """
    Merges  alphamap output with protein positions to the filtered intensity matrix, which does not have protein positions. 
    
    Args: 
        formatted_data: Condensed alphamap output
        localized_matrix: Site-localized DIA-NN matrix
        heavies
    
    Returns: Merged ouput of intensity matrix with protein position information appended.

    """
    localized_matrix['Precursors_Heavies'] = localized_matrix['Modified.Sequence']       #Maintains record or heavies included in sequence, this is confirmed to be a deep copy
    
    for h in heavies:
        h = '(' + h + ')'
        localized_matrix['Modified.Sequence'] = localized_matrix.apply(lambda x: x['Modified.Sequence'].replace(h,''),axis=1) #Removes heavy annotations from this sequence so that it will match the sequences in alphamap output (which cannot contain heavy annotations)
        
    localized_matrix[prot_column] = localized_matrix.apply(lambda x: alphabetize(x['Protein.Ids']),axis=1)   #Alphabetizes protein IDs in matrix to match order of protein IDs in report it will be merged with
    
    merged_report = localized_matrix.merge(formatted_data, "inner", on=["Modified.Sequence", "Protein.Ids"]) #Merges two data frames based on common protein IDs and modified sequences that do not contain heavies

    return(merged_report)


# In[15]:


def generate_collapsed_report_DIANN(report, ptm_phrase, target_ptm, prot_column, experimental, engine, heavies, keep):
    """
    Generates collapsed PTM report where each entry represents a unique phosphosite. A phosphosite
    is defined as a uniquely phosphorylated peptide sequence. Differentially heavy-modified peptides are listed as unique phosphosite entries.

    Args:
        report: Wide-form dataframe ready to be collapsed
        ptm_phrase
        target_ptm
        prot_column
        experimental
        engine
        heavies
        keep

    Returns: Collapsed PTM site report. 

    """       

    target = report.loc[report['Modified.Sequence'].str.contains(target_ptm)]            #Only target PTM sequences
    

    collapse_keys = []
    targets_unannotated = []

    for index, row in target.iterrows():
        
        protein_ids = row[prot_column]
        sequence = remove_experimental_mods(row['Precursors_Heavies'])                    #This column contains sequences containing heavy mods
        pep_prot_positions = row['start']
        
        #Target and heavy PTM positions
        ptm_prot_positions = find_prot_positions(sequence, pep_prot_positions, engine, target_ptm, heavies)

        #Integer-unannotated target positions used by find_flanking function
        unannotated = [[re.sub("\(.*?\)", "", s) for s in pos if target_ptm in s] for pos in ptm_prot_positions]  #Keeps only sequences with target mod (ex: 'UniMod:21'), and removes integer annotations in brackets
        targets_unannotated.append(unannotated)
        
        #S(UniMod:21)69 --> S69

        #Generate collapse key: #Collapse key contains: Protein IDs, PTM protein locations, heavy mod protein locations
        k = (protein_ids,str(ptm_prot_positions))         
        collapse_keys.append(k)

    target['Collapse'] = collapse_keys                                     #Creates new column with collapse key values
    target[ptm_phrase +'ProteinLocations'] = targets_unannotated


    keep.append(ptm_phrase +'ProteinLocations')                            #Add additional column that should be kept in final output
    
    #Prepare columns for grouping
    other_cols = {col: 'first' for col in keep}
    quants = {col: lambda x: x.sum(skipna=False) for col in samples}       #Sum all the values that are not NA

    col_map = {**other_cols, **quants}

    collapsed = target.groupby(['Collapse'], as_index=False).agg(col_map)  #Group by collapse key, sum values by sample

    keep.extend(samples)                                                   #Append file names as column titles, these columns contain file-specific quant values
    collapsed_report = collapsed[keep]                                     #Only keep files that are assigned to be kept 
    
    collapsed_report.to_csv(wd + 'VM_CollapsedReport.tsv', sep = '\t', index = False )
    
    return(collapsed_report)  


# # Spectronaut specific functions.

# In[16]:


def generate_collapsed_report_SN(report, ptm_phrase, prot_column, experimental, samples, keep):
    
    """
    Takes normal report from Spectronaut output and generates collapsed PTM report where each entry represents a unique phosphosite. A phosphosite
    is defined as a uniquely phosphorylated peptide sequence. Differentially heavy-modified peptides are listed as unique phosphosite entries.

    Args:
        report: Spectronaut Normal Report. Must contain the following columns for script to function:
            'R.FileName'
            'FG.IntMID'
            'EG.IntPIMID'
            'PG.ProteinAccessions'
            'EG.ProteinPTMLocations'
            'FG.Quantity'
            'EG.PTMAssayProbability'
            'EG.Cscore'
        ptm_phrase
        prot_column
        experimental
        samples
        keep

    Returns: Collapsed PTM site report.
    """

    target = report.loc[report['EG.ModifiedSequence'].str.contains(ptm_phrase)]    #Only entries with target PTM are filtered

    
    collapse_keys = []
    targets_unannotated = []
    
    for index, row in target.iterrows():
        
        sequence = remove_experimental_mods(row['FG.IntMID'])  
        pep_prot_positions = row['PEP.PeptidePosition']
        protein_ids = row[prot_column]
        
        #Target and heavy PTM positions
        ptm_prot_positions = find_prot_positions(sequence, pep_prot_positions, engine, target_ptm, heavies)
        
        #Integer-unannotated target positions used by find_flanking function
        unannotated = [[re.sub("\[.*?\]", "", s) for s in pos if target_ptm in s] for pos in ptm_prot_positions]  #Keeps only sequences with target mod (ex: '+80'), and removes integer annotations in brackets
        targets_unannotated.append(unannotated)
        
        #Collapse keys
        k = (protein_ids,str(ptm_prot_positions))    #Collapse key contains: Protein IDs, target PTM protein locations (this way even sequences with other modifications will still collapse if they have the same target PTM positions), heavy mod protein locations
        collapse_keys.append(k)

        
    #Add columns to target dataframe
    target['Collapse'] = collapse_keys                                     #Creates new column with collapse key values
    target[ptm_phrase +'ProteinLocations'] = targets_unannotated           #Creates new column with integer-unannotated target PTM protein positions
    
    first_column = target.pop('Collapse')
    target.insert(0, 'Collapse', first_column)
    
    #For rows with the same collapse key in the same run, create a new column with the summed quantity across those rows. That way when picking a representative row, the column will still have summed quant.
    target['QuantSum'] = target.groupby(['Collapse','R.FileName'])['FG.Quantity'].transform('sum')     
    
    #Get a dataframe that's collapsed by the combination of collapse key AND file
    localization = target.groupby(['Collapse', 'R.FileName'])['EG.PTMAssayProbability'].transform(max) == target['EG.PTMAssayProbability']  #Keep row with highest PTM localization confidence
    target = target[localization]

    cscore = target.groupby(['Collapse', 'R.FileName'])['EG.Cscore'].transform(max) == target['EG.Cscore']                                  #Keep row with the highest CScore
    target = target[cscore]
    
    target = target.groupby(['Collapse', 'R.FileName']).first().reset_index()     #If rows in a group have same PTM assay probability and CScore, keep first row                                                           

    #Create a dataframe with wide-form intensities
    pivoted = target.pivot(index='Collapse', columns='R.FileName', values='QuantSum') 
    
    #Create a dataframe that is collapsed just by collapse key
    collapsing_localization = target.groupby(['Collapse'])['EG.PTMAssayProbability'].transform(max) == target['EG.PTMAssayProbability']
    collapsed = target[collapsing_localization]

    collapsing_cscore = collapsed.groupby(['Collapse'])['EG.Cscore'].transform(max) == collapsed['EG.Cscore']
    collapsed = collapsed[collapsing_cscore]

    collapsed = collapsed.groupby(['Collapse']).first().reset_index()
    
    
    #Merge wide-form intensities with collapsed long form report to create a wide form collapsed report
    collapsed_report = collapsed.merge(pivoted, "inner", on=["Collapse"])

    #Select columns to keep in the ouput
    keep.append(ptm_phrase +'ProteinLocations')

    #Append file ptm_phrases as column titles, these columns contain file-specific quant values
    keep.extend(samples)

    #Generate new collapsed report with selected columns only
    collapsed_report = collapsed_report[keep]
    collapsed_report.to_csv(wd + 'VM_CollapsedReport.tsv', sep='\t', index=False)

    return(collapsed_report)


# # At this point, collapsed_report is in the same format for either software and can be summarized.

# In[17]:


def summarize_non_missing(setlist, collapsed_report, samples, sample, threshold, target_ptm, ptm_phrase, experimental):
    """
    Reports summary stats of sequences present across all samples searched together.

    Args:
        setlist: list of sets, one set of sequences per sample
        collapsed_report: Collapsed PTM site report
        samples: List of all names
        sample: Column name that will be used for the summary document
        threshold
        target_ptm
        ptm_phrase
        experimental

    Returns: Dictionary with values for
        % precursors containing target modification 
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique PTM sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed PTM-sites, including heavy mods

    """
    
    intersection = set.intersection(*setlist)                                 #Sequences that are found across all samples
    num_modified_sequences = len(intersection)

    ptm_sequences = [x for x in intersection if target_ptm in x]              #Only PTM-containing precursors
    enrichment = (len(set(ptm_sequences)) / num_modified_sequences) * 100     #PTM-containing precursors / all precursors

    ptm_no_experimental = list(map(remove_experimental_mods, ptm_sequences))  #Remove experimentally-introduced modifications, keep heavy mods
    unique_ptm_peptides = len(set(ptm_no_experimental))                       #This will represent the number of phosphopeptides

    no_missing = collapsed_report.dropna(thresh = threshold, subset = samples)                    #Drop sites that have missing values over the threshold allowed
    ptm_sites = len(no_missing)

    
    ret_df = pd.DataFrame({'Run': sample, '%_Precursors_with_' + ptm_phrase: enrichment,
             'num_modified_sequences': num_modified_sequences,
             'num_' + ptm_phrase + '_peptides': unique_ptm_peptides,
             'num_' + ptm_phrase + '_sites': ptm_sites}, index=[0])

    return (ret_df)


# In[18]:


def summarize_any(full, collapsed_report, sample, engine, target_ptm, all_sequences, ptm_phrase, experimental):
    """
    Generates summary stats for any data frame (could be subset of original report)

    Args:
        full: Data frame that has not been collapsed
        collapsed_report: Collapsed PTM site report
        sample: Column name that will be used for the summary document
        engine
        target_ptm
        all_sequences
        ptm_phrase
        experimental

    Returns: Dictionary with values for
        % precursors containing target modification
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique PTM sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed PTM-sites, including heavy mods

    """
    num_modified_sequences = len(full[all_sequences].unique())                          #All precursors, including experimentally-introduced mods and heavy mods
    
    #Enrichment efficiency
    ptm_sequences = [x for x in full[all_sequences] if target_ptm in x]                 #Only PTM-containing precursors
    enrichment = (len(set(ptm_sequences)) / num_modified_sequences) * 100         #PTM-containing precursors / all precursors
    
    #Number of unique PTM modified sequences
    ptm_no_experimental = list(map(remove_experimental_mods, ptm_sequences))      #Remove experimentally-introduced modifications, keep heavy mods
    unique_ptm_peptides = len(set(ptm_no_experimental))                           #This will represent the number of PTM-peptides

    #Number of phosphosites per sample
    if sample != 'Combined':
        select = collapsed_report[sample].values.tolist()
        res = [i for i in select if math.isnan(i) == False]
        ptm_sites = len(res)

    #Length of data frame is the number of phosphosites identified across all samples searched together
    if sample == 'Combined':
        ptm_sites = len(collapsed_report)
        
    #Create data frame to add to summary report
    ret_df = pd.DataFrame({'Run': sample, '%_Precursors_with_' + ptm_phrase: enrichment,
        'num_modified_sequences': num_modified_sequences,
        'num_' + ptm_phrase + '_peptides': unique_ptm_peptides,
        'num_' + ptm_phrase + '_sites': ptm_sites}, index=[0])
    
    
    return (ret_df)


# In[19]:


def summarize(report, collapsed_report, samples, file_col, all_sequences, ptm_phrase, threshold):
    """
    Reports summary stats across all runs searched together.

    Args:
        report: Initial report from either software (Spectronaut Normal report, DIANN report.tsv)
        collapsed_report: Collapsed PTM site report
        samples: List of all file names
        file_col
        all_sequences
        ptm_phrase
        threshold

    Returns: Datafrane with following values for each sample, combined data set, and intersection across samples:
        % precursors containing target modification (STY)
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique STY sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed phosphosites, including heavy mods

    """
    print('Summarizing data...')
    data = pd.DataFrame(columns= ['Run', '%_Precursors_with_' + ptm_phrase, 'num_modified_sequences','num_' + ptm_phrase + '_peptides','num_' + ptm_phrase + '_sites'])       #Initialize datarframe with columns
    
    summary_list = list(samples).copy()
    summary_list.append('Combined')
    

    non_missing = []                                        #Collect precursors found in each sample individually, will be list of sets
    row_list = []                                           #Collect dataframes to concatenate at the end for final summary dataframe

    for s in summary_list:
        if s != 'Combined':
            one = report[report[file_col] == s]                                                                          #This should be included in a summarize_any and in a single function
            
        else:
            one = report
            
        ret = summarize_any(one, collapsed_report, s, engine, target_ptm, all_sequences, ptm_phrase, experimental)       #Returns a row of the dataframe containing stats for that sample
        row_list.append(ret)

        set_IDs = set(one[all_sequences])  
        non_missing.append(set_IDs)               #Add to non_missing, a list of sets of modified sequences from each sample
        

    #Add row for intersection of data across all samples searched together
    intersection = summarize_non_missing(non_missing, collapsed_report, samples, 'Intersection', len(samples), target_ptm, ptm_phrase, experimental)                      #Sites quantified across all samples in dataset
    fifty_perc_complete = summarize_non_missing(non_missing, collapsed_report, samples, 'Fifty Percent Complete', len(samples)/2, target_ptm, ptm_phrase, experimental)   #Sites quantified across at least 50% of all samples in dataset

    row_list.append(intersection)
    row_list.append(fifty_perc_complete)
    
    summary_df = pd.concat(row_list)                         #Combined all rows into one single dataframe
    
    #Export summary to working directory
    summary_df.to_csv(wd + 'VMSummary_' + str(threshold) + '.tsv', sep='\t', index=False)
    print('Summarizing complete')
    
    return(data)


# In[20]:


def find_flanking_phosphositeplus(collapsed_report, organism): 
    """
    Provides flanking sequences of 7 amino acids around target PTM residue(s) for all protein IDs a sequence has flanking information for.

    Args:
        collapsed_report: Collapsed PTM site report
        organism

    Returns: Collapsed report with flanking sequence information added

    """
    
    meta = pd.read_csv('Z:\\Helium_Tan\\PhosphositePlus\\Phosphorylation_site_dataset', sep = '\t')
    meta = meta[meta['ORGANISM'] == organism]   #Filter for organism, typically human

    dict = {}

    for index, row in meta.iterrows():
        prot_id = row['ACC_ID']
        residue = row['MOD_RSD']
        flank_seq = row['SITE_+/-7_AA'].upper()

        res = residue.split('-')[0]
        k = (prot_id, res)

        if k not in dict:
            dict[k] = None
            dict[k] = flank_seq

        else:
            dict[k] = flank_seq


    flanks_column = []         #Global list
    
    #Iterate through entries and add flanking sequences for each protein ID and site
    for index, row in collapsed_report.iterrows():
        
        prot_id = row[prot_column].split(';')                             #Get all protein IDs a sequence is mapped to
        locations = list(row[ptm_phrase +'ProteinLocations'])             #PTM locations in protein

        row_flanks = {}                                             #Each row/entry will have a dictionary if there are flanking sequences in phosphosite table, with protein ID as key and flanking sequence(s) as value(s)
        
        #Find flanks for each protein ID listed in entry
        for i in range(0, len(prot_id)):                            
            
            prot_flanks = []         #If there are multiple sites on the same sequence, they will all be appended here
            
            sites = locations[i]     #Can contain multiple sites
            id = prot_id[i]

            contains = False
            for s in sites:                                          #Each PTM-site on sequence associated with that protein ID
                rep_k = (id,s)                                       #(Protein ID, PTM Location on that protein ID)

                if rep_k in dict.keys():                             #Phosphosite table has mapping seq
                    contains = True
                    prot_flanks.append(dict[rep_k])
                    
                else:
                    prot_flanks.append('')                                     #Sometimes if a seq is doubly phosphorylated only one site will be mapped in phosphosite table


            if contains:
                if id not in row_flanks:
                    row_flanks[id] = None
                    row_flanks[id] = prot_flanks
                else:
                    row_flanks[id] = prot_flanks

        if len(row_flanks) > 0:
            flanks_column.append(row_flanks)
        else:
            flanks_column.append(None)


    collapsed_report['7AA_Flanking_PhosphositePlus'] = flanks_column           #Add a column for PhosphoSitePlus data to collapsed report


# In[21]:


def find_flanking(collapsed_report, prot_column):
    """
    Adds flanking sequence to the collapsed report based on the user-provided fasta file.
    
    Args:
        collapsed_report: PTM-collapsed report
        prot_column: Software-specific column containing protein IDS
    
    Returns: Collapsed report with flanking sequences added

    """
    
    parsed_fasta = SeqIO.parse(fasta_path, 'fasta')    #Parse the fasta file

    fasta_dict = {}

    for seq_record in parsed_fasta:                    #Create dictionary {Protein ID: protein sequence} based on fasta file
        header = seq_record.id
        prot_id = header.split('|')[1]
        seq = seq_record.seq
    
    
        fasta_dict[prot_id] = seq
    
    #Fasta dictionary is ready, now add flanking sequences
    
    flanks_column = []     #Global list
    
    for index, row in collapsed_report.iterrows():      #Iterate through collapsed report
        prot_ids = row[prot_column].split(';')
        ptm_locations = list(row[ptm_phrase +'ProteinLocations'])


        row_flanks = {}
        for i in range(0, len(prot_ids)):   

            prot_flanks = []                            #Sequences for all target PTM sites in the mapped protein ID

            modified = ptm_locations[i]                 #target PTM protein locations
            id = prot_ids[i]

            for x in modified:
                char = x[0]
                pos = int(x[1:]) - 1                    #0-indexed position
                seq = fasta_dict[id]

                flank_start = pos - 7                 
                if flank_start < 0:                     #If there aren't 7 amino acids preceding the target residue, just make the start index the start of the sequence
                    flank_start = 0                   

                flank_end = pos + 8
                if flank_end > len(seq):
                    flank_end = len(seq) + 1            #If there aren't 7 amino acids proceeding the target residue, just make the ending index the end of the sequence

                flanking = str(seq[flank_start:flank_end])
                prot_flanks.append(flanking)

            row_flanks[id] = prot_flanks
        
        flanks_column.append(row_flanks)
        
    collapsed_report['7AA_Flanking'] = flanks_column
    
    collapsed_report.to_csv(wd + 'VM_Collapsedreport_Flanks.tsv', sep = '\t')
    print('Flanking sequences added')


# # Run all the code 

# In[22]:


#Runs through all code sequentially to walk through steps individually as required for software-specific processing
def run():
    
    if engine == 'DIANN':

        preprocessed = preprocess(report, threshold, file_name, fasta_path)   #return a dictionary instead of a list
        localized_report = preprocessed['localized_report']
        
        alphamap_uncondensed = preprocessed['formatted_data']
        alphamap_condensed = condense_formatted_data(alphamap_uncondensed)
        
        localized_matrix = localize_matrix(localized_report, matrix)
        final_report = merge(alphamap_condensed, localized_matrix, heavies)                      
        collapsed_report = generate_collapsed_report_DIANN(final_report, ptm_phrase, target_ptm, prot_column, experimental, engine, heavies, keep)
        summarize(localized_report, collapsed_report, samples, file_col, all_sequences, ptm_phrase, threshold)
        
    if engine == 'SN':
        
        collapsed_report = generate_collapsed_report_SN(report, ptm_phrase, prot_column, experimental, samples, keep)
        summarize(report, collapsed_report, samples, file_col, all_sequences, ptm_phrase, threshold)
        
    if ptm_phrase =='Phospho':
        find_flanking_phosphositeplus(collapsed_report, organism)
    find_flanking(collapsed_report, prot_column)
        
run()

