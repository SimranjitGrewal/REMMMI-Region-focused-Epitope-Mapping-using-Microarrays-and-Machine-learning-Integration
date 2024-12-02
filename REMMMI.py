##this program contains artifactual elements and temporary functions that have been commented out
# They have been included and kept because in some cases it has been employed to generate some figures

#a minimum requirement to employ this program is an excel file or CSV file containing 2 columns, the first colum

global feature_weights
global peptide_skip
global frame_size #frame size variable plus one is tree depth 
print("a minimum requirement to employ this program is an excel file or CSV file containing 2 columns for each protein")
print("Column 1 must be labelled 'Peptide' in the first row and below it should be peptide sequences")
print("Column 2 must be labelled 'Raw Median' in the first row and below it should be binding scores corresponding to same-row sequence")
print("Before running, generate per-amino acid S4PRED for the concattenated sequence of all peptides in the array for each protein")
print("S4PRED function uses overlap values to generate an equal number of variables corresponding to total peptides")
####IMPORTS###############################################################################################################################
from ast import List
import pandas as pd
from itertools import product
#from sklearn_genetic.space import Continuous, Categorical, Integer
#from sklearn_genetic.plots import plot_fitness_evolution, plot_search_space
import numpy as np
import matplotlib.pyplot as plt
import sklearn
from sklearn import tree
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import statistics
from openpyxl import Workbook
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
import pickle
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_regression
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression
import graphviz
import dtreeviz

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedKFold

from sklearn.neural_network import MLPRegressor
#from sklearn_genetic.schedules import InverseAdapter

#import torch
#import Bio
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import RFE
import statistics


#import subprocess
#import os



#yes= r'"C:\Users\sgrew\source\repos\pytorch\pytorch\env\s4pred", python run_model.py, --device gpu, --save-files, --outdir, "C:\Users\sgrew\Desktop\codefolder", --save-by-idx, "C:\Users\sgrew\Desktop\codefolder\file folder for loaded 7th pep counts corrected averages\FCR3.fasta"'
#os.system(yes)
#print("yes")

####FUNCTIONS#############################################################################################################################

####FUNCTIONS############################# SPECIFICALLY FOR BUILDING PHASE 1 AND 2 DF

#take excel and open it as an array only importing peptides as strings
#takes columns and converts it to a list of strings
#deletes old array
#converts strings to text files in folder with index 0 to X as their file name
#converts


#features! This block of functions generates features from the included dataset
#to swap which features are included or excluded, comment them in or out in the function 'initialize features'

def find_aa_pairs(peptide_data,aa_pairs,peptide_skip): #this will go through every dataframe in our list of data frames and then every column and performs certain functions
    b=0
    peptide_data[aa_pairs]=0 #this had to be added because the values were already zero until I tried looping the function, then I had to force it to be the case for the function call. This is not true of other functions that due not initialize because those are not within the While loop for mutliple function calls
    while b < len(peptide_data):
        d=[] #empty list of all the peptide pairs in the chain
        for c in range(len(peptide_data.iloc[b,peptide_data.columns.get_loc('Peptide')])-1-(int(peptide_skip))):#the range here indicates which peptides should be counted, peptide skip here is a variable apramterer set by user
            d.append(peptide_data.iloc[b,peptide_data.columns.get_loc('Peptide')][c]+peptide_data.iloc[b,peptide_data.columns.get_loc('Peptide')][c+1]) #builds the peptide pair list
        for e in range(len(d)):
            peptide_data.iloc[b,aa_pairs.index(d[e])+peptide_data.columns.get_loc('Raw Median')+1]=peptide_data.iloc[b,aa_pairs.index(d[e])+peptide_data.columns.get_loc('Raw Median')+1]+1  #this will go through d and increment the correct pair column by 1 there is a +1 because of data frame column offset
        b=b+1
    print(peptide_data)
    return peptide_data #aa pairs

def calculate_hydropathy_peptide(optimized_dataset):
    #From the paper our scale was obtained - Characterizing Hydropathy of Amino Acid Side Chain in a Protein Environment by Investigating the Structural Changes of Water Molecules Network
    hydropathy_list=[]
    hydropathy={'A':0.47,'G':0.8,'L':0.79,'M':0.38,'F':0.24,'W':0.65,'K':2.76,'Q':1.26,'E':3.32,'S':1.46,'P':0.55,'V':0.46,'I':0.00,'C':0.59,'Y':0.83,'H':2.29,'N':2.03,'R':1.60,'D':3.20,'T':0.66}
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
     peptide_hydropathy=0
     for b in range(len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])): #b is the index of the amino acid in row 'a' peptide
      peptide_hydropathy= peptide_hydropathy + hydropathy[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]] #pulls each letter forma  pepetide and attaches it to its respective position and then finds on the list of features what the index of that combination is and then increments that value by 1 in the correct locaiton. A explores each peptide and B explore every amino acid in each peptide
     hydropathy_list.append(peptide_hydropathy)
    optimized_dataset["hydropathy"]=hydropathy_list
    return optimized_dataset #hydropathy

def calculate_charge(optimized_dataset): #net charge
    charge_list=[]
    charge= {'A':0,'G':0,'L':0,'M':0,'F':0,'W':0,'K':1,'Q':0,'E':-1,'S':0,'P':0,'V':0,'I':0,'C':0,'Y':0,'H':0,'N':0,'R':1,'D':-1,'T':0}
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
     peptide_charge=0
     for b in range(len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])): #b is the index of the amino acid in row 'a' peptide
      peptide_charge= peptide_charge + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]] #pulls each letter forma  pepetide and attaches it to its respective position and then finds on the list of features what the index of that combination is and then increments that value by 1 in the correct locaiton. A explores each peptide and B explore every amino acid in each peptide
     charge_list.append(peptide_charge)
    optimized_dataset['charge_At_7ph']=charge_list
    return optimized_dataset #charge

def calculate_charge_plus(optimized_dataset): #total positive residue count
    charge_list=[]
    charge= {'A':0,'G':0,'L':0,'M':0,'F':0,'W':0,'K':1,'Q':0,'E':0,'S':0,'P':0,'V':0,'I':0,'C':0,'Y':0,'H':0,'N':0,'R':1,'D':0,'T':0}
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
     peptide_charge=0
     for b in range(len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])): #b is the index of the amino acid in row 'a' peptide
      peptide_charge= peptide_charge + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]] #pulls each letter forma  pepetide and attaches it to its respective position and then finds on the list of features what the index of that combination is and then increments that value by 1 in the correct locaiton. A explores each peptide and B explore every amino acid in each peptide
     charge_list.append(peptide_charge)
    optimized_dataset['positive_charge_At_7ph']=charge_list
    return optimized_dataset #charge positives 1

def calculate_charge_minus(optimized_dataset): #total negative residue count
    charge_list=[]
    charge= {'A':0,'G':0,'L':0,'M':0,'F':0,'W':0,'K':0,'Q':0,'E':-1,'S':0,'P':0,'V':0,'I':0,'C':0,'Y':0,'H':0,'N':0,'R':0,'D':-1,'T':0} #only negatives
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
     peptide_charge=0
     for b in range(len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])): #b is the index of the amino acid in row 'a' peptide
      peptide_charge= peptide_charge + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]] #pulls each letter from a pepetide and attaches it to its respective position and then finds on the list of features what the index of that combination is and then increments that value by 1 in the correct locaiton. A explores each peptide and B explore every amino acid in each peptide
     charge_list.append(peptide_charge)
    optimized_dataset['negative_charge_At_7ph']=charge_list
    return optimized_dataset #charge negatives 1

def calculate_aminoacid_counts(optimized_dataset,peptide_skip): #simple function counts the occurances of each amino acid and glues them into a list and adds that back to the datagrame
    alanine_count_list=[]
    valine_count_list=[]
    arginine_count_list=[]
    histidine_count_list=[]
    lysine_count_list=[]
    asparticacid_count_list=[]
    glutamicacid_count_list=[]
    serine_count_list=[]
    threonine_count_list=[]
    asparagine_count_list=[]
    glutamine_count_list=[]
    cysteine_count_list=[]
    glycine_count_list=[]
    proline_count_list=[]
    isoleucine_count_list=[]
    leucine_count_list=[]
    methionine_count_list=[]
    phenylalanine_count_list=[]
    tyrosine_count_list=[]
    tryptophan_count_list=[]

    charge= {'A':1,'G':1,'L':1,'M':1,'F':1,'W':1,'K':1,'Q':1,'E':1,'S':1,'P':1,'V':1,'I':1,'C':1,'Y':1,'H':1,'N':1,'R':1,'D':1,'T':1}
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
     alanine_count=0
     valine_count=0
     arginine_count=0
     histidine_count=0
     lysine_count=0
     asparticacid_count=0
     glutamicacid_count=0
     serine_count=0
     threonine_count=0
     asparagine_count=0
     glutamine_count=0
     cysteine_count=0
     glycine_count=0
     proline_count=0
     isoleucine_count=0
     leucine_count=0
     methionine_count=0
     phenylalanine_count=0
     tyrosine_count=0
     tryptophan_count=0
     for b in range(len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])-peptide_skip): #b is the index of the amino acid in row 'a' peptide
      if optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='A':
          alanine_count= alanine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='G':
          glycine_count = glycine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='L':
          leucine_count = leucine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='M':
          methionine_count = methionine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='F':
          phenylalanine_count = phenylalanine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='W':
          tryptophan_count = tryptophan_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='K':
          lysine_count = lysine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='Q':
          glutamine_count = glutamine_count +1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='E':
          glutamicacid_count = glutamicacid_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='S':
          serine_count = serine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='P':
          proline_count = proline_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='V':
          valine_count = valine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='I':
          isoleucine_count= isoleucine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='C':
          cysteine_count = cysteine_count +1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='Y':
          tyrosine_count = tyrosine_count+1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='H':
          histidine_count = histidine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='N':
          asparagine_count = asparagine_count + 1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='R':
          arginine_count = arginine_count +1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='D':
          asparticacid_count = asparticacid_count+1
      elif optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=='T':
          threonine_count = threonine_count+1

     alanine_count_list.append(alanine_count)
     valine_count_list.append(valine_count)
     arginine_count_list.append(arginine_count)
     histidine_count_list.append(histidine_count)
     lysine_count_list.append(lysine_count)
     asparticacid_count_list.append(asparticacid_count)
     glutamicacid_count_list.append(glutamicacid_count)
     serine_count_list.append(serine_count)
     threonine_count_list.append(threonine_count)
     asparagine_count_list.append(asparagine_count)
     glutamine_count_list.append(glutamine_count)
     cysteine_count_list.append(cysteine_count)
     glycine_count_list.append(glycine_count)
     proline_count_list.append(proline_count)
     isoleucine_count_list.append(isoleucine_count)
     leucine_count_list.append(leucine_count)
     methionine_count_list.append(methionine_count)
     phenylalanine_count_list.append(phenylalanine_count)
     tyrosine_count_list.append(tyrosine_count)
     tryptophan_count_list.append(tryptophan_count)
    optimized_dataset['alanine_count_list']=alanine_count_list
    optimized_dataset['valine_count_list']=valine_count_list
    optimized_dataset['arginine_count_list']=arginine_count_list
    optimized_dataset['histidine_count_list']=histidine_count_list
    optimized_dataset['lysine_count_list']=lysine_count_list
    optimized_dataset['asparticacid_count_list']=asparticacid_count_list
    optimized_dataset['glutamicacid_count_list']=glutamicacid_count_list
    optimized_dataset['serine_count_list']=serine_count_list
    optimized_dataset['threonine_count_list']=threonine_count_list
    optimized_dataset['asparagine_count_list']=asparagine_count_list
    optimized_dataset['glutamine_count_list']=glutamine_count_list
    optimized_dataset['cysteine_count_list']=cysteine_count_list
    optimized_dataset['glycine_count_list']=glycine_count_list
    optimized_dataset['proline_count_list']=proline_count_list
    optimized_dataset['isoleucine_count_list']=isoleucine_count_list
    optimized_dataset['leucine_count_list']=leucine_count_list
    optimized_dataset['methionine_count_list']=methionine_count_list
    optimized_dataset['phenylalanine_count_list']=phenylalanine_count_list
    optimized_dataset['tyrosine_count_list']=tyrosine_count_list
    optimized_dataset['tryptophan_count_list']=tryptophan_count_list
    return optimized_dataset #amino acid acount 1

#https://onlinelibrary-wiley-com.login.ezproxy.library.ualberta.ca/doi/full/10.1002/prot.10560
def calculate_sidechain_energy(optimized_dataset):
    ENERGY= {'A':-1.35,'G':0,'L':-4.27,'M':-3.88,'F':-4.78,'W':-5.12,'K':-3.06,'Q':-2.12,'E':-2.03,'S':-0.72,'P':-3.01,'V':-3.52,'I':-4.48,'C':-1.78,'Y':-3.93,'H':-3.00,'N':-1.73,'R':-3.37,'D':-1.32,'T':-1.59}
    peptide_sidechain_energy=[]
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
        peptide_NRG=0
        b=0
        while b < len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')]): #b is the index of the amino acid in row 'a' peptide
            peptide_NRG= peptide_NRG + ENERGY[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]]
            b=b+1
        peptide_sidechain_energy.append(peptide_NRG)
    optimized_dataset['sidechain_energy']=peptide_sidechain_energy
    return optimized_dataset


#https://web.expasy.org/protscale/pscale/PolarityGrantham.html 
def calculate_sidechain_polarity(optimized_dataset):
    polarity= {'A':8.1,'G':9.0,'L':4.9,'M':5.7,'F':5.2,'W':5.4,'K':11.3,'Q':12.3,'E':12.3,'S':9.2,'P':8,'V':5.9,'I':5.2,'C':5.5,'Y':6.2,'H':10.4,'N':11.6,'R':10.5,'D':13.00,'T':8.6}  
    peptide_polarity=[]
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
        peptide_POL=0
        b=0
        while b < len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')]): #b is the index of the amino acid in row 'a' peptide
            peptide_POL= peptide_POL + polarity[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]]
            b=b+1
        peptide_polarity.append(peptide_POL)
    optimized_dataset['polarity']=peptide_polarity
    return optimized_dataset

#COHORTCHARGE is an artifactual feature that was found to significantly decrease accuracy because it only counted certain charged residues
#it was strongly correlated to binding because it used charge BUT it ignored much of the data so it baited the algorithm to pick it but was predictively redundant and missing key information
def calculate_cohort_charge(optimized_dataset): 
    charge_cohort_1=[]
    charge_cohort_2=[]
    charge_cohort_3=[]
    charge_cohort_4=[]
    charge= {'A':0,'G':0,'L':0,'M':0,'F':0,'W':0,'K':1,'Q':0,'E':-1,'S':0,'P':0,'V':0,'I':0,'C':0,'Y':0,'H':0,'N':0,'R':1,'D':-1,'T':0}
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
     peptide_charge1=0
     peptide_charge2=0
     peptide_charge3=0
     peptide_charge4=0
     b=0
     while b < len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')]): #b is the index of the amino acid in row 'a' peptide
      peptide_charge1= peptide_charge1 + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+0]] + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+1]]
      peptide_charge2= peptide_charge2 + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+1]] + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+2]]
      peptide_charge3= peptide_charge3 + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+2]] + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+3]]
      if b+4 < len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')]): #this is because there is no 20th index amino acid for the peptide so we have 1 less for the fourth cohort last pairing
       peptide_charge4= peptide_charge4 + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+3]] + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+4]]
      else: 
       peptide_charge4= peptide_charge4 + charge[optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b+3]]
      b=b+4

     charge_cohort_1.append(peptide_charge1)
     charge_cohort_2.append(peptide_charge2)
     charge_cohort_3.append(peptide_charge3)
     charge_cohort_4.append(peptide_charge4)
    optimized_dataset['charge_cohort_1']=charge_cohort_1
    optimized_dataset['charge_cohort_2']=charge_cohort_2
    optimized_dataset['charge_cohort_3']=charge_cohort_3
    optimized_dataset['charge_cohort_4']=charge_cohort_4
    return optimized_dataset #cohort charge 

def insert_secondary_structure_data(peptide_array_df): #uses the S4PRED output and sticks it onto the DF, may need adjustments if formatting differs signficantly

    print("in the order you loaded the arrays, you need to load the secondaries")
    array_count=int(input("how many peptide secondary structure arrays are you loading? : "))

    final_coil_prob_list_peptide=[]
    final_helix_prob_list_peptide=[]
    final_sheet_prob_list_peptide=[]

    for a in range(array_count):#traverse alla rrays and we append individual values at end
      print("Please enter the path of the ", a , "th peptide array  of interest: ")
      file_path=input(": ")
      frame_shift_array=input("please enter the frame shift of this peptide array--for example if the array shifts by 1 residue then enter 1 , etc ")
      file_path= file_path.replace('"',"")
      peptide_array_df=pd.read_excel(file_path,usecols='D,E,F',header=None) #THIS WILL PULL THE COLUMNS FOR THE SECONDARY STRUCTURE BASED ON HOW I FORMATTED THE FILES
     
      #excel_of_weights_model_10 = pd.ExcelWriter(r'C:\Users\sgrew\Desktop\applepie.xlsx')
      #peptide_array_df.to_excel(excel_of_weights_model_10)
      #excel_of_weights_model_10.save()


      #for each of the three provided probabilities, we pull a list and compress it to thge amount of peptides in the array
      coil_prob_list_residue=peptide_array_df[3].values.tolist()
      single_array_coil_prob_list_peptide= build_list_of_peptide_secondaries(coil_prob_list_residue,frame_shift_array)
    

      helix_prob_list_residue=peptide_array_df[4].values.tolist()
      single_array_helix_prob_list_peptide=build_list_of_peptide_secondaries(helix_prob_list_residue,frame_shift_array)

      sheet_prob_list_residue=peptide_array_df[5].values.tolist()
      single_array_sheet_prob_list_peptide=build_list_of_peptide_secondaries(sheet_prob_list_residue,frame_shift_array)

      final_coil_prob_list_peptide=final_coil_prob_list_peptide + single_array_coil_prob_list_peptide
      final_helix_prob_list_peptide=final_helix_prob_list_peptide + single_array_helix_prob_list_peptide
      final_sheet_prob_list_peptide=final_sheet_prob_list_peptide + single_array_sheet_prob_list_peptide


      #then we need to pull them into lists, pull 20 amino acids segments, average the secondarys structure values and then hop over by the frame shift and repeat

    return final_coil_prob_list_peptide,final_helix_prob_list_peptide,final_sheet_prob_list_peptide #secondary structure

def calculate_sulhpidepotential(optimized_dataset): # determines if there are 2 cysteins 3 amino acid or more away
    sulphide_list=[]
    c_index_list=[]
    sulphide_potential=0
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptides
     sulphide_potential=0
     c_index_list=[]

     #print("this is the peptide ",optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])

     for b in range(len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])): #b is the index of the amino acid in row 'a' peptide

      #print ("this is the amino acid ",optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b])

      if optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]=="C":
       c_index_list.append(b)
     #print("this is the c_index_list ", c_index_list)
     if len(c_index_list) > 0:
         for c in range(len(c_index_list)):
            if sulphide_potential==0:
             if abs(c_index_list[0]-c_index_list[c]) >= 3: #i changed this from 4 to 3 if this breaks 
              sulphide_potential=1
             else:
              sulphide_potential=0
     #print("this is the sulphide potential before appending ", sulphide_potential)
     sulphide_list.append(sulphide_potential)
     sulphide_potential=0
     c_index_list=[]
     #input("")
    optimized_dataset['sulphide_potential']=sulphide_list
    return optimized_dataset #sulphide poential 1

def initialize_features():
    features=[]
    features=create_analysis_options()
    features.append("coil_prob")
    features.append("helix_prob")
    features.append("sheet_prob")
    features.append('hydropathy')
    features.append('charge_At_7ph')
    #features.append('charge_cohort_1')
    #features.append('charge_cohort_2')
    #features.append('charge_cohort_3')
    #features.append('charge_cohort_4')
    features.append('positive_charge_At_7ph')
    features.append('negative_charge_At_7ph')
    features.append('alanine_count_list')
    features.append('valine_count_list')
    features.append('arginine_count_list')
    features.append('histidine_count_list')
    features.append('lysine_count_list')
    features.append('asparticacid_count_list')
    features.append('glutamicacid_count_list')
    features.append('serine_count_list')
    features.append('threonine_count_list')
    features.append('asparagine_count_list')
    features.append('glutamine_count_list')
    features.append('cysteine_count_list')
    features.append('glycine_count_list')
    features.append('proline_count_list')
    features.append('isoleucine_count_list')
    features.append('leucine_count_list')
    features.append('methionine_count_list')
    features.append('phenylalanine_count_list')
    features.append('tyrosine_count_list')
    features.append('tryptophan_count_list')
    features.append('sulphide_potential')
    features.append('sidechain_energy')
    features.append('polarity')

    print(features[0:len(features)])


    return features
####FUNCTIONS############################# SPECIFICALLY FOR MACHINE LEARNING AND RELATED FUNCTIONS

def set_data(dataset,features): ################################################################################################################ setting our features and label
    y=dataset["Raw Median"]
    X=dataset[features]
    return X,y

def machine_learning(X,y,frame_size): ################################################################################################################### single instance of model building

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)
    print(X_test)
    #input("a")
    #input("a")
    #input("a")
    #input("a")
    #input("a")
    #input("a")
    #input("a")
    #model = RandomForestRegressor(max_depth=frame_size+1,n_estimators=100,bootstrap=True) #boot strap falses includes whole data set
    model = DecisionTreeRegressor(max_depth=frame_size+1,criterion='squared_error')
    #print("start")
    #model  = MLPRegressor(hidden_layer_sizes=(100), max_iter=1000).fit(X_train, y_train)
    #print("done")
    #model = Ridge()
    #model = KNeighborsRegressor()
    #model = GradientBoostingRegressor()
    #model = IsotonicRegression()
    #model = evolution(model)
    #model = LinearRegression()
    #print(len(model.coefs_),"length total")
    #print(len(model.coefs_[0]),"length 0")
    #print(model.coefs_[0])
    ###############################################model = model.fit(X_train, y_train)
    model = model.fit(X,y)

    ##temporary
    peptide_array_data=load_pep() #this is commented out until I finish looping the pooping

    y_test=peptide_array_data["Raw Median"]
    features=initialize_features()
    X_test=peptide_array_data[features]
    ## temporary above


    X=peptide_array_data[features]
    predictions = model.predict(X_test)
    ##############################################predictions = model.predict(X_test)

    print(sklearn.tree.export_text(model))

    input("")
    #print(predictions)
    #input()
    print("mean absolute error is: ",mean_absolute_error(y_test, predictions))
    print("mean absolute squared error is: ",mean_squared_error(y_test, predictions))
    print("r2 is: ",r2_score(y_test, predictions))
    return mean_absolute_error(y_test, predictions),mean_squared_error(y_test, predictions),r2_score(y_test, predictions),model


def machine_learning_evaluation(X,y,frame_size): ######## broad machine learning function that builds the specified number of models and evaluates theri accuracy
    all_r2_values=[]
    all_absolute_error_means=[]
    all_squared_error_means=[]
    
    r2_sum=0
    mean_absolute_error_sum=0
    mean_squared_error_sum=0
    
    test_count=int(input("how many times would you like the model development to be performed?: "))
    
    smallest_error=10000000
    smallest_squared_error=10000000
    biggest_r2=-10000000
    iteration=0
    while_looper_variable=input("would you like to test multiple instances for hyperparameter optimization? O for yes, anything else for no")
    once_over_variable=1

    while int(while_looper_variable)==0 or once_over_variable==1:

     #while_looper_variable=input("loop? 0 for yes")

     #comment this out when done testing
     smallest_error=10000000
     smallest_squared_error=10000000
     biggest_r2=-10000000
     all_absolute_error_means=[]
     all_r2_values=[]
     all_squared_error_means=[]
     mean_absolute_error_sum=0
     mean_squared_error_sum=0
     r2_sum=0

     for z in range(test_count):

         #print("Model ", z+1, " of ", test_count)
         a,b,c,d=machine_learning(X,y,frame_size)
         #print(z)
         #this code should pull best model from a test set that is optimized
         if a<smallest_error:
             smallest_error=a
             smallest_error_squared_error=b
             minimum_error_model=d
         if b<smallest_squared_error:
             smallest_squared_error=b
             smallest_squared_error_error=a
             minimum_squared_error_model=d
         if c>biggest_r2:
             biggest_r2=c
             biggest_r2_error=a
             biggest_r2_squared_error=b
             max_r2_model=d
        
         #stats 
         all_absolute_error_means.append(a)

         all_r2_values.append(c)

         all_squared_error_means.append(b)

         mean_absolute_error_sum=mean_absolute_error_sum+a
         mean_squared_error_sum=mean_squared_error_sum+b
         r2_sum=r2_sum+c


     #Statistical analysis and evluation of models
     #print("the standard deviation of errors is : " ,(statistics.stdev(all_absolute_error_means)))
     #print("the standard deviation of errors is : " ,(statistics.stdev(all_squared_error_means)))
     #print("mean_squared_error_avg: ",mean_squared_error_sum/test_count)
     #print("mean_absolute_error_avg: ",mean_absolute_error_sum/test_count)
     #print("r2_avg: ",r2_sum/test_count)
     #print("the best of all r2 scores is ", biggest_r2)
     #print("the best of all mean absolute errors is ",smallest_error)
     #print("the best of all squared errors  is ", smallest_squared_error)

     #print(smallest_error,";",biggest_r2,";",smallest_squared_error,";",mean_squared_error_sum/test_count,";",r2_sum/test_count,";",mean_absolute_error_sum/test_count,";",(statistics.stdev(all_absolute_error_means)),";",(statistics.stdev(all_squared_error_means)))

     #iteration counter to count how many times we ran our test group and once over variable flip to make sure it only occurs once if such an option is selected
     iteration=iteration+1
     once_over_variable=0

     #for visualization if selected
     print(sklearn.tree.export_text(d))

     #this block is for different iteration tests
     if iteration==100:
         input("There has been " + iteration + " loops completed, this is the programmed pause point" )
     #if iteration==20:
     #    print("frame size has reached ",frame_size)
    
    all_best_features=""
    continue_variable='0'
    input("pause pause pause")

    #this block does  feature selection


    while continue_variable=='0':
     minimum_features=int(input("how many features this time?"))
     rfecv=RFE(estimator=minimum_squared_error_model,verbose=1,n_features_to_select=minimum_features)
     rfecv.fit(X,y)
     rfecv.transform(X)
     print(rfecv.ranking_)
     print('optimal number of features: {}'.format(rfecv.n_features_))
     print(np.where(rfecv.support_==True)[0])
     #print(np.where(rfecv.support_==False)[0])
     input()
     all_best_features= all_best_features + str(np.where(rfecv.support_==True)[0]) +"[][]spacer[][]"
     continue_variable=input("continue? 0 for yes, any other number for no")
    print(all_best_features)

    #This section is intended to be used for using the best model to test a specific test dataset for which there must be a list of peptides and sequences with the first row being labelled 'Peptides'
    # steph_test=load_pep_specific_test()
    #predictions = minimum_squared_error_model.predict(steph_test)
    #print(predictions)
    print("")
    print("Best?")


    #text_rep=tree.export_text(minimum_squared_error_model)
    #print(text_rep)
    #features= initialize_features()
    #fig=plt.figure(dpi=1200)
    #_= tree.plot_tree(minimum_squared_error_model,feature_names=features)
    #fig.savefig("decisiontree.png",bbox_inches="tight")

    input("done saving tree")
    input("")
    #while continue_variable=='0':
    # variance=int(input("select a variance threshhold"))
    # selector = VarianceThreshold(threshold=variance)
    # selector.fit_transform(X)
    # print(selector.get_support())
    # continue_variable=input("continue? 0 for yes, any other number for no")
    
    #feature_weights=minimum_squared_error_model.best_estimator_.feature_importances_.tolist() #this function pulls the best estimators feature weights 
    
    
    return minimum_error_model,minimum_squared_error_model,smallest_error,smallest_error_squared_error,smallest_squared_error_error,smallest_squared_error #a1:b2 is for storing stats about each model for naming purposes


#method for evolution-based model refinement that may be useful for designing an ideal model with a more complex algorithm or normally distributed dataset
def evolution(model): #an artifact function when experimenting with evolution as a way to improve models
    param_grid = {#'min_weight_fraction_leaf': Continuous(0.001, 0.5, distribution='log-uniform'),
                  #'max_depth': Integer(75,150),
                  'max_leaf_nodes': [500],
                  #'n_estimators':[50]
                  #'max_features': [1],
                  }
    
    cv = RepeatedKFold(n_splits=10, n_repeats=3)
    
    adapter_mutate = InverseAdapter(initial_value=0.9, end_value=0.1, adaptive_rate=0.1)
    adapter_crossover = InverseAdapter(initial_value=0.1, end_value=0.9, adaptive_rate=0.1)
    
    evolved_estimator = GridSearchCV(estimator=model,
                                     cv=cv,
                                     scoring='neg_mean_squared_error', #optimizes for mean squared, uses neg for maximization algo "negative_mean_squared_error"
                                     refit='neg_mean_squared_error',
                                     #population_size=10, 
                                     #generations=5, #iterations
                                     #tournament_size=20, 
                                     #elitism=True,
                                     #crossover_probability=adapter_crossover, #cross over between kept iterations
                                     #mutation_probability=adapter_mutate, #probability of unbiased mutations
                                     param_grid=param_grid, #hyper parameters
                                     #criteria='max',
                                     #algorithm='eaMuPlusLambda',
                                     n_jobs=-2, #gpu utilization -N where number of gpus used is (all+1-N)
                                     verbose=1,
                                     #keep_top_k=4,
                                     return_train_score=True,
                                   )

    return evolved_estimator

#model download function, requires pickle which may not work with all versions of pytorch or scikit learn
def download_model(minimum_error_model,minimum_squared_error_model,a1,a2,b1,b2):
    
    
    filename_add=input("what would you like to name the minumum error model? : ")
    filename1 = r'C:\Users\sgrew\Desktop\minimum_error_model-- '+  filename_add + '-ME-' + str(a1)+ '-MSE-' + str(a2)
    filename2 = r'C:\Users\sgrew\Desktop\minimum_squared_error_model-- '+  filename_add +'-ME-' + str(b1) + '-MSE-' + str(b2)
    with open(filename1, 'wb') as f1:
        pickle.dump(minimum_error_model, f1)  
    with open(filename2, 'wb') as f2:
        pickle.dump(minimum_squared_error_model,f2)  
    input('')    #download model
    
def download_feature_weights(aa_combos,feature_weights):
    ##################### building a data frame after cross val of best features and save excel
    
    print("starting feature weighting")
    feature_weight_dict={"Feature":aa_combos,"Weights":feature_weights}
    phase_1_weights=pd.DataFrame(feature_weight_dict)
    print(phase_1_weights)
    
    excel_of_weights_model_1 = pd.ExcelWriter(r'C:\Users\sgrew\Desktop\file dump for python\Feature_Weights_Model_1.xlsx')
    phase_1_weights.to_excel(excel_of_weights_model_1)
    excel_of_weights_model_1.save()
    
    print("done feature weighting")
    input("") #download feature weights
    
    ####################        



####FUNCTIONS############################# SPECIFICALLY FOR FUNCTIONS RELATING TO MODEL AND ACTIONS OF MODEL
#first functions allows you to select a model
def model_selection():
    filename=input("Please enter the filename of the model you are interested in: ")
    model=pickle.load(open(filename, 'rb'))
    return model

def prediction_model_on_array(model,variables,probe_array):
    predicted_outcomes=model.predict(probe_array[variables])
    print(predicted_outcomes)
    peptide_list = probe_array['Peptide'].tolist()  
    predicted_outcomes=predicted_outcomes.tolist() 
    df_tuple=list(zip(peptide_list,predicted_outcomes))
    df=pd.DataFrame(df_tuple,columns=['Peptide','Raw Median'])
    return df

####FUNCTIONS############################# SPECIFICALLY FOR ANALYZING 1 PEPTIDE BUILDING LIBRARY-DATAFRAME

def build_list_of_peptide_secondaries(coil_prob_list_residue,frame_shift_array):
    
    coil_prob_list_peptide=[]
    for a in range(0,len(coil_prob_list_residue)-19,int(frame_shift_array)):
        sum_of_20=0
        

        for b in range(0,20):

            sum_of_20=sum_of_20 + coil_prob_list_residue[a+b]

        coil_prob_list_peptide.append(sum_of_20)
    return coil_prob_list_peptide


def build_frame_library(peptide,frame_size): #this is inteded for individual peptides to build  a library of each frame_size length chain
    peptide_library=[]
    walking_frame=0
    while len(peptide_library)<(len(peptide)+1-frame_size-peptide_skip): #at the moment this only collect continuous chains but can be improved
        peptide_library.append(peptide[walking_frame+peptide_skip:walking_frame+frame_size+peptide_skip])
        walking_frame=walking_frame+1
    return peptide_library

def create_analysis_options(): #very simple collects all amino acids and pairs them iteratively
    all_amino_acids=['A','V','I','L','M','F','Y','W','S','T','N','Q','C','G','P','R','H','K','D','E']
    iterated_list=[]
    for a in all_amino_acids: #this loop creates a lsit of every pair
        for b in all_amino_acids:
            iterated_list.append(''.join([a,b]))   
    return iterated_list


#you can add more functions below this if you intend to add more features to model
def create_dataframe(iterated_list,library): #creates data frame for just each peptide in the library 
    peptide_dataframe=pd.DataFrame(library,columns=['Peptide'])
    peptide_dataframe['Raw Median']=0
    peptide_dataframe[iterated_list]=0
#under same function, it also modifies the data frame to have accurate info
    for peptide in library: #for each peptide in library
        pairs=[] #initializes list of aa pairs
        for index in range(len(peptide)-1): # for each amino acid in each peptide
            pairs.append(peptide[index]+peptide[index+1]) #stick together the pairs
        for pair in pairs:
            peptide_dataframe.iloc[library.index(peptide),iterated_list.index(pair)+peptide_dataframe.columns.get_loc('AA')]=peptide_dataframe.iloc[library.index(peptide),iterated_list.index(pair)+peptide_dataframe.columns.get_loc('AA')]+1 #takes relative indexes and increments them accordingly, starts at first pair in lsit
    return peptide_dataframe



def load_pep(): #this function will pul just peptide and their corresponding median value

    final_df= pd.DataFrame() #all the arrays will be appended tip to tail here
    print("Please enter the path of the  peptide array  of interest: ")
    file_path=input(":")
    file_path= file_path.replace('"',"")
    final_df=pd.read_excel(file_path)
    #peptide_array_df=pd.read_csv(file_path) #loads just the data in columns in list f
    return final_df #returns only the completed unedited array, the data needs to be incrememnted based upon their order and pairing


def load_pep_phase_1():

    final_df= pd.DataFrame() #all the arrays will be appended tip to tail here
    used_columns=['Peptide','Raw Median'] #this command should be an input later for ease of access and flexibility##########

    array_count=int(input("how many peptide arrays are you loading? : "))
    for a in range(array_count):
      print("Please enter the path of the ", a+1 , "th peptide array  of interest: ")
      file_path=input(":")
      file_path= file_path.replace('"',"")
      peptide_array_df=pd.read_excel(file_path,usecols=used_columns)
      final_df=pd.concat([final_df,peptide_array_df],ignore_index=True)

      #peptide_array_df=pd.read_csv(file_path,usecols=used_columns) #loads just the data in columns in list f
    
    return final_df #returns only the completed unedited array, the data needs to be incrememnted based upon their order and pairing

def load_pep_specific_test():

    final_df= pd.DataFrame() #all the arrays will be appended tip to tail here
    #used_columns=['Peptide'] #this command should be an input later for ease of access and flexibility##########

    array_count=int(input("how many peptide arrays are you loading? : "))
    for a in range(array_count):
      print("Please enter the path of the ", a+1 , "th peptide array  of interest: ")
      file_path=input(":")
      file_path= file_path.replace('"',"")
      peptide_array_df=pd.read_excel(file_path)
      final_df=pd.concat([final_df,peptide_array_df],ignore_index=True)
      print(peptide_array_df)

      #peptide_array_df=pd.read_csv(file_path,usecols=used_columns) #loads just the data in columns in list f
    
    return final_df #returns only the completed unedited array, the data needs to be incrememnted based upon their order and pairing

def make_list_of_secondary_structure_potential_by_ps4pred_array():

####FUNCTIONS############################# SPECIFICALLY FOR BUILDING A NEW DATAFRAME WITH POSITIONS and testing it
 return

def clean_up_data_array(optimized_dataset):
    #this function will only leave behind the HIGHEST value for each respective peptide
    optimized_dataset=optimized_dataset.sort_values('Raw Median', ascending=False).drop_duplicates('Peptide').sort_index()
    
    #this block should fix indeces by smashing the array and putting it back together - its pretty fast too as long as its done before all the features added
    raw_median_list_temp= optimized_dataset['Raw Median'].tolist()
    Peptide_list_temp=optimized_dataset['Peptide'].tolist()
    optimized_dataset=pd.DataFrame({'Peptide':Peptide_list_temp,'Raw Median':raw_median_list_temp})
    
    ###
    #this code will remove all values that are equal to 0 (consider shifting to a standard deviation threshholdor smnth 
    
   
    #for a in range(len(optimized_dataset)):
        #if optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Raw Median')]==0:
            #optimized_dataset.drop(str(a))
    ###
            
    return optimized_dataset
   

    

def build_new_array(optimized_dataset,frame_size):
    list_of_features=[]
    all_amino_acids=['A','V','I','L','M','F','Y','W','S','T','N','Q','C','G','P','R','H','K','D','E']
    for amino_acid in all_amino_acids:
        for a in range(frame_size):
            list_of_features.append(amino_acid+str(a+1))
    optimized_dataset[list_of_features]=0
    
    for a in range(len(optimized_dataset)): #a is index of row in the dataframe to get peptide
        for b in range(len(optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')])): #b is the index of the amino acid in row 'a' peptide
                optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('A1')+list_of_features.index(str((optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]))+str(b+1))] =optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('A1')+list_of_features.index(str((optimized_dataset.iloc[a,optimized_dataset.columns.get_loc('Peptide')][b]))+str(b+1))]+1 #pulls each letter forma  pepetide and attaches it to its respective position and then finds on the list of features what the index of that combination is and then increments that value by 1 in the correct locaiton. A explores each peptide and B explore every amino acid in each peptide
    return optimized_dataset,list_of_features

#### MAIN FUNCTION PHASE 1############################# SPECIFICALLY FOR CREATING A DATA FRAME FORM PEPTIDE ARRAY DATA AND INCLUDING ALL OF THE IMPORTANT FEATURES WE WANT

def phase_1_build_a_array():
    aa_combos=create_analysis_options()
    peptide_skip=int(input("How many amino acids would you like to skip? "))
    frame_size=int(input("How large is your desired frame? "))
    peptide_array_data=load_pep_phase_1()
    peptide_array_data[aa_combos]=0
    peptide_array_data['hydropathy']=0
    peptide_array_data['charge_At_7ph']=0
    peptide_array_data= calculate_cohort_charge(peptide_array_data)
    peptide_array_data= calculate_sulhpidepotential(peptide_array_data)
    peptide_array_data= calculate_charge_plus(peptide_array_data)
    peptide_array_data= calculate_charge_minus(peptide_array_data)
    peptide_array_data= calculate_sidechain_polarity(peptide_array_data)
    peptide_array_data= calculate_sidechain_energy(peptide_array_data)
    coil_list,helix_list,sheet_list=insert_secondary_structure_data(peptide_array_data)
    print(coil_list)
    input("")
    peptide_array_data["coil_prob"]=coil_list
    peptide_array_data["helix_prob"]=helix_list
    peptide_array_data["sheet_prob"]=sheet_list   
    peptide_array_data= calculate_hydropathy_peptide(peptide_array_data)
    peptide_array_data= calculate_charge(peptide_array_data)
    #this is where i would add more features for analysis and i have to append tot he feature list the added features. there are pep analysis softwares that can do alot of crazy stuff
    #pepfold2 is good for structures

    #this block loops to build all skip lengths
    while peptide_skip!=21:
     peptide_array_data= calculate_aminoacid_counts(peptide_array_data,peptide_skip)
     peptide_array_data=find_aa_pairs(peptide_array_data,aa_combos,peptide_skip)

     peptide_array_data=peptide_array_data.sort_values(by='Raw Median')
     peptide_array_data=peptide_array_data.drop_duplicates(subset="Peptide",keep='last')

     file_name=r'C:\Users\sgrew\Desktop\Whole_skip'
     file_name=file_name + str(peptide_skip) + '.xlsx'
     file_name = r"{}".format(file_name)
     print(file_name)

     excel_of_weights_model_10 = pd.ExcelWriter(file_name)
     peptide_array_data.to_excel(excel_of_weights_model_10)
     excel_of_weights_model_10.close()
     peptide_skip=peptide_skip+1
    
    return peptide_array_data

    
    
#### MAIN FUNCTION PHASE 2############################# SPECIFICALLY FOR CREATING THE INGREDIENT MODEL FROM ARRAY DATA SET TO BE USED IN PHASE 3

def phase_2_build_model_1():

    peptide_array_data=load_pep() #this is commented out until I finish looping the pooping

    y=peptide_array_data["Raw Median"]

    features=initialize_features()
    frame_size=int(input("How large is your desired frame? "))
    input("")



    X=peptide_array_data[features]

    minimum_error_model,minimum_squared_error_model,a1,a2,b1,b2=machine_learning_evaluation(X,y,frame_size) #a1:b2 is for storing stats about each model for naming purposes
    


    download_model(minimum_error_model,minimum_squared_error_model,a1,a2,b1,b2)
    

    
    return minimum_squared_error_model
    
            

#### MAIN FUNCTION PHASE 3############################# SPECIFICALLY FOR COMBINING PEPTIDE ANALYSIS AND MODEL_1 WORK TO FIND OPTIMAL FRAME

def phase_3_normalize_pca(model,whole_array_df):#this function should output a DF of raw medians paired to respective optimal frames-- later it will need model as input
    

    
    df_with_selected_frames=pd.DataFrame(columns=['Peptide','Raw Median'])
    
    list_of_best_frames=[] #after these loops run these two lists will be zippered together to build a new dataframe
    list_of_raw_medians=[]

    for a in range(len(whole_array_df)): #for every peptide in an array
        print(a, " of ", len(whole_array_df))
        
        
        
        library=build_frame_library(whole_array_df.iloc[a,whole_array_df.columns.get_loc('Peptide')],frame_size) #a is here for traversing list, could be more efficient with a while maybe
        individual_peptide_df=create_dataframe(aa_combos,library) #takes the library and turns it into a df
        prediction_individual_peptide_df=prediction_model_on_array(model,aa_combos,individual_peptide_df)#we took old df and just pulled peptides paired to Raw median
        #now we just need to pull only the highest raw median value from this list
        #these two variables are the best raw median and its corresponding peptide
        highest_raw_median=0
        best_frame=''
        
        #these two variables will be an ordered list of all highest raw medians sharing an index with corresponding peptides
        
        #this approach is to collect two lists and build a new DF of just the best ones
        #another potential approach is to have another column as an index that can be called on to signify the start peptide - a bit messier but maybe faster processing wise
        for b in range(len(prediction_individual_peptide_df)):
            print(b, " of ", len(prediction_individual_peptide_df))
            if prediction_individual_peptide_df.iloc[b,prediction_individual_peptide_df.columns.get_loc('Raw Median')]>highest_raw_median:
                highest_raw_median=prediction_individual_peptide_df.iloc[b,prediction_individual_peptide_df.columns.get_loc('Raw Median')]
                best_frame=prediction_individual_peptide_df.iloc[b,prediction_individual_peptide_df.columns.get_loc('Peptide')]
            list_of_best_frames.append(best_frame)
            #print(list_of_best_frames)
            print(whole_array_df.iloc[a,whole_array_df.columns.get_loc('Raw Median')])
        
            list_of_raw_medians.append(whole_array_df.iloc[a,whole_array_df.columns.get_loc('Raw Median')]) #a + raw median is the original raw median from array
            #print(list_of_raw_medians)
            
    optimized_dataset=pd.DataFrame({'Peptide':list_of_best_frames,'Raw Median':list_of_raw_medians})
    
    optimized_dataset =                                            clean_up_data_array(optimized_dataset)
    optimized_dataset,list_of_features=                            build_new_array(optimized_dataset,frame_size)    
    excel_of_output_data = pd.ExcelWriter(r'C:\Users\sgrew\Desktop\file dump for python\optimized_Dataframe_of_only_trimmedframes.xlsx')
    optimized_dataset.to_excel(excel_of_output_data)
    excel_of_output_data.save()
    
    return optimized_dataset, frame_size, list_of_features #this data set at the moment is a frame with the real raw median



#### MAIN FUNCTION PHASE 4 ############################# This function is unique to phase 4 because it involves having to skip on adding the secondary structure feature 
#this occurs when the probe set is a mishmash of peptides as opposed to a whole sequence
#in this event, the array is built as usualy except for secondary not being added , it must be added MANUALLY before loading it again

def phase_4_build_a_array_non_protein_probe():
    aa_combos=create_analysis_options()
    peptide_skip=int(input("How many amino acids would you like to skip? "))
    frame_size=int(input("How large is your desired frame? "))
    peptide_array_data=load_pep_phase_1()
    peptide_array_data[aa_combos]=0
    peptide_array_data['hydropathy']=0
    peptide_array_data['charge_At_7ph']=0
    peptide_array_data= calculate_cohort_charge(peptide_array_data)
    peptide_array_data= calculate_sulhpidepotential(peptide_array_data)
    peptide_array_data= calculate_charge_plus(peptide_array_data)
    peptide_array_data= calculate_aminoacid_counts(peptide_array_data,peptide_skip)
    peptide_array_data= calculate_charge_minus(peptide_array_data)
    peptide_array_data= calculate_sidechain_polarity(peptide_array_data)
    peptide_array_data= calculate_sidechain_energy(peptide_array_data)
    peptide_array_data=find_aa_pairs(peptide_array_data,aa_combos,peptide_skip)
    peptide_array_data= calculate_hydropathy_peptide(peptide_array_data)
    peptide_array_data= calculate_charge(peptide_array_data)
    #this is where i would add more features for analysis and i have to append tot he feature list the added features. there are pep analysis softwares that can do alot of crazy stuff
    #pepfold2 is good for structures

    probe_variable=input("If your probe set is a whole protein, enter 1, otherwise it is a mash of unrelated peptides and you will manually append secondary structures later and you can enter anything else")
    if int(probe_variable)==1:
        coil_list,helix_list,sheet_list=insert_secondary_structure_data(peptide_array_data)
        print(coil_list)
        input("")
        peptide_array_data["coil_prob"]=coil_list
        peptide_array_data["helix_prob"]=helix_list
        peptide_array_data["sheet_prob"]=sheet_list   
    
    excel_of_weights_model_10 = pd.ExcelWriter(r'C:\Users\sgrew\Desktop\probe_array.xlsx')
    peptide_array_data.to_excel(excel_of_weights_model_10)
    excel_of_weights_model_10.close()
    
    return peptide_array_data


###################################################################################MAIN######################################################333

def main():

 while True:

  



    start_point=input("Which step would you like to begin from (1,2,3,4) ? : ")
    print("1 will Create feature data set from scratch and then build models")
    print("2 will just build models from a dataset")
    print("3 is a principal component analysis unsupervised (this function contains its own set of intialize features, edit this to use certain features or add new)")
    print("4 is to load a model and test it on novel data")

    if start_point==    '1':
        peptide_array_data                            = phase_1_build_a_array()
        minimum_squared_error_model                   = phase_2_build_model_1()
        #optimized_dataset,frame_size,list_of_features = phase_3_probe_peptide_with_ingredient_model(minimum_squared_error_model,peptide_array_data)
        #a                                             = phase_4_build_recipe_model(optimized_dataset,frame_size,list_of_features)
        input("")
        
    elif start_point == '2':

        minimum_squared_error_model                   = phase_2_build_model_1()
        #optimized_dataset,frame_size,list_of_features = phase_3_probe_peptide_with_ingredient_model(minimum_squared_error_model,peptide_array_data)
        #phase_4_build_recipe_model(optimized_dataset,frame_size,list_of_features)
        input("")
        
    elif start_point == '3':  

        from sklearn.preprocessing import StandardScaler #for normalizing

        print("please input excel datafile contianing features and labwel")
        
        peptide_array_data=load_pep()#we have the whole data set here
        #print(peptide_array_data)
        del peptide_array_data['Peptide']  #we had to remove these as we compile features unsupervised so we dont want non feature data
        del peptide_array_data['index']
        del peptide_array_data['Raw Median']
        #print(peptide_array_data)
        peptide_array_data_normalized=StandardScaler().fit_transform(peptide_array_data)

        #we have a modifiable list of features here and the label itself so we can remake a df if need be
        
        columns=[]
        #columns.append('index')
        #columns.append('Peptide')
        #columns.append('Raw Median')
        aa_pairs=create_analysis_options()
        for a in range(len(aa_pairs)):
            columns.append(aa_pairs[a])
        columns.append('hydropathy')
        columns.append('charge_At_7ph')
        columns.append('charge_cohort_1')
        columns.append('charge_cohort_2')
        columns.append('charge_cohort_3')
        columns.append('charge_cohort_4')
        columns.append('sulphide_potential')
        columns.append('positive_charge_At_7ph')
        columns.append('alanine_count_list')
        columns.append('valine_count_list')
        columns.append('arginine_count_list')
        columns.append('histidine_count_list')
        columns.append('lysine_count_list')
        columns.append('asparticacid_count_list')
        columns.append('glutamicacid_count_list')
        columns.append('serine_count_list')
        columns.append('threonine_count_list')
        columns.append('asparagine_count_list')
        columns.append('glutamine_count_list')
        columns.append('cysteine_count_list')
        columns.append('glycine_count_list')
        columns.append('proline_count_list')
        columns.append('isoleucine_count_list')
        columns.append('leucine_count_list')
        columns.append('methionine_count_list')
        columns.append('phenylalanine_count_list')
        columns.append('tyrosine_count_list')
        columns.append('tryptophan_count_list')
        columns.append('negative_charge_At_7ph')
        columns.append("coil_prob")
        columns.append("helix_prob")
        columns.append("sheet_prob")
        
        #print(peptide_array_data.shape)
        peptide_array_data_normalized = pd.DataFrame(peptide_array_data_normalized,columns=columns) #normalized array without peptide labels or the index colummn
        
        from sklearn.decomposition import PCA
        print(peptide_array_data_normalized)
        pca_components=PCA(n_components=2)
        pca_cool=pca_components.fit_transform(peptide_array_data_normalized)
        #print(pca_components.get_feature_names_out)
        #print(pca_components.n_features_)
        ram_ranch=pca_components.inverse_transform(pca_cool)
        #print(pca_components.score(peptide_array_data))
        #print(pca_components.get_feature_names_out)
        #print(pca_components.get_feature_names_out)
        #print(pca_components.get_feature_names_out)
        print(pd.DataFrame(data =ram_ranch,columns=columns))
        pca_cool=pd.DataFrame(data =pca_cool, columns = ['principal component 1', 'principal component 2'])
        excel_file_temp = pd.ExcelWriter(r'C:\Users\sgrew\Desktop\prinicpalcomponentrawmedian.xlsx')

        pca_cool.to_excel(excel_file_temp)
        excel_file_temp.close()

        print(pca_cool)

        input("")
        
    elif start_point == '4':

        build_or_download_variable=input("If you need to build a new array press 1, if you would like to load an existing one press anything else")
        if int(build_or_download_variable)==1:
            probe_array=phase_4_build_a_array_non_protein_probe() #build array with features from peptides AND saves it. if you need to append secondaries, do that after this step
        else:
            probe_array=load_pep_specific_test() #loads an array to test if you need to

        test_model=model_selection() #load a pre existing model

        variables=initialize_features() #gets active feature list

        results_df=prediction_model_on_array(test_model,variables,probe_array)
        
        excel_file_temp = pd.ExcelWriter(r'C:\Users\sgrew\Desktop\predicted_values.xlsx')

        results_df.to_excel(excel_file_temp)
        excel_file_temp.close()        

    else:
        print("this input was not recognized -- try again")


main()