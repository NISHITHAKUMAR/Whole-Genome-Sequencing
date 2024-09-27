import os
from Bio import SeqIO
import requests
import pandas as pd 
import numpy as np 
from  important_kmers import important_Kmer
from classification_kmers import Classification_metrics
from collections import Counter
# import warnings filter
from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category= FutureWarning)
simplefilter(action='ignore', category= UserWarning)
simplefilter(action='ignore', category= DeprecationWarning)
import warnings
warnings.filterwarnings('ignore')

# Dictionary Antibiotics Abbrivation
antibiotics = {
    'IMIPENEM': 'IPM',
    'GENTAMICIN': 'GEN',
    'NALIDIXIC_ACID': 'NAL',
    'AMPICILLIN': 'AMP',
    'CEFOPERAZONE': 'CZO',
    'CEFUROXIME': 'CXM',
    'MEROPENEM': 'MEM',
    'CEFEPIME': 'CZM',
    'AMIKACIN': 'AMK',
    'TRIMETHOPRIM': 'TMP',
    'TIGECYCLINE': 'TGC',
    'COLISTIN': 'COL',
    'AMOXICILLIN': 'AMC',
    'CEFTRIAXONE': 'CRO',
    'CIPROFLOXACIN': 'CIP',
    'PIPERACILLIN': 'PIP'
}

# print(antibiotics.keys())
# print(antibiotics.get('AMIKACIN'))


## Functions ......................................................
# creating  uodate locus tag from  sample1
def updatetag(genome_path):
    if not os.path.isfile(genome_path):
        print("Updating all genome sequence header..........")
        for root,dir, files in os.walk(genome_path):
            for file in files:
                if not file.endswith("_locus.txt"):
                    if file.endswith(".txt"):
                        # print(root)
                        txt_name  =  file.replace(".txt","")
                        new_txt = file.replace(".txt","_locus.txt")
                        txt_files = os.path.join(root, new_txt)
                        with open(txt_files,"w") as new_txt:
                            for seq_record in SeqIO.parse(os.path.join(root, file), "fasta"):
                                        desid = seq_record.description.find("gene=") +5
                                        if (desid - 5) == -1 :
                                            locus = f""">NAN_{txt_name.replace("_CDS","").split("_")[-1]}_{len(str(seq_record.seq))}\n"""
                                            new_txt.write(locus)
                                            new_txt.write(str(seq_record.seq)+"\n")                     
                                        else:
                                            locus = f""">{seq_record.description[desid:].split("]")[0]}_{txt_name.replace("_CDS","").split("_")[-1]}_{len(str(seq_record.seq))}\n"""
                                            new_txt.write(locus)
                                            new_txt.write(str(seq_record.seq)+"\n")                     

# download files .................................................
def download_file(url, destination):
    # print(url)
    response = requests.get(url)
    print(response)
    if response.status_code == 200:
        with open(destination, 'wb') as file:
            file.write(response.content)
        print(f"File downloaded successfully to {destination}")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

# Input fasta file
def CRT_Fasta_With_CDHIT(input_cdhit_res):
    if not os.path.isdir("AllGenome-fasta"):
         os.mkdir("AllGenome-fasta")
    os.chdir("AllGenome-fasta")     
    # Read the input fasta file
    records = SeqIO.parse(input_cdhit_res, "fasta")
    # Dictionary to store sequences based on the group
    sequence_dict = {}
    # Iterate through each sequence and group them based on "NEC1", "NEC2", etc.
    for record in records:
        # Extract the group identifier (e.g., NEC1, NEC2) from the sequence header
        group_identifier = record.id.split('_')[-2]
        # Add the record to the corresponding group in the dictionary
        if group_identifier not in sequence_dict:
            sequence_dict[group_identifier] = []
        sequence_dict[group_identifier].append(record)
    # Write sequences to separate fasta files for each group
    for group_identifier, group_records in sequence_dict.items():
        output_file = f"{group_identifier}.fasta"
        SeqIO.write(group_records, output_file, "fasta")
    print("Fasta files for each genome using teh cds-hit result  created successfully.")
    os.chdir("..")

def Genome_Tester(Genome_test_src,AllGenome_fasta_path):
    os.makedirs("Results_List", exist_ok=True)
    os.chdir("Results_List")

    # Loop through each fasta file in the directory
    for fasta_file in os.listdir(AllGenome_fasta_path):
        if fasta_file.endswith(".fasta"):
            fasta_path = os.path.join(AllGenome_fasta_path, fasta_file)
            print(fasta_path)
            # Generate output file names based on the fasta file name
            output_prefix = os.path.splitext(fasta_file)[0]
            output_list_file = f"{output_prefix}_{k_value}.list"
            output_txt_file = f"{output_prefix}.txt"

            # commands
            os.system(f"glistmaker {fasta_path} -w 13 -o {output_prefix}")
            
            
    list_dir = os.getcwd()

    os.chdir("..")
    os.makedirs("Results_Kmer", exist_ok=True)
    os.chdir("Results_Kmer")
    for list_file in os.listdir(list_dir):
        if list_file.endswith(".list"):
            list_path = os.path.join(list_dir, list_file)
            output_prefix = os.path.splitext(list_file)[0]
            output_txt_file = f"{output_prefix}.txt"
            # os.system(f"{Genome_test_src}/glistquery {list_path} > {output_txt_file}")
            os.system(f"glistquery {list_path} > {output_txt_file}")

    os.chdir("..")

def kmer_res_proc(Results_Kmer_path):
    # os.chdir()
    # print(os.path.basename())
    if not os.path.isdir("Processed_Kmer"):
        os.mkdir("Processed_Kmer")
    all_seq = []
    all_kmer = []
    genome_names = []
    for file in sorted(os.listdir(Results_Kmer_path)):
        if file.endswith(".txt"):
            txt_file = os.path.join(Results_Kmer_path, file)
            output_prefix = file.split("_")[0:2]
            gemone_name = '_'.join(output_prefix)
            print(gemone_name)
            kmer_data = pd.read_csv(txt_file, sep="\t")
            seq_col = kmer_data.iloc[:, 0]
            kmer_col = kmer_data.iloc[:, 1]
            
            seq_col_list = seq_col.tolist()
            kmer_col_list = kmer_col.tolist()
                    
            seq_col_list.insert(0, gemone_name.replace("_13.txt",""))
            kmer_col_list.insert(0, gemone_name.replace("_13.txt",""))
            
            seq_col = pd.Series(seq_col_list)
            kmer_col = pd.Series(kmer_col_list)
            
            all_seq.append(seq_col)
            all_kmer.append(kmer_col)
            genome_names.append(gemone_name)
            
            
    seq_dataframe = pd.concat(all_seq, axis = 1)

    kmer_dataframe  = pd.concat(all_kmer, axis = 1)
    kmer_dataframe = kmer_dataframe.fillna(0)
    # Transposing the dataframes
    seq_dataframe = seq_dataframe.T
    kmer_dataframe = kmer_dataframe.T
    seq_dataframe.sort_values(by = [seq_dataframe.columns[0]],inplace=True)
    kmer_dataframe.sort_values(by = [kmer_dataframe.columns[0]],inplace=True)

    # print(seq_dataframe.head())
    # print(kmer_dataframe.head())
    print("Data Writing... ")
    seq_dataframe.to_csv("Processed_Kmer/Kmer_Seq_h.txt", sep="\t", index=False)
    kmer_dataframe.to_csv("Processed_Kmer/Kmer_Data_h.txt", sep="\t", index=False, header=False)
    seq_dataframe.iloc[:, 1:].to_csv("Processed_Kmer/Kmer_Seq.txt", sep="\t", index=False)
    kmer_dataframe.iloc[:, 1:].to_csv("Processed_Kmer/Kmer_Data.txt", sep="\t", index=False, header=False)

    print("Data writting completed... ")

def Classification_Result_Excel(clf_results_folder, output_excel_path):
    res_data = []
    for file in os.listdir(clf_results_folder):
        if "result" in file:
            file_path = os.path.join(clf_results_folder, file)
            df = pd.read_csv(file_path, header = [0])
            ant_name = os.path.splitext(file)[0].split("__")[1]
            # print(ant_name+".......................")
            df.insert(0, "ANTIBIOTIC", ant_name)
            df.rename(columns={'Unnamed: 0': "Metrics"}, inplace=True)
            df.reset_index(inplace=True)
            df = df.iloc[:,1:]
            new = pd.pivot_table(df, index='ANTIBIOTIC', columns='Metrics')
            # print(new)
            res_data.append(new)
            
    data_all  = pd.concat(res_data, axis=0)
    data_all.reset_index(inplace=True)

    ant_col = data_all.iloc[:,0]
    data_new = data_all
    df = data_new
    df.columns = ['_'.join(col) for col in df.columns.values]
    columns = df.columns
    for i in range(1, len(columns), 2):
        col1 = columns[i]
        model_name =col1.split('.')[0]
        col2 = columns[i + 1] if i + 1 < len(columns) else None
        # Combine columns and create a new column
        if col2:
            combined_col_name = f'''{model_name}_{col1.replace(model_name,"")}_{col2.split("_",1)[1]}'''
            df[combined_col_name] = df[col1].astype(str) + "+-"+df[col2].astype(str)
    df.drop(columns=columns[1:],inplace=True)
    # Split columns by "_" and create a MultiIndex
    multiindex = pd.MultiIndex.from_tuples([col.split("_", 1) for col in df.columns], names=["First_Level", "Second_Level"])
    df.columns = multiindex
    # print(df)
    grouped = df.groupby(axis=1, level=0)

    algos = {}
    for name, group in grouped:
        if not name =="ANTIBIOTIC":
            acc_mean = np.mean(group.iloc[:, 1].str.split("+",expand=True)[0].astype(float))
            algos[name] = acc_mean
    max_key = max(algos, key=lambda k: algos[k])
    max_df = grouped.get_group(max_key)
    max_df.insert(0, "ANTIBIOTIC", grouped.get_group("ANTIBIOTIC")["ANTIBIOTIC"])

    with pd.ExcelWriter(output_excel_path) as wr:
        max_df.to_excel(wr, sheet_name = "ML Metrics-"+max_key)
        df.to_excel(wr, sheet_name= "ML Metrics - Classifiers")
   
def seq_select(Imp_Kmers_Feature_path,kmer_seq_data_path):
    if not os.path.isdir("Selected_Kmer_seq"):
        os.mkdir("Selected_Kmer_seq")
    kmer_seq_data = pd.read_csv(kmer_seq_data_path+"/Kmer_Seq_h.txt",sep="\t")
    with open("Selected_Kmer_seq/selected_seq.txt","w") as table:
        table.write(f"""gname\tindex\tantibiotic\tsequence\n""")
        for file in os.listdir(Imp_Kmers_Feature_path):
            if file.endswith(".txt"):
                txt_file = os.path.join(Imp_Kmers_Feature_path, file)
                data = [0]
                with open(txt_file, 'r') as file1:
                    for line in file1.readlines():
                        # we added genome column at first position thats why increases index by one 
                        #to select right sequences
                        data.append(int(line.strip())+1)
                # print(data)
                select_col = kmer_seq_data.iloc[:, data]
                ant_name = file.split("__")[-1].split(".")[0]
                ant_name = ant_name.replace(" ", "_")
                with open("Selected_Kmer_seq/"+ant_name+"_selected.fasta","w")as fasta: 
                    for gname in select_col[select_col.columns[0]]:
                        data1 = select_col[ select_col[select_col.columns[0]] == gname]
                        data1 = data1.iloc[:,1:]
                        for i in range(len(data1.columns)):
                            if type(data1.iloc[0,i]) == str:
                                fasta.write(f">{gname}_{data1.columns[i]}_{ant_name}\n")
                                fasta.write(data1.iloc[0,i]+"\n")
                                # -1 index ..........
                                table.write(f"""{gname}\t{int(data1.columns[i])-1}\t{ant_name}\t{data1.iloc[0,i]}\n""")
    
    # os.system(f"""cat Selected_Kmer_seq/*.fasta > Selected_Kmer_seq/All_ImpKmer_seq.fasta""")                



def createKMERExcel(blast_res_dir_path,output_Excel1,ReferenceGeneCatalog,antibiotics_acr,input_amr_path,Selected_Kmer_seq_path):
    data_list = []
    sheet_names = []
    gene_list = []
    gen_dic = {}
    ref_df = pd.read_csv(ReferenceGeneCatalog,sep="\t")
    ref_df.rename(columns = {"gene_family": "subject_acc.ver"}, inplace = True) 
    for file in os.listdir(blast_res_dir_path):
        if file.endswith(".txt"):
            blast_res_path  = os.path.join(blast_res_dir_path, file)
            data_frame  = pd.read_csv(blast_res_path, sep = "\t")
            data_frame.columns = ["query_acc.ver", "subject_acc.ver", "identity", "alignment_length", "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end","evalue", "bit_score"]
            data_frame["kmer_length"] = 13
            data_frame["actual_match_length"] = (data_frame["alignment_length"] -data_frame["mismatches"])-data_frame["gap_opens"] 
            data_frame["Query_coverage"] = (data_frame["actual_match_length"]/data_frame["kmer_length"])*100
            data_frame = data_frame[data_frame["Query_coverage"] >= 90]
            data_frame["subject_acc.ver"] = data_frame["subject_acc.ver"].str.split("|",expand=True)[6]
            data_frame.drop_duplicates(subset=["query_acc.ver","subject_acc.ver"],keep="first",inplace=True)       

            # print(data_frame)
            data_frame = pd.merge(data_frame,ref_df,how="left",on="subject_acc.ver")
            data_frame.drop_duplicates(subset=["query_acc.ver","subject_acc.ver"],keep="first",inplace=True)            
            sheet_name = file.split("_")[0]
            data_list.append(data_frame)
            sheet_names.append(sheet_name)

            gen_dic[sheet_name] = Counter(data_frame["subject_acc.ver"])
            gene_list.append(gen_dic)
    # Create a Pandas Excel writer using XlsxWriter as the engine
    data = {}
    for nested_dict in gene_list:
        for col_name, col_dict in nested_dict.items():
            for idx, value in col_dict.items():
                if col_name not in data:
                    data[col_name] = {}
                data[col_name][idx] = value
    # Create DataFrame from new dictionary
    df = pd.DataFrame.from_dict(data, orient='index')
    # Transpose DataFrame
    df = df.transpose()
    # Reset index
    df.reset_index(inplace=True)      
    df.fillna(0,inplace=True)
    # df.rename(columns=antibiotics_acr,inplace=True)
    print(df.head())
    try:
        df.rename(columns={"index":"subject_acc.ver"},inplace=True)
    except:
        print("trouble in index ")    

    df.set_index('subject_acc.ver',inplace=True)
    df[df > 1] = 1
    df.reset_index()        
    # print(df) 

    newdf = pd.concat(data_list,axis=0)
    newdf.drop(columns=["query_acc.ver","identity","alignment_length","mismatches","gap_opens","q.start",
        "q.end","s.start","s.end","evalue","bit_score","kmer_length","actual_match_length"],inplace=True)
    newdf.drop_duplicates(["subject_acc.ver"],keep="first",inplace=True)
    newdf = pd.merge(newdf,df,how="left",on=["subject_acc.ver"])


    
    # Sav   
    writer = pd.ExcelWriter(output_Excel1, engine='xlsxwriter')
    newdf.to_excel(writer, sheet_name="Selected Genes", index=False)#,header = None)



    AMR_data = pd.read_csv(input_amr_path)
    print(AMR_data.columns)
    # Write each dataframe to a different worksheet
    print(sheet_name)
    for i, df in enumerate(data_list):
        data = df.iloc[:,:2]
        unique_gene = data["subject_acc.ver"].value_counts().reset_index()
        unique_gene.columns = ["subject_acc.ver","Count"]
        unique_gene = pd.merge(data,unique_gene,on="subject_acc.ver",how="left")
        # input_amr_path = "/home/Desktop/Nishitha_MacBook/AMR_AI_MAIN/Inputs/Demo_AMR_data_RSI.csv"
        total_genome = len(AMR_data["GENOME_ID"])
        print("Total number of genomes in the dataset: ", total_genome)
        RS_data_list  = []
        for query_acc in data["query_acc.ver"]:
            index = int(query_acc.split("_")[1])
            antibio = "_".join(query_acc.split("_")[2:])
            with open(Selected_Kmer_seq_path+"/"+antibio+"_selected.fasta") as kmer_seq:
                seq = []
                for  line in kmer_seq:
                    if ">" in line:
                        if int(line.split("_")[1]) == index:
                            seq.append(line.split("_")[0].replace(">",""))
          
            RS_count  = AMR_data[AMR_data["GENOME_ID"].isin(seq)][antibio].to_list()    
            for i in range(len(RS_count)):
                if not  RS_count[i] == "R":
                    RS_count[i] = "S"
            RS_data_list.append([query_acc,total_genome,round(RS_count.count("R"),4),round(RS_count.count("S"),4),round(((RS_count.count("R"))/total_genome)*100,4),round((RS_count.count("S")/total_genome)*100,4),round(((RS_count.count("S")+RS_count.count("R"))/total_genome)*100,4)])
        print(antibio)
        RI_data = pd.DataFrame(data=RS_data_list, columns=['query_acc.ver', 'Total_no_Sample',"Count_R","Count_S","Percentage_of_R","Percentage_of_S","Percentage_of_RS"])
        unique_gene = pd.merge(unique_gene,RI_data,on=["query_acc.ver"],how="left")
        unique_gene.drop_duplicates(subset =["query_acc.ver","subject_acc.ver"])
        print(df.shape,"\n")
        df = pd.merge(df,unique_gene,on=["query_acc.ver","subject_acc.ver"], how="left")
        print(df.shape,"\n")
        df.drop_duplicates(keep ="first", inplace=True)
        print(df.shape,"\n")
        # try:
        #    sheetName = antibiotics_acr.get(sheet_names[i])
        #    print(sheetName)
        # except:
        #     sheetName = sheet_names[i]
        #     print(sheetName,"except")     
               
        df.to_excel(writer, sheet_name="_".join(df["query_acc.ver"][0].split("_")[2:]), index=False)#,header = None)
    # Save the Excel file
    writer._save()



def createCDSxls(blast_res_CDS_path,ReferenceGeneCatalog,input_genome_xls,clstr_file,output):
    input_genome_xls = pd.read_excel(input_genome_xls)
    #.............inout genome xlsx
    ref_df = pd.read_csv(ReferenceGeneCatalog,sep="\t")
    ref_df.rename(columns = {"gene_family": "subject_acc.ver"}, inplace = True) 
    data_frame  = pd.read_csv(blast_res_CDS_path , sep = "\t",comment="#")
    data_frame.columns = ["query_acc.ver", "subject_acc.ver", "identity", "alignment_length", "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end","evalue", "bit_score"]

    data_frame["kmer_length"] =  data_frame["query_acc.ver"].str.split("_",expand=True)[2]
    data_frame["kmer_length"] =  data_frame["kmer_length"].astype(int)
    data_frame["actual_match_length"] = (data_frame["alignment_length"] -data_frame["mismatches"])-data_frame["gap_opens"] 
    data_frame["Query_coverage"] = (data_frame["actual_match_length"]/data_frame["kmer_length"])*100
    data_frame = data_frame[data_frame["Query_coverage"] >= 90]
    data_frame["subject_acc.ver"] = data_frame["subject_acc.ver"].str.split("|",expand=True)[6]
    # print(data_frame)
    data_frame = pd.merge(data_frame,ref_df,how="left",on="subject_acc.ver")
    data_frame.sort_values(["subject_acc.ver"],ascending=False,inplace=True)
    data_frame.drop_duplicates(subset=["query_acc.ver","subject_acc.ver"],keep="first",inplace=True)            
    data =[]
    for _, line in enumerate(open(clstr_file,'r')):
        if not line.strip().startswith(">"):
            indx = line.find(">")
            gene = line[indx+1:].split("_")[0].strip()
            gnome = line[indx+1:].split("_")[1]
            data.append([gene,gnome])
    data = pd.DataFrame(data=data,columns=["first","Accession no & Genome Name"])
    data = data[data["first"] != "NAN"]    
    data = data[data["first"].isin(data_frame["subject_acc.ver"].str.strip().to_list())]
    data = data.pivot_table(index='Accession no & Genome Name', columns='first', aggfunc='size', fill_value=0)
    data[data > 1] = 1
    data = data.reset_index()
    input_genome_xls = pd.merge(input_genome_xls,data,how="left",on=["Accession no & Genome Name"])


    # Create a Pandas Excel writer using XlsxWriter as the engine
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    # Write each dataframe to a different worksheet
    input_genome_xls.to_excel(writer, sheet_name="Summary-Results", index=False)#,header = None)
    data_frame.to_excel(writer, sheet_name="AMR-CDS-Results", index=False)#,header = None)
    writer._save()




# code Execution.............................
print("code execution started.............") 

# 1 create All_cds.fasta 
input_genome_xls = os.getcwd()+"/Inputs/MIC values of NEC Genomes.xlsx"
genome_path = "Genome"
All_cds = "AMR-Results/Genome/All_cds.fasta"
to_download = ['https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/data/latest/AMR_CDS',"https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/data/latest/ReferenceGeneCatalog.txt"]
thr = 3
memory = 250000
# Set the path to the directory containing your fasta files
fasta_file = "cd-hit-results"
Genome_test_src = ""
# Genome_test_src = "home/Nishitha_MacBook/NGS_tools/GenomeTester4-master/src""

# Set the parameters for glistmaker
k_value = 13
Root_dir = os.getcwd()



# creating AMR folder first ..........

if not os.path.isdir("AMR-Results"):
     print("creating AMR-Results folder ................")
     os.mkdir("AMR-Results")
# updating all the tags from Genome folders files
updatetag(genome_path)        



# STEP.1 create  Genome/All_cds.fasta  if not presents ............
if not os.path.isdir("AMR-Results/Genome"):
     os.mkdir("AMR-Results/Genome")

if not os.path.isfile(All_cds):
        print("creating ............. AMR-Results/Genome/All_cds.fasta")
        txt_files = ""
        for root,dir, files in os.walk(genome_path):
            for file in files:
                # print(files)
                if file.endswith("_locus.txt"):
                        txt_files += f"{os.path.join(root, file)} "
        os.system(f"""cat {txt_files} >{All_cds}""")
else:
    print("File Already Exists.......AMR-Results/Genome/All_cds.fasta\n")


#STEP.2 creating cd-hit-results from All_cds.fasta 
if not os.path.isdir("AMR-Results/cd-hit-results"):
    os.mkdir("AMR-Results/cd-hit-results")
    os.system(f"""cd-hit -i {All_cds} -o AMR-Results/cd-hit-results/cd-hit-results -T {thr} -M {memory}""")
if not os.path.isfile("AMR-Results/cd-hit-results/cd-hit-results"):    
       print("creating ............. cd-hit-results/cd-hit-results")
       os.system(f"""cd-hit -i  {All_cds} -o AMR-Results/cd-hit-results/cd-hit-results -T {thr} -M {memory}""")
else:
    print("File already created in  AMR-Results/cd-hit-results/cd-hit-results")

input_cdhit_res = os.path.join(Root_dir, "AMR-Results/cd-hit-results/cd-hit-results")

#STEP.3 downloading  AMR_CDS and  ReferenceGeneCatalog.txt
if not os.path.isdir("AMR_gene_SEQ"):
    os.mkdir("AMR_gene_SEQ")
for link in to_download:
    if os.path.basename(link) == "ReferenceGeneCatalog.txt":
        destination_file = Root_dir +"/AMR_gene_SEQ/"+os.path.basename(link)
    else:
        destination_file = Root_dir +"/AMR_gene_SEQ/"+os.path.basename(link)+".fasta"
    if not  os.path.isfile(destination_file):
        download_file(link, destination_file)
        print(destination_file)


# STEP.4 
#.............................................................................................................
os.chdir("AMR-Results")
print("program directory : ",os.getcwd())
if not os.path.isdir("Kmer-find"):
     os.mkdir("Kmer-find")
os.chdir("Kmer-find")  
print("Accesing Kmer-find directories...............")

CRT_Fasta_With_CDHIT(input_cdhit_res = input_cdhit_res)
print("CRT_Fasta_With_CDHIT run successfully")
Genome_Tester(Genome_test_src,AllGenome_fasta_path = os.getcwd()+"/AllGenome-fasta")

kmer_res_proc(Results_Kmer_path = os.getcwd()+"/Results_Kmer")

os.chdir("..")
#.............................................................................................................
if not os.path.isdir("ML_Results"):
     os.mkdir("ML_Results")
os.chdir("ML_Results")  

important_Kmer(input_amr_data = Root_dir+"/Inputs/Demo_AMR_data_RSI.csv" ,kmer_data =  Root_dir+"/AMR-Results/Kmer-find/Processed_Kmer/Kmer_Data.txt")


Imp_Kmers_Feature_path = Root_dir+"/AMR-Results/ML_Results/results/Imp_Kmers_Feature"
kmer_seq_data_path = Root_dir+"/AMR-Results/Kmer-find/Processed_Kmer"
print("..........................")
seq_select(Imp_Kmers_Feature_path,kmer_seq_data_path)

blast_seq_path = Root_dir+"/AMR-Results/ML_Results/Selected_Kmer_seq"

# Classification Results Metrics
#STEP.2 creating cd-hit-results from All_cds.fasta
if not os.path.isdir("results/Classification_Results"):
    os.mkdir("results/Classification_Results")
clf_results_folder = Root_dir+"/AMR-Results/ML_Results/results/Classification_Results"
Antibiotic_data_path = Root_dir+"/Inputs/Demo_AMR_data_RSI.csv"
kmer_path = Root_dir+"/AMR-Results/Kmer-find/Processed_Kmer/Kmer_Data.txt"
feature_dir_path = Root_dir+"/AMR-Results/ML_Results/results/Imp_Kmers_Feature"
Classification_metrics(clf_results_folder,Antibiotic_data_path,kmer_path,feature_dir_path)


#create classification excel....//...here

os.chdir("..")  

# create database ...........
AMR_CDS_file = Root_dir + "/AMR_gene_SEQ/AMR_CDS.fasta"    
# AMR_CDS_file = Root_dir + "/CARD_Gene_Seq/nucleotide_fasta_protein_homolog_model.fasta"    

if not os.path.isdir("AMR_Blast-Results-CDS"):
    os.mkdir("AMR_Blast-Results-CDS")
print("Creating Blast Database.................................................")
os.system(f"""makeblastdb -in {AMR_CDS_file} -dbtype nucl -out AMR_Blast-Results-CDS/AMR""")
os.system(f"""blastn -db AMR_Blast-Results-CDS/AMR -query {input_cdhit_res} -outfmt 7 -num_threads 16 -out AMR_Blast-Results-CDS/AMR_CDS_BLASTN_results.txt""")

print("Performing Blast for ML Important Kmers..................................")
if not os.path.isdir("AMR_Kmer_BLASTN-results"):
    os.mkdir("AMR_Kmer_BLASTN-results")

for file in os.listdir(blast_seq_path):
    if file.endswith(".fasta"):
        output = file.replace(".fasta",".txt")
        file_path = os.path.join(blast_seq_path, file)
        #print(f"""blastn -db AMR_Blast-Results-CDS/AMR -query {file_path} -num_threads 16 -outfmt 7 -evalue 1000 -word_size 13 -gapopen 0 -gapextend 0 -out AMR_Kmer_BLASTN-results/{output}\n""")   
        os.system(f"""blastn -db AMR_Blast-Results-CDS/AMR -query {file_path} -num_threads 16 -outfmt 6 -evalue 1000 -word_size 13 -gapopen 0 -gapextend 0 -out AMR_Kmer_BLASTN-results/{output}""")


#create excel..............................
if not os.path.isdir("Excel-Results"):
    os.mkdir("Excel-Results") 

blast_res_dir_path = Root_dir + "/AMR-Results/AMR_Kmer_BLASTN-results"
output_Excel1 = "Excel-Results/Table-2-AMR-Antibiotic-Kmer-Results.xlsx"
ReferenceGeneCatalog = Root_dir +"/AMR_gene_SEQ/ReferenceGeneCatalog.txt"
# ReferenceGeneCatalog = Root_dir +"/CARD_Gene_Seq/Annotation"

createKMERExcel(blast_res_dir_path,output_Excel1,ReferenceGeneCatalog,antibiotics_acr= antibiotics,input_amr_path=Antibiotic_data_path,Selected_Kmer_seq_path=blast_seq_path)

clstr_file = Root_dir +"/AMR-Results/cd-hit-results/cd-hit-results.clstr"

blast_res_CDS_path = Root_dir +"/AMR-Results/AMR_Blast-Results-CDS/AMR_CDS_BLASTN_results.txt"
output = "Excel-Results/Table-1-AMR-CDS-Results.xlsx"
print(os.getcwd())
createCDSxls(blast_res_CDS_path,ReferenceGeneCatalog,input_genome_xls,clstr_file,output)

output_excel_path = Root_dir+"/AMR-Results/Excel-Results/Table-3-Kmer-AI-Results.xlsx"
Classification_Result_Excel(clf_results_folder, output_excel_path)



os.chdir(Root_dir)

