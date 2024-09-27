import pandas as pd

# add path here

NEC1_gene_result = "NEC1_gene_result.txt"   # blast results
NEC1_All_Variations = "NEC1_All_Variations.vcf"  # variation file
df_path = "NEC1_DrugResistant.xlsx"

des_path = "/home/Desktop/nishitha/Project/final/ECOLI/Gene_description.xlsx"
file_path = "/home/Desktop/nishitha/Project/final/ECOLI/bacteria.gb"


# file = open(NEC1_gene_result,"r").read()

#save as excel file
header = ["qseqid","sseqid","pident","length","mismatch","gaps","qstart","qend","sstart","send","evalue","bitscore","qlen"]

# save txt data as it is in excel
df = pd.read_csv(NEC1_gene_result,sep="\t",names=header)


# df.to_excel('output.xlsx', index=False)
# shee2 with the pident is greater than equals to 99
df1 = df[df["pident"] >= 99]

#split sseqid into two more columns
df1[["dbname","gene"]] = df1["sseqid"].str.split('|', expand=True)
df1["gene"] = df1["gene"].str.upper()

#separte by _ select second  remove - :
cleaned_gene = []
for i in df1["gene"]:
   if "_" in str(i):
      l = i.split("_")
      cleaned_gene.append(l[1])
   elif '-' in str(i):
      i = i.replace('-','')
      print(f"find {i}")
      i = ''.join(i)
      cleaned_gene.append(i)
   else:
      cleaned_gene.append(i)     
   
df1["gene"] = cleaned_gene     

# select only those row which start with 'MEG_'
MEG = df1[df1['dbname'].str.startswith('MEG_')]
MEG_final = MEG.sort_values("length")
MEG_final = MEG_final.sort_values('pident').drop_duplicates('gene',).sort_index()

# select only those row which start with 'CARD_gb'
CARD_gb = df1[df1['dbname'].str.startswith('CARD_gb')]
CARD_gb_final = CARD_gb.sort_values("length")
CARD_gb_final = CARD_gb_final.sort_values('pident').drop_duplicates('gene',).sort_index()

# select only those row which start with 'ResFinder'
ResFinder = df1[df1['dbname'].str.startswith('ResFinder')]
ResFinder_final = ResFinder.sort_values("length")
ResFinder_final = ResFinder_final.sort_values('pident').drop_duplicates('gene',).sort_index()
#print(ResFinder.head())



#for final need to append discription 
#open the discription df
#des_path = "Gene_description.xlsx"
description = pd.read_excel(des_path)

description["Description"] = description["Description"].str.replace(">",'')

description.rename(columns={"Gene" : "sseqid"},inplace=True)
con = [MEG_final,CARD_gb_final,ResFinder_final]
final = pd.concat(con)
dfinal = description.merge(final, on="sseqid", how = 'inner')

meg = []
MEG_gene = [i.lower() for i in MEG_final["gene"]]
dfinal_gene = [i.lower() for i in dfinal["gene"]]


for i in dfinal_gene:
        if i in  MEG_gene :
           meg.append("yes")
        else:
            meg.append("no")    

card = []
card_gene = [i.lower() for i in CARD_gb_final["gene"]]

for i in dfinal_gene:
        if i in  card_gene :
           card.append("yes")
        else:
            card.append("no")    


rest = []
rest_gene = [i.lower() for i in ResFinder_final["gene"]]

for i in dfinal_gene:
        if i in  rest_gene :
           rest.append("yes")
        else:
            rest.append("no")    

#final sheet
dfinal['MEG'] = meg
dfinal['CARD'] = card
dfinal['RESTFINDER'] = rest
dfinal = dfinal.sort_values('pident').drop_duplicates('gene',).sort_index()
dfinal = dfinal.drop_duplicates(subset=["qstart","qend","qend","length"]).sort_index()


new_cols = ['qseqid','sseqid', 'Description', 'pident', 'length', 'mismatch',
       'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
       'qlen', 'dbname', 'gene', 'MEG', 'CARD', 'RESTFINDER']
dfinal=dfinal[new_cols]


# read 'file ................


# path = "NEC1_All_Variations.vcf"
col = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NEC1"]

#NEC1_All_Variations = "NEC1_All_Variations.vcf"
file = pd.read_csv(NEC1_All_Variations, sep="\t", comment='#',names=col)
# print(file.iloc[1,7])
# filter data based on upsream and down
new_file = file[~file["INFO"].str.contains("upstream_gene_variant")] 
new_file2 = new_file[~new_file["INFO"].str.contains("downstream_gene_variant")]
#print(new_file2.iloc[:5,7])

# Create an empty DataFrame to store the matched values
df_matched = pd.DataFrame(columns=['Description','sseqid', 'qseqid', 'pident', 'length', 'mismatch',
       'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
       'qlen', 'dbname', 'gene', 'MEG', 'CARD', 'RESTFINDER',"CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NEC1"])

# Iterate over each row in the first dataset
for index, row in new_file2.iterrows():
    value = row['POS']
    
    # Check if the value is between any range in the second dataset
    matched_rows = dfinal[(value >= dfinal['qstart']) & (value <= dfinal['qend'])]
    
    # If there are matching rows, append them to the matched DataFrame
    if not matched_rows.empty:
        # Assign values from the current row to the matched rows using .loc
        matched_rows.loc[:, 'CHROM'] = row['CHROM']
        matched_rows.loc[:, 'POS'] = row['POS']
        matched_rows.loc[:, 'ID'] = row['ID']
        matched_rows.loc[:, 'REF'] = row['REF']
        matched_rows.loc[:, 'ALT'] = row['ALT']
        matched_rows.loc[:, 'QUAL'] = row['QUAL']   
        matched_rows.loc[:, 'FILTER'] = row['FILTER']
        matched_rows.loc[:, 'INFO'] = row['INFO'] 
        matched_rows.loc[:, 'FORMAT'] = row.loc['FORMAT']
        matched_rows.loc[:, 'NEC1'] = row['NEC1']
        
        df_matched = pd.concat([df_matched, matched_rows])  # Concatenate matched rows to df_matched
        

        
# Reset the index of df_matched and drop the previous index
df_matched.reset_index(drop=True, inplace=True)





# add mutation column in dfinal.....................................................................

mutation = []
qstart = dfinal['qstart'].to_list()
qend = dfinal["qend"].to_list()

for i in range(0,len(dfinal["qstart"])):
   if qstart[i] in df_matched["qstart"].values and qend[i] in df_matched["qend"].values:
      mutation.append("yes")
   else:
      mutation.append("no")
      
dfinal["Mutation"] = mutation      

#gen categories.......................................................................................
#/Users/nishitha/Dropbox/My Mac (Nishitha MacBook Pro)/Desktop/Mutation in AMR genes/ECOLI/Final_Ecoli_gene_category_Final.xlsx
gen_path = "/Users/nishitha/Dropbox/My Mac (Nishitha MacBook Pro)/Desktop/Mutation in AMR genes/ECOLI/Final_Ecoli_gene_category_Final.xlsx.xlsx"
Gen_categories = pd.read_excel(gen_path)
Gen_categories["gene"] = [i.upper() for i in Gen_categories['Gene']]

newfinal = Gen_categories.merge(dfinal, on="gene", how = 'inner')


# try:
#    newfinal.drop(columns="Gene",inplace=True)
# except Exception as e:
#    print("Gene Not Found.....")   
print(newfinal.columns)
new_cols = ['qseqid','sseqid', 'Description', 'pident', 'length', 'mismatch',
       'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
       'qlen', 'dbname', 'gene', 'MEG', 'CARD', 'RESTFINDER',"Mutation",'Gene_Category',"category","Color"]

newfinal=newfinal[new_cols]
      







#creating excel file 
with pd.ExcelWriter('NEC1_DrugResistant.xlsx') as writer:  # doctest: +SKIP
   df.to_excel(writer, sheet_name='toexcel',index=False)
   df1.to_excel(writer, sheet_name='pident>=99',index=False)
   MEG_final.to_excel(writer, sheet_name='MEG',index=False)
   CARD_gb_final.to_excel(writer, sheet_name='CARD_gb9',index=False)
   ResFinder_final.to_excel(writer, sheet_name='ResFinder',index=False)
   df_matched.to_excel(writer, sheet_name='Mutation',index=False)

   newfinal.to_excel(writer, sheet_name='final',index=False)
   print("Excel file Saved succusesfuly .......")
   

#path
# file_path = "/Users/nishitha/Dropbox/My Mac (Nishitha MacBook Pro)/Desktop/Mutation in AMR genes/ECOLI/bacteria.gb"
# df_path = "NEC1_DrugResistant.xlsx"

# #access the genbank file first .......
# with open(file_path) as f:
#     f_out = None
#     for line in f:
#         if line.startswith('ORIGIN'):      # we need a new output file       
#             if f_out:
#                 f_out.close()
#             f_out = open(f'{"ORIGIN"}.gb', 'w')
#         if f_out:
#             f_out.write(line)
#     if f_out:
#         f_out.close()
 
import pandas as pd
# open axcel file
df_path = "/home/Desktop/nishitha/Project/struct/NEC1/NEC1_DrugResistant.xlsx"
xls =pd.ExcelFile(df_path)
df = pd.read_excel(xls,"final")


newdf = df[df["Mutation"].str.contains("yes")]
newdf["Gene_Category"].str.strip()
print(newdf.iloc[:5,21])
#sorted........................................
newdf = newdf.sort_values(by=['qstart','qend'], ascending=[True, True])
#find out the length
length = (newdf["qend"] - newdf["qstart"]).to_list()

qstart = [50]

x = length[0]
qend = [2000]


for i in range(0,len(length)-1):
   s = qend[i]
   s = s+100
   qstart.append(s)
   
   qe = s + (length[i+1]*2)
   qend.append(qe)
   
# create  a list , for complement or not ....................
comp = []
for i in range(0,len(newdf["Mutation"])):
    if newdf.iloc[i,9] >= newdf.iloc[i,10]:
        comp.append(f"""{qstart[i]}..{qend[i]}""")
    else:
        comp.append(f"""complement({qstart[i]}..{qend[i]})""") 
       
#open gb file 
file = open("newbact.gb","w")
        

locus =f"""LOCUS       NEC1_Assembly               {qend[-1]-12} bp    DNA     linear   UNC 19-JUN-2023
FEATURES             Location/Qualifiers"""

file.write(locus)



     
for i in range(0,len(newdf["Mutation"])):
     cds =f"""
     {newdf.iloc[i,21]}             {comp[i]}
                     /note="{newdf.iloc[i,15]}"
                     /gene="{newdf.iloc[i,15]}"
                     /inference="XXXXXXX"
                     /locus_tag="{newdf.iloc[i,1]}"
                     /product="{newdf.iloc[i,15]}" """
     file.write(cds)           

#access the genbank file first .......
with open(file_path) as f:
    f_out = None
    for line in f:
        if line.startswith('ORIGIN'):      # we need a new output file       
            if f_out:
                f_out.close()
            f_out = open(f'{"ORIGIN"}.gb', 'w')
        if f_out:
            f_out.write(line)
    if f_out:
        f_out.close()  
import os

print(round(50819/60)+1)

# Define the maximum line count last qend 
max_lines = round((qend[-1]/60)+1)


file_path = "ORIGIN.gb" # Replace with the actual file path

# Read the contents of the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Split the lines into two parts
part1 = lines[:max_lines]
# Write the first part to a new file
with open('neworigin.gb', 'w') as file:
    file.writelines(part1)
    file.write("//")
 
 # Reading data from first file 
with open('newbact.gb') as fp: 
    data = fp.read() 
with open('neworigin.gb') as fp: 
        data2 = fp.read() 
# Merging two files into one another file 
data += "\n"
data += data2 
with open ('Final.gb', 'w') as fp: 
   fp.write(data)                      
    
                



os.system("rm ORIGIN.gb")
os.system("rm newbact.gb")
os.system("rm neworigin.gb")
print("Final.gb is created")

context = dict(zip(newdf.category,newdf.Color))


pyfile = open("plot.py","w")
plot = f"""


from dna_features_viewer import BiopythonTranslator

class CustomTranslator(BiopythonTranslator):
 # Label fields indicates the order in which annotations fields are
 # considered to determine the feature's label
  label_fields = ["label", "note", "name", "gene"]
  def compute_feature_legend_text(self, feature):
    return feature.type

  def compute_feature_color(self, feature):
    return{context}[feature.type]

  def compute_feature_box_color(self, feature):
   return "white"

  def compute_feature_box_linewidth(self, feature):
   return 0
translator = CustomTranslator()
graphic_record = translator.translate_record("Final.gb")
ax, _ = graphic_record.plot(figure_width=30,figure_height = 5,strand_in_label_threshold=7)
graphic_record.plot_legend(ax=ax, loc=1, ncol=3, frameon=False)
ax.figure.savefig("A_linear_plot.svg", bbox_inches="tight") 


"""
pyfile.write(plot)
pyfile.close()
os.system("python3 /home/Desktop/nishitha/Project/struct/NEC1/plot.py") 
   
   
