import polars as pl
import os
import shutil
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor
import subprocess

id_path = 'id.csv' # step 1. query id path (NEC)
meta_path = '/home/nishitha/Desktop/Pylogroups/Metadata.csv' # step 1. downloaded from microreact


blast_df_path = 'result.txt' # step 3. blast result 



genomes_files_path = 'Genomes'



# ------------------------------extra commands-------------------------------
# creating a blast database -- makeblastdb -in genomes.fna -dbtype nucl -out blastdatabase
# concatenateing file --- cat *.fna > concat.fna
# blast cmd -- blastn -query query.fasta -db blastdatabase -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen" -out result.txt
# installing ANIclustermap ---- pip install aniclustermap
# ANIcluster cmd -------  ANIclustermap -i /home/nishitha/Desktop/Pylogroups/best_match -o /home/nishitha/Desktop/Pylogroups/ANI --fig_width 15 --cmap_colors white,orange,red --dendrogram_ratio 0.25




################################################# --- download the reference data ---- ######################################

def downlad_reference_data(meta_path):
    df =pl.read_csv(meta_path,ignore_errors = True)


    def download(i):

        try : 
            
            cmd = 'curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/'+i+'/download?include_annotation_type=GENOME_FASTA > '+i+".zip" 
            ## cmd = curl https://www.ebi.ac.uk/ena/browser/api/fasta/+i+?download=true&gzip=true > i ,,
            subprocess.run(cmd,shell=True)
            print(f'\n\n Downloaded ... {i}\n\n')
        except:
            print(f'\nFailed to download {i}\n')


    def unzip():
        os.makedirs('unzipped',exist_ok=True)
        zipfiles = [i for i in os.listdir(os.getcwd()) if i.endswith(".zip")]

        for i in zipfiles:
            dir = i[:i.find('.zip')]
            cmd = f'unzip -d unzipped/{dir} {i}'
            subprocess.run(cmd,shell=True)
        
        subprocess.run('rm -f *.zip',shell=True)
        
    def genomes(path):
        os.makedirs('Genomes',exist_ok=True)
        fna = []
        for root, dirs, files in os.walk(path):

            for i in files:
                if i.endswith('.fna'):
                    fna.append(os.path.join(root,i))
        for j in fna:
            shutil.copy(j,f'Genomes/')

        subprocess.run('rm -rf unzipped', shell=True)
        
        
        
    def process_data_down(df):

        with ThreadPoolExecutor() as executor:
            futures = []

            for i in df['id'].sort():
                futures.append(executor.submit(download, i) )


    

    def collect(path, df):
        new_df = pl.DataFrame()

        files = [i for i in os.listdir(path) if i.endswith('.fna')]

        for i in files:
            temp = pl.DataFrame()
            temp = temp.with_columns(id = pl.lit(i[:i.find('_',i.find('_')+1)]))
            new_df = pl.concat([new_df,temp])
        
        print(new_df)

        print(df.join(new_df, on='id', how='anti'))
        not_pre = df.join(new_df, on='id', how='anti')
        not_pre.write_csv('not_downloaded.csv')
        
        print('\n\n............... reference data downloaded ...............\n\n')
        
        
    process_data_down(df)
    unzip()
    genomes('unzipped')
    collect('Genomes',df)
################################################################################################################################



def query_genome_download(id_path):
    df =pl.read_csv(id_path,ignore_errors = True)
    os.makedirs('query_genomes',exist_ok=True)
    for i in df['id']:
        with open(f'query_genomes/{i}.fna', 'w') as fasta_file:

            fetch_handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
            record = fetch_handle.read()
            fasta_file.write(record)
            fetch_handle.close()
    print('\n\n............... query genome downloaded ...............\n\n')


######################################### header processing ########################################################
def process_header(i, input_dir, output_dir):

    
    
    with open(f'{input_dir}/'+i,'r') as IN:
        file = IN.read()
        file = file.replace(">", ">"+i[:i.find('_',i.find('_')+1)]+'|')
    with open(f'{output_dir}/'+i,'w') as OUT:
        OUT.write(file)
        
    



def process_data_header(genomes_files_path,output_dir,extention):
    fasta = [ i for i in os.listdir(genomes_files_path) if i.endswith(f'{extention}')]
    os.makedirs(f'{output_dir}',exist_ok=True)
    with ThreadPoolExecutor() as executor:
        futures = []

        for i in fasta:

            futures.append(executor.submit(process_header, i, genomes_files_path, output_dir) )

    print('\n\n............... Header renaming completed ...............\n\n')


####################################################################################################################
    
def concat_genomes(genomes_files_path,extention,out):
    os.makedirs('blast_run/',exist_ok=True)
    cmd = f'cat {genomes_files_path}/*{extention} > blast_run/{out}.fna'
    subprocess.run(cmd,shell=True)
    subprocess.run('rm -rf contig_renamed_genomes',shell=True)

    print('\n\n............... Geneomes concatenated ...............\n\n')

def hit_extraction(blast_df_path, meta_path , genomes_files_path):
    
    nec_meta = pl.DataFrame()
    genomes = [os.path.join(genomes_files_path,i) for i in os.listdir(genomes_files_path)]

    files = pl.DataFrame({'files' : genomes})

    df = pl.read_csv(blast_df_path, separator='\t',has_header=False,ignore_errors=True)
    print(df)

    meta = pl.read_csv(meta_path,ignore_errors=True)
    meta = meta.sort('Phylogroup','SequenceScore', descending=True)
    
    for i in df['column_1'].unique().sort():
    
        temp = df.filter(pl.col('column_1') == i)
        temp = temp.sort('column_4',descending=True)[:1]
        acc = temp['column_2'].to_list()[0]
        acc = acc[:acc.find('|')]

        phylo = meta.filter(pl.col('id') == acc)['Phylogroup'].to_list()[0]

        extracted = meta.filter((pl.col('id') == acc))

        print(f'{i} - {phylo}')

        temp_nec_meta = extracted.select(('id','Phylogroup'))
        temp_nec_meta = temp_nec_meta.with_columns(nec = pl.lit(i))

        nec_meta = pl.concat([nec_meta,temp_nec_meta])


    nec_meta.write_csv('NEC_data.tsv',separator='\t')

    print(nec_meta.select(('id')).unique())

    os.makedirs('best_match',exist_ok=True)

    for i in nec_meta['id'].unique():
        g = files.filter(pl.col('files').str.contains(i))['files'].to_list()[0]
        print(g)
        shutil.copy(f'{g}',f'best_match/')
        



def main():
    # step 1
    # downlad_reference_data(meta_path)
    # process_data_header(genomes_files_path,'contig_renamed_genomes','.fna')
    # concat_genomes(genomes_files_path='contig_renamed_genomes',extention='.fna',out='genomes')
    # query_genome_download(id_path)
    
    #step 2 --- manual
    ## rename query genome manually
    ## blast
     
    
    #step 3
    hit_extraction(blast_df_path, meta_path , genomes_files_path)
    
    # step 4 --- manual
    #ANIcluster map generation
main()

