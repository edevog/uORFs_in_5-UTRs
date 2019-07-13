import numpy as np
import pandas as pd



class DataWrangler:
    
    def exon_data(self, file_path):
        '''
        file_path = the file path for exon data from genome browser in BED format
        puts exon data files into data frames, separates annotation and exon names
        into separate columns, groups the data frame by annotation
        '''
        exon_path = file_path
        exon_headers = ['chr', 'pos_start', 'pos_end', 'annotation', '0s', 'dir']
        exon_df = pd.read_csv(exon_path, sep="\t", header=None, names=exon_headers)
        exon_df["annotation"] = exon_df["annotation"].str.split('_', 1)
        exon_df[['annotation','exon_name']] = pd.DataFrame(exon_df.annotation.values.tolist(), index= exon_df.index)
        
        #removes the column of 0's and reorders the columns
        exon_df = exon_df[['chr', 'pos_start', 'pos_end', 'annotation','exon_name', 'dir']]
        exon_df_by_ann = exon_df.groupby('annotation')
        
        return exon_df_by_ann

    def gencode_data(self, file_path):
        exon_path = file_path
        exon_headers = ['chr', 'pos_start', 'pos_end', 'annotation', '0s', 'dir']
        exon_df = pd.read_csv(exon_path, sep="\t", header=None, names=exon_headers)
        exon_df["annotation"] = exon_df["annotation"].str.split('.', 1)
        exon_df[['annotation','exon_name']] = pd.DataFrame(exon_df.annotation.values.tolist(), index= exon_df.index)
        
        #removes the column of 0's and reorders the columns
        exon_df = exon_df[['chr', 'pos_start', 'pos_end', 'annotation','exon_name', 'dir']]
        exon_df_by_ann = exon_df.groupby('annotation')
        
        return exon_df_by_ann
    
    def ribo_seq_data(self,GWIPS_path):
        '''
        puts position and count data from GWIPS genome browser into a dataframe
        '''
        GWIPS_headers = ['chr', 'pos_start', 'pos_end', 'count']
        GWIPS_df = pd.read_csv(GWIPS_path, sep="\t", header=None, names=GWIPS_headers)
        return GWIPS_df
    
    def protein_transcript_data(self, protein_transcript_path):
        '''
        creates a data frame of the protein:transcript data file where the transcript id's are 
        in a list
        '''
        protein_transcript_headers = ['protein', 'protein_id', 'transcript_id']
        protein_transcript_df = pd.read_csv(protein_transcript_path, sep="\t", header=None, names=protein_transcript_headers, engine='python')
        #splits the transcript id's by ; and puts them into a list
        protein_transcript_df["transcript_id"] = protein_transcript_df["transcript_id"].str.split(";")
        
        return protein_transcript_df

    def pcawg_ref_uorfs(self, file_path):
        
        headers = ['chr', 'start_pos', 'end_pos', 'protein', '.s', 'dir']
        df = pd.read_csv(file_path, sep="\t", header=None, names=headers)
        df["protein"] = df["protein"].str.split(':', 1)
        # print(df)

        df[['protein','5_p_exon_count']] = pd.DataFrame(df.protein.values.tolist(), index= df.index)
        
        #reorders the columns
        df = df[['chr', 'start_pos', 'end_pos', 'protein','5_p_exon_count', '.s', 'dir']]

        return df
    
    def ribo_seq_add_zeros(self,GWIPS_df,interval,chromosome):
        '''
        finds missing positions in the range and adds a row with that position
        and a count of zero
        assuming missing positions are 0-counts????
        '''
#         print(GWIPS_df['pos_start'].tolist())
        for i in range(interval[0],interval[1]+1):
            if i not in GWIPS_df['pos_start'].tolist():
#                 print(i)
                GWIPS_df = GWIPS_df.append(pd.DataFrame({'chr':[chromosome],'pos_start':[i],'count':[0]}),ignore_index = True)
        GWIPS_df = GWIPS_df.sort_values(by=["pos_start"])
        return GWIPS_df
                
    
    def flatten_list(self, list_o):
        '''
        inputs a list of lists, outputs a single list
        '''
        flat_list = []
        for sublist in list_o:
            for item in sublist:
                flat_list.append(item)
        return flat_list
    
    def print_grouped_df(self,data_frame):
        for key,item in data_frame:
            print(key)
            print(data_frame.get_group(key), "\n\n")

    def make_consensus_dict(self, primary_transcript_df, coding_exons_df, fp_exons_df):
        consensus_transcripts_dict = {}
        missing_transcripts = []
        genome_browser_list = []
        fp_transcript_list = []

        for key, item in coding_exons_df:
            genome_browser_list.append(key)

        for key,item in fp_exons_df:
            fp_transcript_list.append(key)

        for i in range(len(primary_transcript_df['transcript_id'])): #for the full list of the transcript ids for a single protein
            for transcript_id in primary_transcript_df['transcript_id'][i]: #a single transcript id for a protein
                if transcript_id in genome_browser_list and transcript_id in fp_transcript_list:
                    consensus_transcripts_dict[primary_transcript_df['protein'][i]]=transcript_id
                    break
                else:
                    missing_transcripts.append(transcript_id)
        return consensus_transcripts_dict

    def make_transcript_chromosome_dict(self,coding_exons_df):
        # coding_exons_transcripts = []
        transcript_chromosome_dict = {}

        for key,item in coding_exons_df:
            # coding_exons_transcripts.append(key)
            transcript_chromosome_dict[key] = item['chr'].iloc[0]
        
        return transcript_chromosome_dict




