'''
This program works in conjunction with the identify_uORFs_data_wrangler.py file.
It uses ribo-seq data in order to identify peaks which could be potential uORFs
and gives each potential a uORF a normalized score based on the peak in its 
associated ORF start codon ribo-seq peak.

s = normal start score (main ORF) 
c_i = course reads at position_i in UTR 
u_i = score for a position_i in 5' UTR
'''
from find_uORFs_data_wrangler import *
from scipy.signal import find_peaks



class IdentifyuORFs:
    
    def __init__(self, GWIPS_df, fp_exons_df, coding_exons_df):
        self.GWIPS_df = GWIPS_df
        self.fp_exons_df = fp_exons_df
        self.coding_exons_df = coding_exons_df
        self.dw = DataWrangler()
    
    def find_sig_peaks(self, interval, data):
        '''
        finds peaks of the dataset over a certain interval
        filters out peaks that are not within 1 standard deviation
        of the mean
        returns the height of the peaks and its position in the genome
        '''
        pos = np.array(range(interval[0], interval[1]+1))
        count = np.array(data)
        
        peaks_0, _ = find_peaks(data)
        if len(peaks_0) == 0:
            peaks_0, = np.where(count==max(count))
#         print("peaks")
#         print(peaks_0)
        
        where_sig = (count[peaks_0] > np.mean(count[peaks_0]) + np.std(count[peaks_0]))
        
        sig_peaks = count[peaks_0[where_sig]]
        sig_peaks_pos = pos[peaks_0[where_sig]]
        
        if len(sig_peaks) == 0:
            sig_peaks = count[peaks_0]
            sig_peaks_pos = pos[peaks_0]
        
        return sig_peaks, sig_peaks_pos
    
    def _plot_peaks_(self,interval,data):
        '''
        plots the counts at all position, marks significant peaks
        with red x's
        '''
        pos = np.array(range(interval[0], interval[1]+1))
        count = np.array(data)
        plt.plot(pos, count)
        
        sig_peaks, sig_peaks_pos = self.find_sig_peaks(interval, data)
        plt.plot(sig_peaks_pos, sig_peaks, "xr")

        plt.show()
        
    
    def find_start_codon(self, chromosome, annotation="uc003ysh", bp_range=3):
        '''
        finds the start codon in the coding region exon
        scores the start codon by taking the average of the counts
        in the range of adding and subtracting the bp_range 
        (ie. pos 100 will have range 97-103) from the position of the
        start codon
        
        inputs:
            annotation - the name of the annoation we are using
            bp_range - the number of basepairs to right and left of the 
                        start codon we want to take the average counts
                        for the start codon score
        outputs:
            start_codon_score - the score assigned to the start codon in the coding exon
            start_codon_pos - the estimated position of the start codon in the coding exon
        '''
        #checks the direction of the sequence
        #---will all annotations have reads in the same direction???---
        exon_range_list = []
        df = self.coding_exons_df.get_group(annotation)
        for i in range(len(df['pos_start'].tolist())):
                exon_range_list.append([df['pos_start'].tolist()[i],df['pos_end'].tolist()[i]])
        
#         print(exon_range_list)
        if "+" in df['dir'].tolist():
            #in pos dir we look for the first exon
            first_exon = min(exon_range_list)
        else:
            #in the neg dir we look for the last exon
            first_exon = max(exon_range_list)
        
        first_exon_ribo_seq_df = self.GWIPS_df[["chr","pos_start","count"]][(self.GWIPS_df["pos_start"] >= first_exon[0]) & 
                                                                            (self.GWIPS_df["pos_end"] <= first_exon[1]+1) & 
                                                                            (self.GWIPS_df['chr']==chromosome)]
                
        #including positions with zero counts
        first_exon_ribo_seq_complete_df = self.dw.ribo_seq_add_zeros(first_exon_ribo_seq_df,first_exon,chromosome)

        first_exon_count = self.dw.flatten_list(first_exon_ribo_seq_complete_df[['count']].values.tolist())
        
        #debugging
#         print("FIRST EXON")
#         print(first_exon)
#         print(first_exon_ribo_seq_complete_df.head(3))
#         print(first_exon_ribo_seq_complete_df.tail(3))
#         print(len(first_exon_count))
        #debugging end
    
        #handles exception when a protein does not have ribo-seq data
        if all(count == 0 for count in first_exon_count):
            start_codon_score = 0
            start_codon_position = ''
            return start_codon_score, start_codon_position
        
        #finds significant peaks and their position
        sig_peaks, sig_peaks_pos = self.find_sig_peaks(first_exon, first_exon_count)
    
        #plotting the peaks for reference, remove in final code
#         self._plot_peaks_(first_exon, first_exon_count)
        
        
        
        #limits the data frame to the bp_range of the start codon position
#         print("sig peaks pos")
#         print(sig_peaks_pos)
#         print(sig_peaks)
        start_codon_score_range_df = first_exon_ribo_seq_complete_df[['pos_start','count']][(first_exon_ribo_seq_complete_df['pos_start']>= sig_peaks_pos[0]-bp_range)&
                                                                                            (first_exon_ribo_seq_complete_df['pos_start']<= sig_peaks_pos[0]+bp_range)]
        
        #takes the average of the counts over the bp_range of the start codon to find the score
        start_codon_score = np.mean(start_codon_score_range_df['count'])
        start_codon_pos = sig_peaks_pos[0]
        
#         print("Start Codon Position = ", start_codon_pos)
        return start_codon_score, start_codon_pos
    
    def score_UTR(self, chromosome, annotation="uc003ysh", bp_range=3):
        '''
        scores all positions in each 5' UTR for the annotation
        returns a dataframe with the start poisition and score
        
        score = pos_height / start_codon_score
        '''
        scores = []
        positions = []
#         print("score_UTR")
        
        start_codon_score, start_codon_pos = self.find_start_codon(chromosome, annotation, bp_range)
        if start_codon_score == 0:
            GWIPS_score_df = "No ribosome data available in genome browser"
            return GWIPS_score_df,start_codon_score
        
        
        
        df = self.fp_exons_df.get_group(annotation)
        
#         print("5' UTR Positions")
#         print(df)
        
        init_pos = min(df['pos_start'].tolist())
        term_pos = max(df['pos_end'].tolist())
        
#         print("init pos =", init_pos)
#         print("term pos = ", term_pos)
        
        fp_ribo_seq_df = self.GWIPS_df[['chr', 'pos_start', 'pos_end', 'count']][(self.GWIPS_df['pos_start'] >= init_pos) & 
                                                                                 (self.GWIPS_df['pos_start'] <= term_pos) & 
                                                                                 (self.GWIPS_df['chr'] == chromosome)]
                
        pos_start_list = df['pos_start'].tolist()
        pos_end_list = df['pos_end'].tolist()
        fp_ribo_seq_pos_list = fp_ribo_seq_df['pos_start'].tolist()
        
        
        for i in range(len(df['pos_start'])):#iterates over the position of each start position in the 5' UTR dataframe
            for pos in range(pos_start_list[i],pos_end_list[i]):#iterates over the range of the start and end positions of each 5' UTR
                #the every gene position in all 5'UTRs to the positions list
                if pos in fp_ribo_seq_pos_list:#if the position is in the riboseq dataframe, calculate the score, else the score is 0
                    scores.append(int(fp_ribo_seq_df.loc[(fp_ribo_seq_df['pos_start'] == pos),'count'])/start_codon_score)
                    positions.append(pos)
#                 else:
#                     scores.append(0)
                
        GWIPS_score_df = pd.DataFrame({"start_pos" : positions, "score" :scores},columns = ["start_pos","score"])#contains all the positions in all 5'UTRs and their associate score
        
        return GWIPS_score_df,start_codon_score
    
    def compare_coding_region_f_p_UTR(self, chromosome, annotation="uc003ysh", bp_range=3):
        '''
        runs the find start codon method
        finds the riboseq data for all 5' UTRs 
        returns the start codon riboseq height and the riboseq data for the 5' UTR
        '''
        scores = []
        positions = []
#         print("score_UTR")
        
        start_codon_score, start_codon_pos = self.find_start_codon(chromosome, annotation, bp_range)
        if start_codon_score == 0:
            GWIPS_score_df = "No ribosome data available in genome browser"
            return GWIPS_score_df,start_codon_score
        
        
        
        df = self.fp_exons_df.get_group(annotation)
        
#         print("5' UTR Positions")
#         print(df)
        
        init_pos = min(df['pos_start'].tolist())
        term_pos = max(df['pos_end'].tolist())
        
#         print("init pos =", init_pos)
#         print("term pos = ", term_pos)
        
        fp_ribo_seq_df = self.GWIPS_df[['chr', 'pos_start', 'pos_end', 'count']][(self.GWIPS_df['pos_start'] >= init_pos) & 
                                                                                 (self.GWIPS_df['pos_start'] <= term_pos) & 
                                                                                 (self.GWIPS_df['chr'] == chromosome)]
        
        
        return fp_ribo_seq_df, start_codon_score
    
    def find_potential_uORFs(self, chromosome, annotation, bp_range=3):
        '''
        From the scores of all the UTRs for a protein, returns the genome positions that have the 
        '''
        GWIPS_score_df,start_codon_score = self.score_UTR(chromosome, annotation, bp_range) #returns a dataframe with the scores for each position in all UTRs for a protein
        if start_codon_score == 0:
            uORF_candidate_pos = "No ribosome data available in genome browser"
            print("No ribosome data available in genome browser")
            return uORF_candidate_pos
            
#         print("Score Summary \n",GWIPS_score_df['score'].describe())        
        
        #==Used to check if there are any scores greater than 0 in the score dataframe==
#         greater_than_0_df = GWIPS_score_df[['start_pos','score']][GWIPS_score_df['score'] > 0]# creates a datafram with all positions that have a score greater than 0
        
#         plt.hist(greater_than_0_df['start_pos'],bins=50)
#         print("HIST OF GREATER THAN 0")
#         plt.show()
        #====
    
        #sets the score threshold to be greater than the 75% value, if the 75% value is 0, then sets the threshold to 1/(the start codon score) --- the minimum possible score for a UTR
        if GWIPS_score_df['score'].describe()["75%"] > 0:
            score_threshold = GWIPS_score_df['score'].describe()["75%"]
        else:
            score_threshold = 1/start_codon_score
        
        
        UTR_candidates_df = GWIPS_score_df[['start_pos','score']][GWIPS_score_df['score'] >= score_threshold] #creates a dataframe of positions and scores where the scores are larger than the threshold
        
        #we use histogram below to generalize the peaks over more positions. this 
#         plt.plot(UTR_candidates_df['start_pos'],UTR_candidates_df['score'])
#         plt.show()
#         print("=====>")
#         plt.hist(UTR_candidates_df['start_pos'],bins=50)
#         plt.show()
        
        #creates a list of positions when the positions are binned
        pos = np.histogram(UTR_candidates_df['start_pos'],bins=50)[1]
        #creates a list of how many scores for each of the bins
        count = np.histogram(UTR_candidates_df['start_pos'],bins=50)[0]

        if max(count) == 0:
            print("No potential uORFs")
            uORF_cadidate_pos = []
            return uORF_cadidate_pos
        
        #identifies the peaks in the histogram
        peaks, _ = find_peaks(np.histogram(UTR_candidates_df['start_pos'],bins=50)[0])
        
        #if there are no peaks found, then all positions with a count greater than 0 are considered peaks
        if len(peaks) == 0:
            ps = []
            for i in range(len(np.histogram(UTR_candidates_df['start_pos'],bins=50)[0])):
                if np.histogram(UTR_candidates_df['start_pos'],bins=50)[0][i] > 0:
                    ps.append(i)
            peaks = np.array(ps)

        #creates a condition where_sig is a value greater than one standard deviation from the mean
        where_sig = (count[peaks] > np.mean(count[peaks]) + np.std(count[peaks]))
        #defines the list of counts for significant peaks
        sig_peaks = count[peaks[where_sig]]
        #defines the list of positions for significant peaks, these are the uORF candidate positions
        sig_peaks_pos = pos[peaks[where_sig]]
        plt.hist(UTR_candidates_df['start_pos'],bins=50)
        if len(sig_peaks) == 0:
            plt.plot(pos[peaks], count[peaks], "xr")
            plt.show()
            uORF_candidate_pos = pos[peaks]

        else:
            plt.plot(sig_peaks_pos, sig_peaks, "xr")
            plt.show()
            uORF_candidate_pos = sig_peaks_pos
        
        return uORF_candidate_pos
    
    def validation_plot(self, chromosome, annotation, bp_range=3):
        exon_range_list = []
        exon_start_pos = []
        exon_stop_pos = []
        df = self.coding_exons_df.get_group(annotation)
        for i in range(len(df['pos_start'].tolist())):
                exon_range_list.append([df['pos_start'].tolist()[i],df['pos_end'].tolist()[i]])
                exon_start_pos.append(df['pos_start'].tolist()[i])
                exon_stop_pos.append(df['pos_end'].tolist()[i])
        
        print(exon_range_list)
        print(max(exon_stop_pos))
        print(min(exon_start_pos))
        
        #plots the scores/ribosome-counts for each position in the 5'UTRs --- same as genome browser, 
        #can be used to check if there are ribosomes on the 5'UTRs / if thereare any potential uORFs in the 5'UTRs
        plt.plot(GWIPS_score_df['start_pos'],GWIPS_score_df['score'])
        plt.show()
        
        #compares the UTR score scatter plot histogram to see how peaks in the score data are genearlized
        plt.plot(UTR_candidates_df['start_pos'],UTR_candidates_df['score'])
        plt.show()
        print("=====>")
        plt.hist(UTR_candidates_df['start_pos'],bins=50)
        plt.show()

    
  
        
