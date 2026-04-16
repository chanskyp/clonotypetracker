import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import random 
import os
import asyncio

from shiny import App, render, ui, reactive, req

# User interface (UI) definition
app_ui = ui.page_fluid(
    # Add a title to the page with some top padding
    ui.panel_title(ui.h2("Acquired Clonotype Tracking Plotter", class_="pt-5"),
                   #ui.h4("Acquired Clonotype Tracking Plotter", class_="text-muted")
                   ),
    ui.tags.p("TCR-beta chains. Upload your product and patient sample data to visualize the frequency of acquired clonotypes over time post-infusion. Samples should be uploaded in the following order: product, day of infusion, early post-infusion, late post-infusion. Data download functionality not optimized for Chrome, try Safari if issues arise!"),
    # User input file source/data source
    ui.input_select("data_source", "Select a data source:",choices = ["immunoSEQ TCR-beta","10X Genomics VDJ"]),
    
    # User input for using AA seq vs NT seq
    ui.input_select("sequence_type", "Select if you want to compare clonotypes using the amino acid sequence or nucleotide sequence:",choices = ["Amino Acid","Nucleotide"]),
    
    # User input files
    ui.input_file(id="file_upload", label="Choose CSV or TSV Files", accept=[".csv",".tsv"], multiple=True),

    # A numeric input for the number of days post-infusion the first sample was taken
    ui.input_numeric("timepoint1", "Days Post Infusion for Sample 1 (e.g.,35)", value=35),
    
    # A second numeric input for the number of days post-infusion the second sample was taken
    ui.input_numeric("timepoint2", "Days Post Infusion for Sample 2 (e.g.,180)", value=180),
    
    
        # A container for plot output
        
    ui.layout_columns(
        ui.output_data_frame("table_results"),
        ui.output_plot("plot_results")),
    
    ui.download_button("download_data", "Download the acquired clonotypes table"),
 
)

# Using AA seq-- update to a user option for AA or NT
def server(input,output,session): #, output, session
    
   
    @reactive.calc
    def clone_tracking():
        # input.file1() returns a list of dictionaries or None
        file_info = req(input.file_upload())
        if file_info is None:
            return None
        
    
        aa_or_nt = input.sequence_type()
        print(file_info[0]["type"]) # Check the file type of the first uploaded file
        # Check if the file input is CSV or TSV-- this is typically corresponding to 10X vs Immunoseq but some users may convert their tsv to csv already.
        if file_info[0]["type"] == "text/csv": # confirm this is what gets returned for a csv file
            # Access the first uploaded file's temporary path
            # 2. Read the CSV using the 'datapath'
            # file_info is a list of dicts, get the first one
            df_product= pd.read_csv(file_info[0]["datapath"])
            df1= pd.read_csv(file_info[1]["datapath"])
            df2= pd.read_csv(file_info[2]["datapath"])
            df3= pd.read_csv(file_info[3]["datapath"])
        elif file_info[0]["type"] == "text/tab-separated-values":
            df_product= pd.read_csv(file_info[0]["datapath"],sep='\t')
            df1= pd.read_csv(file_info[1]["datapath"],sep='\t')
            df2= pd.read_csv(file_info[2]["datapath"],sep='\t')
            df3= pd.read_csv(file_info[3]["datapath"],sep='\t')
        else: print("Invalid file input, please upload either .csv or .tsv files.")
        
        def transform_10x_data(df):
            # First edit the raw files ! Check if they are 10X or adaptive
            temp_df = df.copy() # When I didn't make a copy I had data loss problems, but I shouldn't need to do this. Recheck later.
            # Create a fake row containing 4 chains so that when we split we don't get an error. We'll drop this later before matching
            #clonotype_id	frequency	proportion	cdr3s_aa	cdr3s_nt	inkt_evidence	mait_evidence
            empty_row = pd.DataFrame({'clonotype_id': ['NA'],'frequency': [0],'proportion': [0.0],'cdr3s_aa': ['NA:NA;NA:NA;NA:NA;NA:NA'], 'cdr3s_nt': ['NA:NA;NA:NA;NA:NA;NA:NA'],'inkt_evidence': ['NA'],'mait_evidence': ['NA']})
            temp_df = pd.concat([empty_row, temp_df], ignore_index=True)
            # Split the cdr3s_aa OR cdr3s_nt on ;
            if aa_or_nt == "Amino Acid":
                temp_df[['VDJ_cdr3s_1', 'VDJ_cdr3s_2','VDJ_cdr3s_3','VDJ_cdr3s_4']] = temp_df['cdr3s_aa'].str.split(';', expand=True)
            
            elif aa_or_nt == "Nucleotide":
                temp_df[['VDJ_cdr3s_1', 'VDJ_cdr3s_2','VDJ_cdr3s_3','VDJ_cdr3s_4']] = temp_df['cdr3s_nt'].str.split(';', expand=True)
            # Chains
            temp_df[['VDJ_chain1','VDJ_cdr3s_1']] = temp_df['VDJ_cdr3s_1'].str.split(':',expand=True)
            temp_df[['VDJ_chain2','VDJ_cdr3s_2']] = temp_df['VDJ_cdr3s_2'].str.split(':',expand=True)
            temp_df[['VDJ_chain3','VDJ_cdr3s_3']] = temp_df['VDJ_cdr3s_3'].str.split(':',expand=True)
            temp_df[['VDJ_chain4','VDJ_cdr3s_4']] = temp_df['VDJ_cdr3s_4'].str.split(':',expand=True)
            
            # Drop the empty row
            temp_df = temp_df[temp_df['clonotype_id'] != 'NA']
            return temp_df
        
        def transform_adpt_data(df):
            # We don't have multiple chains resolved when using the adaptive sequencing platform so we'll set them to NaN
            temp_df = df.copy()
            print(temp_df.columns)
            if aa_or_nt == "Amino Acid":
                temp_df[['VDJ_cdr3s_1']] = temp_df[['cdr3_amino_acid']] # Shouldn't need double brackets here but I had a ValueError saying Columns must be same length as key...
            elif aa_or_nt == "Nucleotide":
                temp_df[['VDJ_cdr3s_1']] = temp_df[['cdr3_nucleotide']]
            temp_df[['VDJ_chain1']] = 'TRB'
            temp_df[['VDJ_cdr3s_2','VDJ_cdr3s_3','VDJ_cdr3s_4','VDJ_chain2','VDJ_chain3','VDJ_chain4']] = np.nan
            # By default, there is no clonotype_id column in immunoSEQ output, so we'll add one here using the row number.
            temp_df['clonotype_id'] = f"clonotype'{temp_df.index.astype(str)}" 
            return temp_df
        
        if input.data_source() == "10X Genomics VDJ":
            df_product = transform_10x_data(df_product)
            df1 = transform_10x_data(df1)
            df2 = transform_10x_data(df2)
            df3 = transform_10x_data(df3)
        elif input.data_source() == "immunoSEQ TCR-beta":
            df_product = transform_adpt_data(df_product)
            df1 = transform_adpt_data(df1)
            df2 = transform_adpt_data(df2)
            df3 = transform_adpt_data(df3)
        
        # Here we only look at the TCR-beta chain. Future iteration to allow user to select TRA vs TRB vs combo for 10X
        def trb_only(df):
            df = df[df.VDJ_chain1 == 'TRB']
            return df
        
        df_product = trb_only(df_product)
        df1 = trb_only(df1)
        df2 = trb_only(df2)
        df3 = trb_only(df3)

        # Inner join the product and patient samples
        # We join on CDR3 AA seq-- future iteration to allow user to choose either NT or AA sequence
        day0_product_overlap = df1.merge(df_product, how='inner',on=['VDJ_cdr3s_1'],suffixes=('_patient', '_product')) # The suffixes can't be too specific or the additional_overlap_checks function won't work.
        df2_product_overlap = df2.merge(df_product, how='inner',on=['VDJ_cdr3s_1'], suffixes=('_patient', '_product')) # Early post
        df3_product_overlap = df3.merge(df_product, how='inner',on=['VDJ_cdr3s_1'], suffixes=('_patient', '_product')) # Late post
        
        def additional_overlap_checks(sample_trb_overlap): #use any of the prior dataframes as the input
            # Check if additional chains and sequences match
            # _x is now _patient and _y is now _product
            #Now check if VDJ_chain2_x is NaN? continue. if VDJ_chain2_y is NaN? continue. elif,
            sample_trb_overlap_chain2_nan= sample_trb_overlap.loc[(sample_trb_overlap['VDJ_chain2_patient'].isnull()) | (sample_trb_overlap['VDJ_chain2_product'].isnull()) | (sample_trb_overlap['VDJ_chain2_patient']==sample_trb_overlap['VDJ_chain2_product'])]

            # VDJ_chain2_x == VDJ_chain2_y
            sample_trb_overlap2_chain2 = sample_trb_overlap_chain2_nan.loc[(sample_trb_overlap_chain2_nan['VDJ_chain2_patient'].isnull()) | (sample_trb_overlap_chain2_nan['VDJ_chain2_product'].isnull()) | (sample_trb_overlap_chain2_nan['VDJ_chain2_patient']==sample_trb_overlap_chain2_nan['VDJ_chain2_product'])]

            #       Check if VDJ_cdr3s_2_x == VDJ_cdr3s_2_y.
            sample_trb_overlap2_seq2 = sample_trb_overlap2_chain2.loc[(sample_trb_overlap2_chain2['VDJ_cdr3s_2_patient'].isnull()) | (sample_trb_overlap2_chain2['VDJ_cdr3s_2_product'].isnull()) | (sample_trb_overlap2_chain2['VDJ_cdr3s_2_patient']==sample_trb_overlap2_chain2['VDJ_cdr3s_2_product'])]

            #Now check if VDJ_chain3_x is NaN? continue. if VDJ_chain3_y is NaN? continue. elif,
            sample_trb_overlap2_chain3_nan = sample_trb_overlap2_seq2.loc[(sample_trb_overlap2_seq2['VDJ_chain3_patient'].isnull()) | (sample_trb_overlap2_seq2['VDJ_chain3_product'].isnull()) | (sample_trb_overlap2_seq2['VDJ_chain3_patient']==sample_trb_overlap2_seq2['VDJ_chain3_product'])]

            #   VDJ_chain3_x == VDJ_chain3_y
            sample_trb_overlap2_chain3 = sample_trb_overlap2_chain3_nan.loc[(sample_trb_overlap2_chain3_nan['VDJ_chain3_patient'].isnull()) |sample_trb_overlap2_chain3_nan['VDJ_chain3_product'].isnull() | (sample_trb_overlap2_chain3_nan['VDJ_chain3_patient']==sample_trb_overlap2_chain3_nan['VDJ_chain3_product'])]

            #       Check if VDJ_cdr3s_3_x == VDJ_cdr3s_3_y.
            sample_trb_overlap2_seq3 = sample_trb_overlap2_chain3.loc[(sample_trb_overlap2_chain3['VDJ_cdr3s_3_patient'].isnull()) |sample_trb_overlap2_chain3['VDJ_cdr3s_3_product'].isnull() | (sample_trb_overlap2_chain3['VDJ_cdr3s_3_patient']==sample_trb_overlap2_chain3['VDJ_cdr3s_3_product'])]

            #Now check if VDJ_chain4_x is Nan? Continue. if VDJ_chain4_y is NaN? continue. elif,
            sample_trb_overlap2_chain4_nan = sample_trb_overlap2_seq3.loc[(sample_trb_overlap2_seq3['VDJ_chain4_patient'].isnull()) | (sample_trb_overlap2_seq3['VDJ_chain4_product'].isnull()) | (sample_trb_overlap2_seq3['VDJ_chain4_patient']==sample_trb_overlap2_seq3['VDJ_chain4_product'])]

            #   VDJ_chain4_x == VDJ_chain4_y
            sample_trb_overlap2_chain4 = sample_trb_overlap2_chain4_nan.loc[(sample_trb_overlap2_chain4_nan['VDJ_chain4_patient'].isnull())|sample_trb_overlap2_chain4_nan['VDJ_chain4_product'].isnull() | (sample_trb_overlap2_chain4_nan['VDJ_chain4_patient']==sample_trb_overlap2_chain4_nan['VDJ_chain4_product'])]

            #       Check if VDJ_cdr3s_4_x == VDJ_cdr3s_4_y.
            sample_trb_overlap2_seq4 = sample_trb_overlap2_chain4.loc[(sample_trb_overlap2_chain4['VDJ_cdr3s_4_patient'].isnull()) | (sample_trb_overlap2_chain4['VDJ_cdr3s_4_product'].isnull()) | (sample_trb_overlap2_chain4['VDJ_cdr3s_4_patient']==sample_trb_overlap2_chain4['VDJ_cdr3s_4_product'])]
            
            
            return sample_trb_overlap2_seq4 # These patient clones overlap with the product

        day0_product_associated = additional_overlap_checks(day0_product_overlap)
        post1_product_associated = additional_overlap_checks(df2_product_overlap)
        post2_product_associated = additional_overlap_checks(df3_product_overlap)

        # looking for the clonotype_id_y (which is the product clone id)-- if this shows up we want to remove the corresponding record
        clones_to_drop = day0_product_associated['clonotype_id_product']  #_y from _product
        def remove_day0prod_clones(post_product_associated, clones_to_drop):
            # if post_product_associated['clonotype_id_y] is in day0_product_associated['clonotype_id_y'] remove the entire row
            mask = ~post_product_associated['clonotype_id_product'].isin(clones_to_drop)
            # Filter the DataFrame using the inverted mask
            df_filtered = post_product_associated[mask]
            return df_filtered
        acquired_clones_post1 = remove_day0prod_clones(post1_product_associated,clones_to_drop)
        acquired_clones_post2 = remove_day0prod_clones(post2_product_associated,clones_to_drop)
        
        acquired_tracked_clones = acquired_clones_post1.merge(acquired_clones_post2, how='inner',on=['VDJ_cdr3s_1'],suffixes=('_firstPost','_secondPost')) # Could be NT too!!
        # Check if additional chains and sequences match
        #Now check if VDJ_chain2_x is NaN? continue. if VDJ_chain2_y is NaN? continue. elif,
        acquired_tracked_clones_chain2_nan= acquired_tracked_clones.loc[(acquired_tracked_clones['VDJ_chain2_patient_firstPost'].isnull()) | (acquired_tracked_clones['VDJ_chain2_patient_secondPost'].isnull()) | (acquired_tracked_clones['VDJ_chain2_patient_firstPost']==acquired_tracked_clones['VDJ_chain2_patient_secondPost'])]

        # VDJ_chain2_x == VDJ_chain2_y
        acquired_tracked_clones_chain2 = acquired_tracked_clones_chain2_nan.loc[(acquired_tracked_clones_chain2_nan['VDJ_chain2_patient_firstPost'].isnull()) | (acquired_tracked_clones_chain2_nan['VDJ_chain2_patient_secondPost'].isnull()) | (acquired_tracked_clones_chain2_nan['VDJ_chain2_patient_firstPost']==acquired_tracked_clones_chain2_nan['VDJ_chain2_patient_secondPost'])]

        #       Check if VDJ_cdr3s_2_x == VDJ_cdr3s_2_y.
        acquired_tracked_clones_seq2 = acquired_tracked_clones_chain2.loc[(acquired_tracked_clones_chain2['VDJ_cdr3s_2_patient_firstPost'].isnull()) | (acquired_tracked_clones_chain2['VDJ_cdr3s_2_patient_secondPost'].isnull()) | (acquired_tracked_clones_chain2['VDJ_cdr3s_2_patient_firstPost']==acquired_tracked_clones_chain2['VDJ_cdr3s_2_patient_secondPost'])]

        #Now check if VDJ_chain3_x is NaN? continue. if VDJ_chain3_y is NaN? continue. elif,
        acquired_tracked_clones_chain3_nan = acquired_tracked_clones_seq2.loc[(acquired_tracked_clones_seq2['VDJ_chain3_patient_firstPost'].isnull()) | (acquired_tracked_clones_seq2['VDJ_chain3_patient_secondPost'].isnull()) | (acquired_tracked_clones_seq2['VDJ_chain3_patient_firstPost']==acquired_tracked_clones_seq2['VDJ_chain3_patient_secondPost'])]

        #   VDJ_chain3_x == VDJ_chain3_y
        acquired_tracked_clones_chain3 = acquired_tracked_clones_chain3_nan.loc[(acquired_tracked_clones_chain3_nan['VDJ_chain3_patient_firstPost'].isnull()) |acquired_tracked_clones_chain3_nan['VDJ_chain3_patient_secondPost'].isnull() | (acquired_tracked_clones_chain3_nan['VDJ_chain3_patient_firstPost']==acquired_tracked_clones_chain3_nan['VDJ_chain3_patient_secondPost'])]

        #       Check if VDJ_cdr3s_3_x == VDJ_cdr3s_3_y.
        acquired_tracked_clones_seq3 = acquired_tracked_clones_chain3.loc[(acquired_tracked_clones_chain3['VDJ_cdr3s_3_patient_firstPost'].isnull()) |acquired_tracked_clones_chain3['VDJ_cdr3s_3_patient_secondPost'].isnull() | (acquired_tracked_clones_chain3['VDJ_cdr3s_3_patient_firstPost']==acquired_tracked_clones_chain3['VDJ_cdr3s_3_patient_secondPost'])]

        #Now check if VDJ_chain4_x is Nan? Continue. if VDJ_chain4_y is NaN? continue. elif,
        acquired_tracked_clones_chain4_nan = acquired_tracked_clones_seq3.loc[(acquired_tracked_clones_seq3['VDJ_chain4_patient_firstPost'].isnull()) | (acquired_tracked_clones_seq3['VDJ_chain4_patient_secondPost'].isnull()) | (acquired_tracked_clones_seq3['VDJ_chain4_patient_firstPost']==acquired_tracked_clones_seq3['VDJ_chain4_patient_secondPost'])]

        #   VDJ_chain4_x == VDJ_chain4_y
        acquired_tracked_clones_chain4 = acquired_tracked_clones_chain4_nan.loc[(acquired_tracked_clones_chain4_nan['VDJ_chain4_patient_firstPost'].isnull())|acquired_tracked_clones_chain4_nan['VDJ_chain4_patient_secondPost'].isnull() | (acquired_tracked_clones_chain4_nan['VDJ_chain4_patient_firstPost']==acquired_tracked_clones_chain4_nan['VDJ_chain4_patient_secondPost'])]

        #       Check if VDJ_cdr3s_4_x == VDJ_cdr3s_4_y.
        acquired_tracked_clones_seq4 = acquired_tracked_clones_chain4.loc[(acquired_tracked_clones_chain4['VDJ_cdr3s_4_patient_firstPost'].isnull()) | (acquired_tracked_clones_chain4['VDJ_cdr3s_4_patient_secondPost'].isnull()) | (acquired_tracked_clones_chain4['VDJ_cdr3s_4_patient_firstPost']==acquired_tracked_clones_chain4['VDJ_cdr3s_4_patient_secondPost'])]
        # add a column for day 0 being 0
        acquired_tracked_clones_seq4['proportion_day0'] = 0
        
        #return render.DataTable(acquired_tracked_clones_seq4)
        acquired_tracked_clones_seq4_df = pd.DataFrame(acquired_tracked_clones_seq4)
        return acquired_tracked_clones_seq4_df
    
    @output
    @render.data_frame
    def table_results():
        return clone_tracking()
    
    @render.plot
    def plot_results():
        acquired_tracked_clones_seq4_df = clone_tracking()
        acquired_tracked_clones_seq4_df = pd.DataFrame(acquired_tracked_clones_seq4_df)
        # Get the timepoints
        timepoint1 = input.timepoint1()
        timepoint2 = input.timepoint2()
        
        # Generate random colors
        # Function to generate a list of N random RGB colors
        def generate_random_colors(n_colors):
            colors = [(random.random(), random.random(), random.random()) for _ in range(n_colors)]
            return colors
        num_categories = len(acquired_tracked_clones_seq4_df['clonotype_id_patient_firstPost'].unique()) # edit -- clonotype_id_x
        if num_categories == 0:
            fig=plt.figure()
            ax1=fig.add_subplot()
            p = plt.stackplot(1,1) #("")
            ax1.set_title('No tracked clones identified.')
        elif num_categories > 0:
            random_palette = generate_random_colors(num_categories)
            # Set this as the current palette
            colorsagain2 = sns.set_palette(random_palette)
            y = acquired_tracked_clones_seq4_df[['proportion_day0','proportion_patient_firstPost','proportion_patient_secondPost']]
            x = (0,timepoint1,timepoint2) # User input for timepoints
            ys = []
            ys = y.values.tolist()
            fig=plt.figure()
            ax1=fig.add_subplot()
            p = plt.stackplot(x,*ys,colors=colorsagain2)
            ax1.set_xlabel('Days Post Infusion')
            ax1.set_ylabel('Clonotype Frequency')
            ax1.set_title('Acquired Clonotype Frequency')
        #ax1.set_ylim(top=0.085) # Could add a user input for limits on the y-axis range if we want
        return p
        
    # Create download button
    @render.download(filename=lambda:"acquired_clonotype_tracking.csv")
    async def download_data(): # async tested -- didn't use await so it failed. 
        await asyncio.sleep(0.25)
        #yield "one,two,three\n"
        acquired_tracked_clones_seq4_df = clone_tracking()
        acquired_tracked_clones_seq4_df = pd.DataFrame(acquired_tracked_clones_seq4_df)
        #path = acquired_tracked_clones_seq4_df.to_csv("acquired_clonotype_tracking.csv")# Save the DataFrame to a CSV file
        # Return dataframe to user 
        #return path
        yield acquired_tracked_clones_seq4_df.to_csv(index=False) # Return the CSV data as a string for download, without the index column


app = App(app_ui, server)