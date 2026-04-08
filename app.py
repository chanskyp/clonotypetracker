import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import random 


from shiny import App, render, ui

# User interface (UI) definition
app_ui = ui.page_fluid(
    # Add a title to the page with some top padding
    ui.panel_title(ui.h2("Acquired Clonotype Tracking Plotter", class_="pt-5"),
                   #ui.h4("Acquired Clonotype Tracking Plotter", class_="text-muted")
                   ),
    ui.tags.p("Upload your product and patient sample data to visualize the frequency of acquired clonotypes over time post-infusion. Samples should be uploaded in the following order: product, day of infusion, early post-infusion, late post-infusion."),
    #
    # User input files
    ui.input_file(id="file_upload", label="Choose CSV Files", accept=[".csv"], multiple=True),

    # A numeric input for the number of days post-infusion the first sample was taken
    ui.input_numeric("timepoint1", "Days Post Infusion for Sample 1 (e.g.,35)", value=35),
    
    # A second numeric input for the number of days post-infusion the second sample was taken
    ui.input_numeric("timepoint2", "Days Post Infusion for Sample 2 (e.g.,180)", value=180),
    
    # A container for plot output
    ui.output_plot("contents"),
)

# Using AA seq-- update to a user option for AA or NT
def server(input,output,session): #, output, session
    #@render.table
    @output
    @render.plot
    def contents():
        # input.file1() returns a list of dictionaries or None
        file_info = input.file_upload()
        if file_info is None:
            return None
        
        # Get the timepoints
        timepoint1 = input.timepoint1()
        timepoint2 = input.timepoint2()
        # Access the first uploaded file's temporary path
        # 2. Read the CSV using the 'datapath'
        # file_info is a list of dicts, get the first one
    
        df_product= pd.read_csv(file_info[0]["datapath"])
        df1= pd.read_csv(file_info[1]["datapath"])
        df2= pd.read_csv(file_info[2]["datapath"])
        df3= pd.read_csv(file_info[3]["datapath"])

        # Here we only look at the TCR-beta chain. Future iteration to allow user to select TRA vs TRB vs combo
        def trb_only(df):
            df = df[df.VDJ_chain1 == 'TRB']
            return df
        
        df_product = trb_only(df_product)
        df1 = trb_only(df1)
        df2 = trb_only(df2)
        df3 = trb_only(df3)

        # Inner join the product and patient samples
        # We join on CDR3 AA seq-- future iteration to allow user to choose either NT or AA sequence
        day0_product_overlap = df1.merge(df_product, how='inner',on=['VDJ_cdr3s_aa1']) #VDJ_cdr3s_nt1
        df2_product_overlap = df2.merge(df_product, how='inner',on=['VDJ_cdr3s_aa1']) # Early post
        df3_product_overlap = df3.merge(df_product, how='inner',on=['VDJ_cdr3s_aa1']) # Late post 
        
        def additional_overlap_checks(sample_trb_overlap): #use any of the prior dataframes as the input
            # Check if additional chains and sequences match
            #Now check if VDJ_chain2_x is NaN? continue. if VDJ_chain2_y is NaN? continue. elif,
            sample_trb_overlap_chain2_nan= sample_trb_overlap.loc[(sample_trb_overlap['VDJ_chain2_x'].isnull()) | (sample_trb_overlap['VDJ_chain2_y'].isnull()) | (sample_trb_overlap['VDJ_chain2_x']==sample_trb_overlap['VDJ_chain2_y'])]

            # VDJ_chain2_x == VDJ_chain2_y
            sample_trb_overlap2_chain2 = sample_trb_overlap_chain2_nan.loc[(sample_trb_overlap_chain2_nan['VDJ_chain2_x'].isnull()) | (sample_trb_overlap_chain2_nan['VDJ_chain2_y'].isnull()) | (sample_trb_overlap_chain2_nan['VDJ_chain2_x']==sample_trb_overlap_chain2_nan['VDJ_chain2_y'])]

            #       Check if VDJ_cdr3s_nt2_x == VDJ_cdr3s_nt2_y.
            sample_trb_overlap2_seq2 = sample_trb_overlap2_chain2.loc[(sample_trb_overlap2_chain2['VDJ_cdr3s_nt2_x'].isnull()) | (sample_trb_overlap2_chain2['VDJ_cdr3s_nt2_y'].isnull()) | (sample_trb_overlap2_chain2['VDJ_cdr3s_nt2_x']==sample_trb_overlap2_chain2['VDJ_cdr3s_nt2_y'])]

            #Now check if VDJ_chain3_x is NaN? continue. if VDJ_chain3_y is NaN? continue. elif,
            sample_trb_overlap2_chain3_nan = sample_trb_overlap2_seq2.loc[(sample_trb_overlap2_seq2['VDJ_chain3_x'].isnull()) | (sample_trb_overlap2_seq2['VDJ_chain3_y'].isnull()) | (sample_trb_overlap2_seq2['VDJ_chain3_x']==sample_trb_overlap2_seq2['VDJ_chain3_y'])]

            #   VDJ_chain3_x == VDJ_chain3_y
            sample_trb_overlap2_chain3 = sample_trb_overlap2_chain3_nan.loc[(sample_trb_overlap2_chain3_nan['VDJ_chain3_x'].isnull()) |sample_trb_overlap2_chain3_nan['VDJ_chain3_y'].isnull() | (sample_trb_overlap2_chain3_nan['VDJ_chain3_x']==sample_trb_overlap2_chain3_nan['VDJ_chain3_y'])]

            #       Check if VDJ_cdr3s_nt3_x == VDJ_cdr3s_nt3_y.
            sample_trb_overlap2_seq3 = sample_trb_overlap2_chain3.loc[(sample_trb_overlap2_chain3['VDJ_cdr3s_nt3_x'].isnull()) |sample_trb_overlap2_chain3['VDJ_cdr3s_nt3_y'].isnull() | (sample_trb_overlap2_chain3['VDJ_cdr3s_nt3_x']==sample_trb_overlap2_chain3['VDJ_cdr3s_nt3_y'])]

            #Now check if VDJ_chain4_x is Nan? Continue. if VDJ_chain4_y is NaN? continue. elif,
            sample_trb_overlap2_chain4_nan = sample_trb_overlap2_seq3.loc[(sample_trb_overlap2_seq3['VDJ_chain4_x'].isnull()) | (sample_trb_overlap2_seq3['VDJ_chain4_y'].isnull()) | (sample_trb_overlap2_seq3['VDJ_chain4_x']==sample_trb_overlap2_seq3['VDJ_chain4_y'])]

            #   VDJ_chain4_x == VDJ_chain4_y
            sample_trb_overlap2_chain4 = sample_trb_overlap2_chain4_nan.loc[(sample_trb_overlap2_chain4_nan['VDJ_chain4_x'].isnull())|sample_trb_overlap2_chain4_nan['VDJ_chain4_y'].isnull() | (sample_trb_overlap2_chain4_nan['VDJ_chain4_x']==sample_trb_overlap2_chain4_nan['VDJ_chain4_y'])]

            #       Check if VDJ_cdr3s_nt4_x == VDJ_cdr3s_nt4_y.
            sample_trb_overlap2_seq4 = sample_trb_overlap2_chain4.loc[(sample_trb_overlap2_chain4['VDJ_cdr3s_nt4_x'].isnull()) | (sample_trb_overlap2_chain4['VDJ_cdr3s_nt4_y'].isnull()) | (sample_trb_overlap2_chain4['VDJ_cdr3s_nt4_x']==sample_trb_overlap2_chain4['VDJ_cdr3s_nt4_y'])]
            
            
            return sample_trb_overlap2_seq4 # These patient clones overlap with the product

        day0_product_associated = additional_overlap_checks(day0_product_overlap)
        post1_product_associated = additional_overlap_checks(df2_product_overlap)
        post2_product_associated = additional_overlap_checks(df3_product_overlap)

        # looking for the clonotype_id_y (which is the product clone id)-- if this shows up we want to remove the corresponding record
        clones_to_drop = day0_product_associated['clonotype_id_y']
        def remove_day0prod_clones(post_product_associated, clones_to_drop):
            # if post_product_associated['clonotype_id_y] is in day0_product_associated['clonotype_id_y'] remove the entire row
            mask = ~post_product_associated['clonotype_id_y'].isin(clones_to_drop)
            # Filter the DataFrame using the inverted mask
            df_filtered = post_product_associated[mask]
            return df_filtered
        acquired_clones_post1 = remove_day0prod_clones(post1_product_associated,clones_to_drop)
        acquired_clones_post2 = remove_day0prod_clones(post2_product_associated,clones_to_drop)
        
        acquired_tracked_clones = acquired_clones_post1.merge(acquired_clones_post2, how='inner',on=['VDJ_cdr3s_aa1']) # Could be NT too!!
        # Check if additional chains and sequences match
        #Now check if VDJ_chain2_x is NaN? continue. if VDJ_chain2_y is NaN? continue. elif,
        acquired_tracked_clones_chain2_nan= acquired_tracked_clones.loc[(acquired_tracked_clones['VDJ_chain2_x_x'].isnull()) | (acquired_tracked_clones['VDJ_chain2_x_y'].isnull()) | (acquired_tracked_clones['VDJ_chain2_x_x']==acquired_tracked_clones['VDJ_chain2_x_y'])]

        # VDJ_chain2_x == VDJ_chain2_y
        acquired_tracked_clones_chain2 = acquired_tracked_clones_chain2_nan.loc[(acquired_tracked_clones_chain2_nan['VDJ_chain2_x_x'].isnull()) | (acquired_tracked_clones_chain2_nan['VDJ_chain2_x_y'].isnull()) | (acquired_tracked_clones_chain2_nan['VDJ_chain2_x_x']==acquired_tracked_clones_chain2_nan['VDJ_chain2_x_y'])]

        #       Check if VDJ_cdr3s_aa2_x == VDJ_cdr3s_aa2_y.
        acquired_tracked_clones_seq2 = acquired_tracked_clones_chain2.loc[(acquired_tracked_clones_chain2['VDJ_cdr3s_aa2_x_x'].isnull()) | (acquired_tracked_clones_chain2['VDJ_cdr3s_aa2_x_y'].isnull()) | (acquired_tracked_clones_chain2['VDJ_cdr3s_aa2_x_x']==acquired_tracked_clones_chain2['VDJ_cdr3s_aa2_x_y'])]

        #Now check if VDJ_chain3_x is NaN? continue. if VDJ_chain3_y is NaN? continue. elif,
        acquired_tracked_clones_chain3_nan = acquired_tracked_clones_seq2.loc[(acquired_tracked_clones_seq2['VDJ_chain3_x_x'].isnull()) | (acquired_tracked_clones_seq2['VDJ_chain3_x_y'].isnull()) | (acquired_tracked_clones_seq2['VDJ_chain3_x_x']==acquired_tracked_clones_seq2['VDJ_chain3_x_y'])]

        #   VDJ_chain3_x == VDJ_chain3_y
        acquired_tracked_clones_chain3 = acquired_tracked_clones_chain3_nan.loc[(acquired_tracked_clones_chain3_nan['VDJ_chain3_x_x'].isnull()) |acquired_tracked_clones_chain3_nan['VDJ_chain3_x_y'].isnull() | (acquired_tracked_clones_chain3_nan['VDJ_chain3_x_x']==acquired_tracked_clones_chain3_nan['VDJ_chain3_x_y'])]

        #       Check if VDJ_cdr3s_aa3_x == VDJ_cdr3s_aa3_y.
        acquired_tracked_clones_seq3 = acquired_tracked_clones_chain3.loc[(acquired_tracked_clones_chain3['VDJ_cdr3s_aa3_x_x'].isnull()) |acquired_tracked_clones_chain3['VDJ_cdr3s_aa3_x_y'].isnull() | (acquired_tracked_clones_chain3['VDJ_cdr3s_aa3_x_x']==acquired_tracked_clones_chain3['VDJ_cdr3s_aa3_x_y'])]

        #Now check if VDJ_chain4_x is Nan? Continue. if VDJ_chain4_y is NaN? continue. elif,
        acquired_tracked_clones_chain4_nan = acquired_tracked_clones_seq3.loc[(acquired_tracked_clones_seq3['VDJ_chain4_x_x'].isnull()) | (acquired_tracked_clones_seq3['VDJ_chain4_x_y'].isnull()) | (acquired_tracked_clones_seq3['VDJ_chain4_x_x']==acquired_tracked_clones_seq3['VDJ_chain4_x_y'])]

        #   VDJ_chain4_x == VDJ_chain4_y
        acquired_tracked_clones_chain4 = acquired_tracked_clones_chain4_nan.loc[(acquired_tracked_clones_chain4_nan['VDJ_chain4_x_x'].isnull())|acquired_tracked_clones_chain4_nan['VDJ_chain4_x_y'].isnull() | (acquired_tracked_clones_chain4_nan['VDJ_chain4_x_x']==acquired_tracked_clones_chain4_nan['VDJ_chain4_x_y'])]

        #       Check if VDJ_cdr3s_aa4_x == VDJ_cdr3s_aa4_y.
        acquired_tracked_clones_seq4 = acquired_tracked_clones_chain4.loc[(acquired_tracked_clones_chain4['VDJ_cdr3s_aa4_x_x'].isnull()) | (acquired_tracked_clones_chain4['VDJ_cdr3s_aa4_x_y'].isnull()) | (acquired_tracked_clones_chain4['VDJ_cdr3s_aa4_x_x']==acquired_tracked_clones_chain4['VDJ_cdr3s_aa4_x_y'])]
        # add a column for day 0 being 0
        acquired_tracked_clones_seq4['proportion_day0'] = 0
        # Generate random colors
        # Function to generate a list of N random RGB colors
        def generate_random_colors(n_colors):
            colors = [(random.random(), random.random(), random.random()) for _ in range(n_colors)]
            return colors
        num_categories = len(acquired_tracked_clones_seq4['clonotype_id_x_x'].unique()) # edit -- clonotype_id_x
        random_palette = generate_random_colors(num_categories)
        # Set this as the current palette
        colorsagain2 = sns.set_palette(random_palette)
        y = acquired_tracked_clones_seq4[['proportion_day0','proportion_x_x','proportion_x_y']]
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



app = App(app_ui, server)


