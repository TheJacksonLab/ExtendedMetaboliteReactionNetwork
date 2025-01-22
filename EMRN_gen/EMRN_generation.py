import pandas as pd
import EMRN_gen_fns as eg

id_reac = pd.read_csv('id_reac_final.csv')
met = pd.read_csv('met_all_noMW_common.csv')

# Extract USPTO reactant-reaction id dataframe and canonicalized series of reactants
unique_id_reac_df, reac_series_can = eg.extract_uspto_reactants(id_reac)

# Run five rounds of reactions by matching USPTO reactants with complete set of potential reactants for each round:
# which are the set of metabolites for round 1, and the complete set of product molecules from round n-1 for round n for the subsequent rounds
# This is followed by filtering of the reactions based on the molecular weight proxy

round1 = eg.emrn_reactant_matching(met['canonical smiles'], reac_series_can, unique_id_reac_df, id_reac)
sus_round1 = eg.mw_sus_filtering(round1, id_reac)

round2 = eg.emrn_reactant_matching(sus_round1['Product'], reac_series_can, unique_id_reac_df, id_reac)
sus_round2 = eg.mw_sus_filtering(round2, id_reac)

round3 = eg.emrn_reactant_matching(sus_round2['Product'], reac_series_can, unique_id_reac_df, id_reac)
sus_round3 = eg.mw_sus_filtering(round3, id_reac)

round4 = eg.emrn_reactant_matching(sus_round3['Product'], reac_series_can, unique_id_reac_df, id_reac)
sus_round4 = eg.mw_sus_filtering(round4, id_reac)

round5 = eg.emrn_reactant_matching(sus_round4['Product'], reac_series_can, unique_id_reac_df, id_reac)
sus_round5 = eg.mw_sus_filtering(round5, id_reac)

# Combine reaction data into a dataframe with each reaction tagged with the round they belong to, without repetition
# For eg. if a given reaction is present in both round 1 and round 4, it will be tagged with round 1 only, since that comes first

# Add 'Round' column to tag each dataframe
sus_round1['Round'] = 'Round1'
sus_round2['Round'] = 'Round2'
sus_round3['Round'] = 'Round3'
sus_round4['Round'] = 'Round4'
sus_round5['Round'] = 'Round5'

# Combine the first dataframe
combined_df = sus_round1.copy()

# Add each dataframe while checking for duplicates
combined_df = eg.add_non_duplicates('reaction', combined_df, sus_round2, 'Round2')
combined_df = eg.add_non_duplicates('reaction', combined_df, sus_round3, 'Round3')
combined_df = eg.add_non_duplicates('reaction', combined_df, sus_round4, 'Round4')
combined_df = eg.add_non_duplicates('reaction', combined_df, sus_round5, 'Round5')

#####################################################################################################################
# Dataframe to store all the molecules generated in the EMRN, tagged by which round they make their first appearance in

#Create a dataframe from metabolites and tag them as 'Metabolite'
met_df = pd.DataFrame(met['canonical smiles'].unique(), columns=['smiles'])
met_df['tag'] = 'Metabolite'

# Create dataframes for each round's unique targets and tag them accordingly
round1_df = pd.DataFrame(sus_round1['Product'].unique(), columns=['smiles'])
round1_df['tag'] = 'Round1'

round2_df = pd.DataFrame(sus_round2['Product'].unique(), columns=['smiles'])
round2_df['tag'] = 'Round2'

round3_df = pd.DataFrame(sus_round3['Product'].unique(), columns=['smiles'])
round3_df['tag'] = 'Round3'

round4_df = pd.DataFrame(sus_round4['Product'].unique(), columns=['smiles'])
round4_df['tag'] = 'Round4'

round5_df = pd.DataFrame(sus_round5['Product'].unique(), columns=['smiles'])
round5_df['tag'] = 'Round5'

combined_df2 = met_df.copy()  # Start with metabolites
combined_df2 = eg.add_non_duplicates('nonreaction', combined_df2, round1_df, 'Round1')
combined_df2 = eg.add_non_duplicates('nonreaction', combined_df2, round2_df, 'Round2')
combined_df2 = eg.add_non_duplicates('nonreaction', combined_df2, round3_df, 'Round3')
combined_df2 = eg.add_non_duplicates('nonreaction', combined_df2, round4_df, 'Round4')
combined_df2 = eg.add_non_duplicates('nonreaction', combined_df2, round5_df, 'Round5')


combined_df2.to_csv('metr_sus_df.csv')
combined_df.to_csv('sus_rxn5.csv')

# To access a given reaction index in the original USPTO database based on Rxn_idx in the combined_df:

# uspto = pd.read_csv(USPTO_FULL.csv)  # Original USPTO grants database
# uspto.iloc[rxn_idx]  # rxn_idx is the reaction you want to look at