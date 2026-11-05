# Summarizes key edits to the iIsor model for RBA modeling (not including added biomass rxns)
# See "#TODO" comments for potential future edits to conside/reconsider
import cobra,json,os,sys,re

model_root_path = '../../../'
sys.path.append(model_root_path)
from pycore import *

# change gene names to avoid RBA names being too long
for gene in 'JL09_g1290','JL09_g1433','JL09_g1518','JL09_g1864','JL09_g223','JL09_g2627','JL09_g265','JL09_g2996','JL09_g3290','JL09_g3484','JL09_g35','JL09_g3504','JL09_g3584','JL09_g3947','JL09_g4008','JL09_g4308','JL09_g4546','JL09_g5116','JL09_g5268','JL09_g5874','JL09_g6026','JL09_g641','JL09_g755','JL09_g800','JL09_g971','JL09_g3308':
    rename_gene(model, gene, gene.split('_')[1])

# remove duplicate genes that cause issues with RBA model building (seem to be alleles of other existing genes)
genes_to_remove = ['JL09_g2650','JL09_g4804','JL09_g4510','JL09_g4745','JL09_g982', 'JL09_g6147']
for gene in genes_to_remove:
    if gene in model.genes:
        cobra.manipulation.remove_genes(model, [gene])

# add cobalt2 to e and c compartments
for compartment in ['e','c']:
    met_id = f'cobalt2_{compartment}'
    if met_id not in model.metabolites:
        cobalt2 = cobra.Metabolite(
            id=met_id,
            name='Cobalt(II) ion',
            formula='Co1',
            charge=2,
            compartment=compartment
        )
        cobalt2.annotation = {
            "biocyc":"META:CO+2",
            'chebi': 'CHEBI:48828',
            "kegg.compound":"C00175",
            "pubchem.compound":"104729",
            "sbo":"SBO:0000247",
            "seed.compound":"cpd00149"
        }
        model.add_metabolites([cobalt2])
    # keep ni2 section in if needed later; but for now, nickel is not included in the model due to lack of evidence and lack of necessity for RBA modeling (no nickel-requiring enzymes identified)
    # met_id = f'ni2_{compartment}'
    # if met_id not in model.metabolites:
    #     ni2 = cobra.Metabolite(
    #         id=met_id,
    #         name='Nickel(II) ion',
    #         formula='Ni1',
    #         charge=2,
    #         compartment=compartment
    #     )
    #     ni2.annotation = {
    #         "biocyc":"META:NI+2",
    #         'chebi': 'CHEBI:49786',
    #         "kegg.compound":"C00244",
    #         "pubchem.compound":"934",
    #         "sbo":"SBO:0000247",
    #         "seed.compound":"cpd20862"
    #     }
    #     model.add_metabolites([ni2])

# add drafts of genes from mtDNA
mtDNA_genes = {'ATP6':'A0A2U9RAK4','ATP8':'A0A2U9RAV5','ATP9':'A0A2U9RBW3','COB':'A0A2U9RBK5','COX1':'A0A2U9RCF9','COX2':'A0A2U9RAR9','COX3':'A0A2U9RBJ3','NAD1':'A0A2U9RCC6','NAD2':'A0A2U9RB22','NAD4':'A0A2U9RCH4','NAD4L':'A0A2U9RAU9','NAD5':'A0A2U9RAN3','NAD6':'A0A2U9RB35'}
for gene_id, uniprot_id in mtDNA_genes.items():
    if "iom" + gene_id not in model.genes:
        new_gene = cobra.Gene("iom" + gene_id, name=gene_id)
        new_gene.notes['notes'] = 'Draft gene added from mtDNA; annotation may be incomplete. UniProt ID is from homologous gene in strain CBS 573.'
        new_gene.annotation['uniprot'] = uniprot_id
        model.genes.append(new_gene)
        # find genes with gene_id in their names and add this new gene to their GPRs as an "or" option
        alt_genes = [gene for gene in model.genes if re.search(r'\b' + re.escape(gene_id) + r'\b', gene.name) and gene.id != new_gene.id]
        if not alt_genes:
            print(f"Warning: no gene counterparts found for gene {gene_id}")
        for alt_gene in alt_genes:
            for rxn in alt_gene.reactions:
                rxn.gene_reaction_rule = re.sub(r'\b' + re.escape(alt_gene.id) + r'\b', f"({alt_gene.id} or {new_gene.id})", rxn.gene_reaction_rule)

# add other genes
genes_to_add = [
{
"id":"JL09_g3719",
"name":"PPA2",
"annotation":{
"sc_id":"YMR267W",
"ncbiprotein":"KGK37152.1",
"sbo":"SBO:0000243",
"uniprot":"A0A099NZ07"
}
},
{
"id":"JL09_g420",
"name":"ALD4",
"annotation":{
"ncbiprotein":"KGK40528.1",
"sbo":"SBO:0000243",
"sc_id":"YOR374W",
"uniprot":"A0A099P647"
}
},
{
"id":"JL09_g946",
"name":"ADE3",
"annotation":{
"ncbiprotein":"KGK39904.1",
"sbo":"SBO:0000243",
"uniprot":"A0A099P6Q8"
}
},
{
"id":"io248592",
"name":"HMG aldolase",
"annotation":{
"sbo":"SBO:0000243",
"sc_id":"YER010C",
"note":"gene not found on UniProt as of 2026-03-03; uniprot value is for S. cerevisiae homolog",
"uniprot":"P40011"
}
},
{
"id":"JL09_g3775",
"name":"ATO2,ADY2",
"notes":{
"notes":"Added manually based on atom mapping model GPRs. Its protein corresponds to protein 303834 in Issorie2 on JGI."
},
"annotation":{
"ncbiprotein":"KGK37087",
"sbo":"SBO:0000243",
"sc_id":"YNR002C,YCR010C",
"uniprot":"A0A099NYU0"
}
},
{
"id":"JL09_g191",
"name":"ATO2,ADY2",
"notes":{
"notes":"Added manually based on atom mapping model GPRs. Its protein corresponds to protein 317917 in Issorie2 on JGI."
},
"annotation":{
"ncbiprotein":"KGK40727",
"sbo":"SBO:0000243",
"sc_id":"YNR002C,YCR010C",
"uniprot":"A0A099P8V4"
}
},
{
"id":"JL09_g192",
"name":"ATO2,ADY2",
"notes":{
"notes":"Added manually based on atom mapping model GPRs. Its protein corresponds to protein 334551 in Issorie2 on JGI."
},
"annotation":{
"ncbiprotein":"KGK40728",
"sbo":"SBO:0000243",
"sc_id":"YNR002C,YCR010C",
"uniprot":"A0A099P6N8"
}
},
{
"id":"JL09_g2166",
"name":"ORT1",
"notes":{
"notes":"Added manually to address missing transporter in proteomics data."
},
"annotation":{
"sbo":"SBO:0000243",
"sc_id":"YOR130C",
"uniprot":"A0A099P163"
}
},
]
for gene_data in genes_to_add:
    if gene_data['id'] not in model.genes:
        new_gene = cobra.Gene(gene_data['id'], name=gene_data['name'])
        new_gene.annotation = gene_data['annotation']
        model.genes.append(new_gene)

model.reactions.ORNtpa_m.gene_reaction_rule = "JL09_g5501 or JL09_g2166"

## fill in missing proteins in GPRs using draft mtDNA genes
model.reactions.ATPS_m.gene_reaction_rule = model.reactions.ATPS_m.gene_reaction_rule.replace("UNKNOWN and UNKNOWN","iomATP8 and iomATP9")
## add COX1 and COX3 to FECOOR_m since no counterparts exist in the model? No, to avoid overestimating protein requirements for RBA (since these genes are likely needed but it's unclear where in the proteomics data their counterparts are)
# model.reactions.FECOOR_m.gene_reaction_rule = model.reactions.FECOOR_m.gene_reaction_rule + " and iomCOX1 and iomCOX3"

# JL09_g746 is listed on UniProt as the large chain of a carbamoyl phosphate synthase complex, for which it also lists that it's a "Heterodimer composed of 2 chains; the small (or glutamine) chain promotes the hydrolysis of glutamine to ammonia, which is used by the large (or ammonia) chain to synthesize carbamoyl phosphate". Thus, the small chain (which both JL09_g4693 and JL09_g4696 are described as being, on their respective UniProt pages) seems to be required.
model.reactions.CBPS_c.gene_reaction_rule = "JL09_g338 or ((JL09_g4693 or JL09_g4696) and JL09_g746)"
model.reactions.CBPSQ_m.gene_reaction_rule = "(JL09_g4693 or JL09_g4696) and JL09_g746"
model.reactions.GLUN_c.gene_reaction_rule = "JL09_g3834"
model.reactions.PAILS_rm.gene_reaction_rule = "JL09_g976"

# add the following genes to all GPRs where their other homolog is present, since they are both homologs of the same S. cerevisiae gene
for rxn in model.reactions:
    rxn_genes_ids = [gene.id for gene in rxn.genes]
    if "JL09_g3946" in rxn_genes_ids and "JL09_g3945" not in rxn_genes_ids:
        rxn.gene_reaction_rule = re.sub(r'\bJL09_g3946\b', "(JL09_g3946 or JL09_g3945)", rxn.gene_reaction_rule)
    if "JL09_g3402" in rxn_genes_ids and "JL09_g100" not in rxn_genes_ids:
        rxn.gene_reaction_rule = re.sub(r'\bJL09_g3402\b', "(JL09_g3402 or JL09_g100)", rxn.gene_reaction_rule)

model.reactions.EX_fe3_e.lower_bound = -1000

model.metabolites.pheme_c.formula = 'C34H30N4O4Fe' # fixed typo where F was included.
model.metabolites.cytc_c.formula = 'C34H30N4O4ApocytcFe' # fixed typo where F was included.

# add exchange reactions for cobalt2
for rxn in ('EX_cobalt2_e','EX_ni2_e'):
    met = rxn.split('EX_')[1]
    if met not in model.metabolites:
        continue
    if rxn not in model.reactions:
        model.add_boundary(metabolite=model.metabolites.get_by_id(rxn.split('EX_')[1]), type='exchange', lb=-1000, ub=1000)
# assume diffusion of metal ions into the mitochondria, since UniProt predicts mitochondrial localization for some of them
# assume AMP diffusion since it's allowed in other compartments
# assume ppa_c diffusion based on r_4686 in yeast-GEM v9.1.1
for met, compartments in {'ppa':{'m':'c'},'cu2':{'c':'mm','mm':'m'}, 'cobalt2':{'e':'c'}, 'ni2':{'e':'c'}, 'amp':{'c':'en'}, 'coa':{'c':'en'}, 'ppi':{'c':'en'}}.items():
    for from_comp, to_comp in compartments.items():
        met1_id = f'{met}_{from_comp}'
        if met1_id not in model.metabolites:
            continue
        met1 = model.metabolites.get_by_id(met1_id)
        # create metabolite in other compartment if it doesn't exist
        met2_id = f'{met}_{to_comp}'
        model = generate_metabolite_in_diff_compartment(model, met2_id, to_comp, None, met1.id)
        met2 = model.metabolites.get_by_id(met2_id)
        # create diffusion reaction if it doesn't exist
        rxn_id = f'{met.upper()}t_{to_comp}_{from_comp}'
        # check if reaction exists, based on stoichiometry
        rxn_exists = False
        for rxn in model.reactions:
            if rxn.metabolites == {met1:-1, met2:1} or rxn.metabolites == {met1:1, met2:-1}:
                rxn_exists = True
                break
        if not rxn_exists:
            diffusion_rxn = cobra.Reaction(rxn_id)
            diffusion_rxn.name = f'{met1.name} transport (or diffusion?)'
            diffusion_rxn.add_metabolites({met1:1, met2:-1})
            diffusion_rxn.lower_bound = -1000
            diffusion_rxn.upper_bound = 1000
            diffusion_rxn.annotation = {
                "sbo":"SBO:0000185"
            }
            if met in ('amp','coa'):
                #TODO: Reexamine these someday, as they're potentially no longer needed for addressing RBA parameterization
                diffusion_rxn.notes['notes'] = f"Reaction added to allow {met1.name} transport between compartments, since RBA parameterization requires them for FACOAL240_en (which is crucial for kapp calculations that fit the proteomics data accurately)."
            model.add_reactions([diffusion_rxn])

# add FMNH2 to the model, and fumarate reductase reaction that uses it
if 'fmnh2_c' not in model.metabolites:
    fmnh2_c = cobra.Metabolite(
        id='fmnh2_c',
        name='Reduced flavin mononucleotide',
        formula='C17H21N4O9P',
        charge=-2,
        compartment='c'
    )
    fmnh2_c.annotation = {
        "biocyc":"META:FMNH2",
        'chebi': 'CHEBI:57618',
        "kegg.compound":"C01847",
        "pubchem.compound":"44229161",
        "sbo":"SBO:0000247",
        "seed.compound":"cpd01270"
    }
    model.add_metabolites([fmnh2_c])
    # add reaction
    frdf_c = cobra.Reaction('FRDFMN_c', 
                            name = 'fumarate reductase (FMN)',
                            subsystem = 'Unassigned'
    )
    fumarate_c = model.metabolites.get_by_id('fum_c')
    succinate_c = model.metabolites.get_by_id('succ_c')
    fmn_c = model.metabolites.get_by_id('fmn_c')
    frdf_c.add_metabolites({
        fumarate_c: -1,
        fmnh2_c: -1,
        succinate_c: 1,
        fmn_c: 1,
        model.metabolites.get_by_id('h_c'): 1
    })
    frdf_c.lower_bound = -1000
    frdf_c.upper_bound = 1000
    frdf_c.gene_reaction_rule = 'JL09_g265'
    frdf_c.notes['notes'] = {
        'yeastkegg_id':'r_1000',
        'notes':'Reaction added to include FMNH2 in the model, for RBA modeling purposes. GPR based on BLAST search of YEL047C (S. cerevisiae fumarate reductase) against iIsor proteome, which identified JL09_g265 as the best hit.'
    }
    frdf_c.annotation = {
        "ec-code":"1.3.1.6",
        "kegg.reaction":"R00408",
        "metanetx.reaction":"MNXR129857",
        "sbo":"SBO:0000176"
    }
    model.add_reactions([frdf_c])
    # check if mass and charge balanced
    imbal = check_mass_balance_cobra(frdf_c,model)
    if any(imbal.values()):
        print(f"Warning: FRDFMN_c reaction is not mass/charge balanced: {imbal}")
    # test if model can make fmnh2_c when a demand reaction is added
    dm_fmnh2_c = cobra.Reaction('DM_fmnh2_c')
    dm_fmnh2_c.name = 'Demand for reduced flavin mononucleotide'
    dm_fmnh2_c.add_metabolites({fmnh2_c:-1})
    dm_fmnh2_c.lower_bound = 0
    dm_fmnh2_c.upper_bound = 1000
    with model as m:
        m.add_reactions([dm_fmnh2_c])
        m.objective = 'DM_fmnh2_c'
        solution = m.optimize()
        if solution.status != 'optimal' or solution.objective_value < 1e-6:
            print("Warning: model cannot produce fmnh2_c after adding FRDFMN_c reaction.")

mets_to_include = [
{
"id":"glu5sa_m",
"name":"L-glutamic 5-semialdehyde",
"compartment":"m",
"charge":0.0,
"formula":"C5H9N1O3",
"annotation":{
"bigg.metabolite":"glu5sa",
"biocyc":"META:L-GLUTAMATE_GAMMA-SEMIALDEHYDE",
"chebi":"CHEBI:17232",
"kegg.compound":"C01165",
"pubchem.compound":"193305",
"sbo":"SBO:0000247",
"seed.compound":"cpd00858",
"ymdb":"YMDB00918"
}
},
{
"id":"ppal_c",
"name":"propanal",
"compartment":"c",
"charge":0.0,
"formula":"C3H6O",
"notes":{
"yeastgem.compound":"s_4232",
},
"annotation":{
"bigg.metabolite":"ppal",
"biocyc":"META:PROPANAL",
"chebi":"CHEBI:17153",
"kegg.compound":"C00479",
"sbo":"SBO:0000247",
"seed.compound":"cpd00858",
"ymdb":"YMDB00918"
}
}
]
for met_data in mets_to_include:
    if met_data['id'] not in model.metabolites:
        new_met = cobra.Metabolite(id=met_data['id'])
        model.add_metabolites([new_met])
        for attr in met_data.keys() - {'id'}:
            if hasattr(new_met, attr) and isinstance(getattr(new_met, attr), dict) and isinstance(met_data[attr], dict):
                getattr(new_met, attr).update(met_data[attr])
            else:
                setattr(new_met, attr, met_data[attr])

rxns_to_include = [
# Added to help RBA model utilize ACACTm40_m since JL09_g3902 isn't in the proteomics data, and omitting these results in unrealistically high kapp values for ACACTm40_m
{
"id":"PYRDC3_c",
"name":"pyruvate decarboxylase (aldedyde-forming)",
"metabolites":{
"2obut_c":-1,
"co2_c":1,
"h_c":-1,
"ppal_c":1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g20 or JL09_g917",
"notes":{
"notes":"Added to allow RBA model to use ACACTm40_m for kapp calculations, since its isozymes are not in the proteomics data and their absence would hinder proteome fitting. Plus, yeast-GEM v9.1.1 includes this rxn.",
"manually_curated":"yes",
"yeastgem.reaction":"r_4651",
},
"annotation":{
"ec-code":"4.1.1.1",
"sbo":"SBO:0000176"
},
},
{
"id":"PPALDD2x_c",
"name":"aldehyde dehydrogenase (propanal, NAD)",
"metabolites":{
"ppal_c":-1,
"h2o_c":-1,
"h_c":2.0,
"nad_c":-1,
"nadh_c":1,
"ppa_c":1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g420 or JL09_g3633",
"notes":{
"notes":"Made irreversible just to minimize the risk of disrupting the model too much, even though its yeast-GEM version is reversible. Added to allow RBA model to use ACACTm40_m for kapp calculations, since its isozymes are not in the proteomics data and their absence would hinder proteome fitting. Plus, yeast-GEM v9.1.1 includes this rxn.",
"manually_curated":"yes",
"yeastgem.reaction":"r_4678",
},
"annotation":{
"ec-code":"1.2.1.5",
"sbo":"SBO:0000176"
},
},
# added since PYDXS_c has the glutaminase subunit missing from proteomics data, implying only 1 half of the rxn takes place
{
"id":"PYDXS1_c",
"name":"pyridoxal 5'-phosphate synthase (glutamine hydrolysing)",
"metabolites":{
"g3p_c":-1,
"nh4_c":-1,
"h2o_c":4,
"h_c":1,
"pi_c":1,
"pydx5p_c":1,
"r5p_c":-1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g3833",
"subsystem":"Vitamin B6 metabolism (ec00750)",
"annotation":{
"kegg.reaction":"R10088",
"sbo":"SBO:0000176"
}
},
# added since kapp calculations were struggling w/ absence of these rxns due to their enzymes not being in the omics data
{
"id":"FACOAL240_en",
"name":"fatty-acid--CoA ligase (lignoceric acid)",
"metabolites":{
"amp_en":1,
"atp_en":-1,
"coa_en":-1,
"ppi_en":1,
"ttccoa_en":1,
"ttcosa_en":-1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g1663",
"subsystem":"Fatty acid degradation",
"notes":{
"notes":"Added to allow RBA model to use fatty-acid--CoA ligase (lignoceric acid) for kapp calculations, since its isozymes are not in the proteomics data and their absence would hinder proteome fitting. Plus, yeast-GEM v9.1.1 includes this rxn.",
"manually_curated":"yes",
"yeastgem76.reaction":"r_2211",
"yeastkegg_id":"r_2211"
},
"annotation":{
"ec-code":"6.2.1.3",
"kegg.reaction":"R00390",
"sbo":"SBO:0000176"
}
},
# added to match atom mapping model (AMM); can paste rxns in from JSON for convenience
{
"id":"ACt_c_m",
"name":"acetate transport",
"metabolites":{
"ac_c":-1,
"ac_m":1
},
"lower_bound":-1000.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g1692 or JL09_g3775 or JL09_g191 or JL09_g192",
"subsystem":"Transport",
"notes":{
"yeastgem76.reaction":"r_1106",
"yeastkegg_id":"r_1106"
},
"annotation":{
"sbo":"SBO:0000185"
}
},
{
"id":"PPA_m",
"name":"inorganic diphosphatase",
"metabolites":{
"h2o_m":-1.0,
"h_m":1.0,
"pi_m":2.0,
"ppi_m":-1.0
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g3719",
"annotation":{
"ec-code":"3.6.1.1",
"metanetx.reaction":"MNXR100808",
"pubmed":"1648084",
"reference":"1648084",
"sbo":"SBO:0000649"
}
},
{
"id":"P5CD2_m",
"name":"delta-1-Pyrroline-5-carboxylate dehydrogenase",
"metabolites":{
"glu5sa_m":-1,
"glu__L_m":1,
"h2o_m":-1,
"h_m":2.0,
"nad_m":-1,
"nadh_m":1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g4292",
"subsystem":"Arginine and proline metabolism"
},
{
"id":"G5SADr_m",
"name":"L-glutamate 5-semialdehyde dehydratase",
"metabolites":{
"1pyr5c_m":1,
"glu5sa_m":-1,
"h2o_m":1,
"h_m":1
},
"lower_bound":-1000.0,
"upper_bound":1000.0,
"gene_reaction_rule":"SPONT",
"annotation":{
"kegg.reaction":"R03314",
"sbo":"SBO:0000176"
}
},
{
"id":"ALDD2y_c",
"name":"aldehyde dehydrogenase (acetylaldehyde, NADP)",
"metabolites":{
"ac_c":1,
"acald_c":-1,
"h2o_c":-1,
"h_c":2.0,
"nadp_c":-1,
"nadph_c":1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g420",
"annotation":{
"ec-code":"1.2.1.4",
"kegg.reaction":"R00711",
"sbo":"SBO:0000176"
}
},
{
"id":"MTHFD_c",
"name":"methylenetetrahydrofolate dehydrogenase (NADP)",
"metabolites":{
"methf_c":1,
"mlthf_c":-1,
"nadp_c":-1,
"nadph_c":1
},
"lower_bound":-1000.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g946",
"notes":{
"Confidence Level":"2",
"yeastgem.reaction":"r_0732"
},
"annotation":{
"ec-code":"1.5.1.5; 3.5.4.9; 6.3.4.3",
"kegg.reaction":"R01220",
"metanetx.reaction":"MNXR101749",
"sbo":"SBO:0000649"
}
},
{
"id":"MTHFC_c",
"name":"methenyltetrahydrofolate cyclohydrolase",
"metabolites":{
"10fthf_c":1,
"h2o_c":-1,
"h_c":1,
"methf_c":-1
},
"lower_bound":-1000.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g946",
"notes":{
"Confidence Level":"2",
"yeastgem.reaction":"r_0725"
},
"annotation":{
"ec-code":"1.5.1.5; 3.5.4.9; 6.3.4.3",
"kegg.reaction":"R01655",
"metanetx.reaction":"MNXR101748",
"sbo":"SBO:0000649"
}
},
{
"id":"GHMT2r_m",
"name":"glycine hydroxymethyltransferase",
"metabolites":{
"gly_m":1,
"h2o_m":1,
"mlthf_m":1,
"ser__L_m":-1,
"thf_m":-1
},
"lower_bound":-1000.0,
"upper_bound":1000.0,
"gene_reaction_rule":"JL09_g5520",
"annotation":{
"ec-code":"2.1.2.1",
"kegg.reaction":"R00945",
"sbo":"SBO:0000176"
}
},
{
"id":"OAADC_c",
"name":"oxaloacetate decarboxylase",
"metabolites":{
"co2_c":1,
"h_c":-1,
"oaa_c":-1,
"pyr_c":1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"io248592",
"subsystem":"Pyruvate metabolism",
"annotation":{
"ec-code":"4.1.1.112",
"kegg.reaction":"R00217"
}
},
{
"id":"G3PT_c",
"name":"glycerol-3-phosphatase",
"metabolites":{
"glyc3p_c":-1,
"glyc_c":1,
"h2o_c":-1,
"pi_c":1
},
"lower_bound":0.0,
"upper_bound":1000.0,
"gene_reaction_rule":"UNKNOWN",
"notes":{
"notes":"Added to gapfill and fix non-growth on either glucose or another substrate",
"Confidence Level":"2",
"yeastgem.reaction":"r_0489"
},
"annotation":{
"ec-code":"3.1.3.21",
"kegg.reaction":"R00841",
"metanetx.reaction":"MNXR99894",
"sbo":"SBO:0000649"
}
},
]
for rxn_data in rxns_to_include:
    if rxn_data['id'] not in model.reactions:
        # set attributes based on rxn_data
        new_rxn = cobra.Reaction(rxn_data['id'])
        if new_rxn.id not in model.reactions:
            model.add_reactions([new_rxn]) # putting it here allows old attributes to also be preserved when not overwritten
        new_rxn.add_metabolites({model.metabolites.get_by_id(met_id): coeff for met_id, coeff in rxn_data['metabolites'].items()})
        for attr in rxn_data.keys() - {'id','metabolites'}:
            # if attr already exists and is a dict, then update it rather than overwrite it
            if hasattr(new_rxn, attr) and isinstance(getattr(new_rxn, attr), dict) and isinstance(rxn_data[attr], dict):
                getattr(new_rxn, attr).update(rxn_data[attr])
            else:
                setattr(new_rxn, attr, rxn_data[attr])
        # check if mass and charge balanced
        imbal = check_mass_balance_cobra(new_rxn,model)
        if any(imbal.values()):
            print(f"Warning: {new_rxn.id} reaction is not mass/charge balanced: {imbal}")

model.reactions.GLYCL_m.bounds = (-1000, 1000) # make reversible based on AMM
model.reactions.FORt_c_m.bounds = (-1000, 1000) # make reversible based on AMM
model.reactions.ATPS_m.reaction = "3.0 adp_m + 10.0 h_c + 3.0 pi_m --> 3.0 atp_m + 3.0 h2o_m + 7.0 h_m" # update to match AMM
model.reactions.FECOOR_m.reaction = "2.0 focytC_m + 4.0 h_m + 0.5 o2_m --> 2.0 ficytC_m + h2o_m + 2.0 h_c" # update to match AMM

model.reactions.ACALDt_c_m.notes['notes'].replace("shutdown acetaldehyde transport to mitochondria to prevent cycling. no ethanol or acetaldehyde secretion observed in vivo (PMID:26479039,27678117,25888986)","")
msg = "[Made reversible based on yeast 9.0.2, and to resolve conflicts with IO's atom mapping model; while this may risk including cycling, it's more accurate to allow transport in both directions]"
if msg not in model.reactions.ACALDt_c_m.notes['notes']:
    model.reactions.ACALDt_c_m.notes['notes'] += msg
model.reactions.ACALDt_c_m.bounds = (-1000, 1000)

for gene in model.genes:
    # if 'uniprot' attribute doesn't exist in annotation, skip
    if "uniprot" not in gene.annotation.keys():
        continue
    # update GPRs for malate dehydrogenase; only A0A099NXM1 should catalyze MDH_x, A0A099P3W9 should catalyze MDH_m and A0A099P797 should catalyze MDH_c. See DeepLoc-2.1_MDH.pdf for details.
    if gene.annotation['uniprot'] == 'A0A099NXM1':
        for rxn in gene.reactions:
            if rxn.id == 'MDH_x':
                rxn.gene_reaction_rule = gene.id
        gene.notes.setdefault('notes', "")
        msg = "Assigned as the sole catalyst of MDH_x based on DeepLoc-2.1 subcellular localization prediction."
        if msg not in gene.notes['notes']:
            gene.notes['notes'] += "Assigned as the sole catalyst of MDH_x based on DeepLoc-2.1 subcellular localization prediction."
    elif gene.annotation['uniprot'] == 'A0A099P3W9':
        for rxn in gene.reactions:
            if rxn.id == 'MDH_m':
                rxn.gene_reaction_rule = gene.id
        gene.notes.setdefault('notes', "")
        msg = "Assigned as the sole catalyst of MDH_m based on DeepLoc-2.1 subcellular localization prediction."
        if msg not in gene.notes['notes']:
            gene.notes['notes'] += "Assigned as the sole catalyst of MDH_m based on DeepLoc-2.1 subcellular localization prediction."
    elif gene.annotation['uniprot'] == 'A0A099P797':
        for rxn in gene.reactions:
            if rxn.id == 'MDH_c':
                rxn.gene_reaction_rule = gene.id
        gene.notes.setdefault('notes', "")
        msg = "Assigned as the sole catalyst of MDH_c based on DeepLoc-2.1 subcellular localization prediction."
        if msg not in gene.notes['notes']:
            gene.notes['notes'] += "Assigned as the sole catalyst of MDH_c based on DeepLoc-2.1 subcellular localization prediction."
    # update GPR for glycine hydroxymethyltransferase based on AMM
    if gene.annotation['uniprot'] == 'A0A099P6F5':
        for rxn in gene.reactions:
            if rxn.id == 'GHMT2r_c':
                rxn.gene_reaction_rule = gene.id
        
for met in model.metabolites:
    if met.name == 'TDP': # some TDP versions have inconsistent CHEBI IDs
        met.annotation['chebi'] = 'CHEBI:58937'
    elif met.name == 'Potassium ion':
        met.annotation['chebi'] = 'CHEBI:29103'
    elif met.name == 'Manganese(II) ion':
        met.annotation['kegg.compound'] = 'C19610'
        met.annotation['pubchem.compound'] = '27854'
    elif met.name == 'Copper(II) ion':
        met.annotation['kegg.compound'] = 'C22424'
save_cobra_model(model,SM_path)