# changes made between 01_extract_enzymes_info_uniprot.ipynb and 02_build_model.ipynb, for reproducibility
# set paths relative to build_model/, since this script is run from 01_extract_enzymes_info_uniprot.ipynb
import json,os,sys, pandas as pd
# model_root_path = '../../'
# sys.path.append(model_root_path)
df_pro_path = 'input/PROTEIN_stoich_curation.tsv'
df_enz_rxns.index = df_enz_rxns['enz'].to_list()
df_pro = pd.read_csv(df_pro_path, index_col=0, sep='\t')
# let certain proteins appear in compartments where their catalyzed reactions occur, but that were not assigned in the previous step
df_pro.at['A0A099NYY2','subloc_assigned'] = ','.join(list(dict.fromkeys(['rm','l'] + df_pro.at['A0A099NYY2','subloc_assigned'].split(','))))
df_pro.at['A0A099P6G0','subloc_assigned'] = ','.join(list(dict.fromkeys(['vm'] + df_pro.at['A0A099P6G0','subloc_assigned'].split(','))))
df_pro.at['A0A099P5L4','subloc_assigned'] = ','.join(list(dict.fromkeys(['vm'] + df_pro.at['A0A099P5L4','subloc_assigned'].split(','))))
## While DeepLoc 2.1 doesn't seem to predict 'mm' as location for A0A099P7B6, it's needed since it catalyzes a reaction in the mitochondrial matrix where no alternatives are listed in the GSM model, and is a membrane-bound protein
df_pro.at['A0A099P7B6','subloc_assigned'] = ','.join(list(dict.fromkeys(['mm'] + df_pro.at['A0A099P7B6','subloc_assigned'].split(','))))
df_enz_rxns.at['JL09_g944','compartment-protein_maps'] = {'m': {'JL09_g944': 'JL09_g944_mm'}}
## While DeepLoc 2.1 doesn't seem to decisively predict 'vm' as a possible location, it's needed since it catalyzes a reaction in vm where no alternatives are listed in the GSM model, and is a membrane-bound protein
df_pro.at['A0A099P5C2','subloc_assigned'] = ','.join(list(dict.fromkeys(['vm'] + df_pro.at['A0A099P5C2','subloc_assigned'].split(','))))
## While DeepLoc 2.1 doesn't seem to decisively predict 'l' as a possible location (doesn't even include 'l' in the locations it assesses), these proteins are membrane-bound and catalyze a rxn there, so add 'l' to reflect limited space on its membrane
for pro in ['A0A099NWY9','A0A099NX66','A0A099NX19']:
	df_pro.at[pro,'subloc_assigned'] = ','.join(list(dict.fromkeys(['l'] + df_pro.at[pro,'subloc_assigned'].split(','))))
df_pro.at['A0A099NWY9','subloc_assigned'] = ','.join(list(dict.fromkeys(['rm'] + df_pro.at['A0A099NWY9','subloc_assigned'].split(','))))
## Added 'en' since it's membrane-bound and catalyzes a reaction in the ER membrane
df_pro.at['A0A099P715','subloc_assigned'] = ','.join(list(dict.fromkeys(['en'] + df_pro.at['A0A099P715','subloc_assigned'].split(','))))
## Adding compartments predicted by DeepLoc 2.1
df_pro.at['A0A099P136','subloc_assigned'] = ','.join(list(dict.fromkeys(['rm'] + df_pro.at['A0A099P136','subloc_assigned'].split(','))))
## Since it's listed as a "Multi-pass membrane protein", put JL09_g2999 in gm despite its rxn being in g
df_pro.at['A0A099NZ49','subloc_assigned'] = ','.join(list(dict.fromkeys(['gm'] + df_pro.at['A0A099NZ49','subloc_assigned'].split(','))))
## remove associated warnings for these changes
for gene in ['JL09_g3691','JL09_g1091','JL09_g609','JL09_g2999','JL09_g3545']:
	for warning in warnings_summary.copy():
		if f'Protein {gene} has location conflicts' in warning:
			warnings_summary.remove(warning)

# update mtDNA-encoded proteins based on BLASTn results against SD108 mtDNA
df_pro.at['A0A2U9RAL6','subloc_assigned'] = ','.join(list(dict.fromkeys(['c'] + df_pro.at['A0A2U9RAL6','subloc_assigned'].split(','))))
proteins = ['A0A2U9RAK4','A0A2U9RAL6','A0A2U9RAN3','A0A2U9RAR9','A0A2U9RAU9','A0A2U9RAV5','A0A2U9RB22','A0A2U9RB35','A0A2U9RBJ3','A0A2U9RBK5','A0A2U9RBW3','A0A2U9RCC6','A0A2U9RCF9','A0A2U9RCH4']
for pro in proteins:
	if pro == 'A0A2U9RAL6':
		df_pro.at[pro,'comments'] = '; '.join(list(filter(lambda x: str(x) not in ('','nan','NaN'), [str(df_pro.at[pro,'comments']), "c included as location since it's a cytosolic ribo subunit; seq inferred via BLASTn of mtDNA-encoded genes from strain CBS573 (uniprot ID from homolog) against strain SD108 mtDNA, followed by conversion into protein sequences"])))
	else:
		df_pro.at[pro,'subloc_assigned'] = ','.join(list(dict.fromkeys(['mm'] + df_pro.at[pro,'subloc_assigned'].split(','))))
		df_pro.at[pro,'comments'] = '; '.join(list(filter(lambda x: str(x) not in ('','nan','NaN'), [str(df_pro.at[pro,'comments']), 'mm included as location since NADHcplxI_c_m and ATPS_m involve transport across mitochondrial membrane; seq inferred via BLASTn of mtDNA-encoded genes from strain CBS573 (uniprot ID from homolog) against strain SD108 mtDNA, followed by conversion into protein sequences'])))
df_pro.at['A0A2U9RAK4','sequence'] = 'MKEENRMSMMFMNSPLEQFKMNTFMSLSSPMMDLSAMNVTTFSLYMVLVLLSMMFMLGLGLDRKLMMGSNWLMAMEATVDTMVNMVKGQMGGTMYGRYVPLIYTLFMFMLMGNLMGMMPYNFAMATSMVFIISMSLSLWMGLTMLGLYLHNVVFFSLFVPVGTPLALVPLLVLIELLSYSARAISLGLRLAANTLSGHLLMSMLGNLVKTFMGISKMTFIIGLLPIAGIFAIVILEFAISCMQAYVFAILTSSYLKDSLFLH'
df_pro.at['A0A2U9RAL6','sequence'] = 'MMKNMMKYMMMYSYKEKSEKRRWMMKREMMKYLELLNMRMNKWVKRMMMSRKRLLSRMNRLEKKNHMNDYKLMSYNFNKNSMMMKSMLSKLMMNLLENVMNGVDLKSGKSKMSGGNMMMSKPVMKENLNTVNMLFYYFLPNNKSYKYFNRMNMYLNKHMKNYKKLVKLNKNYSKSSYLLQDKLMMSNMFNNMMFNRNNMKYMNNNVNNLMSSLNMNSNNKSLYLSLLNKTLNNNLSNLYSNNNINKSIINNTLVNMSMLNNNYNNNNTTFNINNNLFNLLNLLNNNYNNNNNYNNMMSNTNMLLNNYNNNNNNNNNNKCNMHSKMKLERKNMMLNYLLKNEMMNKNQMWMSKMENNKMNNINNNNNMLKMNKENDKTNMFGYMMSYMDMLLGNMMNKSKKSNDMKMMMSKYFGLKEVNMTGMNLKYEFNNTEMLLKLMRKGMSKRKRTLSRMFRFRLKNRMPLLNDKGMLKNKMSNNLLKNLALNNNILNVEYNIKNVLFNSINNNNINNNNNNNMYNDMKDYYNNMSLENKKDLSLNNELLYKNMVGWSLLLKGKVGARKGKNRSNRMLMTKGSFKNNNLYIYNIFDDNSNNNGYTKDRLRLNYMKNSNFISYMDKSTNNGKLGMTLKVNIL'
df_pro.at['A0A2U9RAN3','sequence'] = 'MVNMMVFSTLMGSLMSGLMGRYMGMLYSKMMVCTSVMMSLLATYMLYYDVMMLDQEYRNMMMSWMNMEYLNMDWGFLMDKVSMSVLMPVVTMSSLVHMYAVVYMSHDPHQQRFFSYLSFFTFGMMVLVTGDNLLMLFLGWELMGCASYLLMSFWFTRMSAGKSGLNSLLLNRFGDTFFVMGLSLTMYLVGSLNFDTLFSLNSYLSTDMLTMMLMCMLMGCASKSVQFGLHTWLLNSMEGPTPVSSLLHAACLVMAGVYLLVRCSYMVEYSPMALMLMLVLGGMTTLMSGLMAVVSNDIKKIMALSTMSQVGIMMLGVGMSAYNLAMYHLICHSFFKALLFMSAGAIIHAVINEYQDMRTYGGFHKFLPLSYICIFIASLSLMALPGLTGYYSKDIIIESLYGSYTFTGYIIYWFAIASATLTSVYSIRLMYYIFFNIPNNPKYVYMHLNENPIIMYIPMIILATFSIFIGYLAKDLYLGLGATIYNNTFIHPNHLIIIDTEFSLSSFIKLLPLITSVIFSTILLIFYELFYHHTFIFNNKYTRQFYNFFNQKLFFDQILNNYGFRSWLYISYLLNTYIDKGVLIFLGPYGLSALSLKISNAINQLVFFNLGLILELILFSIFVFIFFNSYTSLYIIIAIIIFIFI'
df_pro.at['A0A2U9RAR9','sequence'] = 'MNMFMNKMWNDVPTPWGMYFQDSATPNMEGMMELHNNMMFYLCMMLGLVSYMFFNVMTVYNNSKMAYKYMNHGQFMEMVWTTFPAMMLLMMAFPSFMLLYMCDEVMAPAMTMKAMGLQWYWKYEYSDFINEKGETMEFESYMVPEDLLEEGQLRMLDVDASVVCPVDTHIRFIVTSADVIHDFCMPSLGVKVDASPGRLNQTSALIQREGVYYGQCSELCGVMHSAMPIKIEAVSLGEFLSWIDEQ'
df_pro.at['A0A2U9RAU9','sequence'] = 'MELILLVLGSMGLNMMDIYMMLEMLMMGCTVNSMWMSNMNNDMMGLLYSLMQMMMAGMESAMGLSMLVSFNKMRGTDEMLRWL'
df_pro.at['A0A2U9RAV5','sequence'] = 'MPQLVPFFFLNTLFYGYMGLFVMLLLVSSMMLPYMLKMSMMRNMMVKL'
df_pro.at['A0A2U9RB22','sequence'] = 'MLSMSMLMFVVLTSMMEKSNECKKWLKRVCLMSLMNVLLMNLSMNSTEMLMLGHSMYNSLLYMDSMNNYMMELMLLVGMMYMLLLSCNMTTMKVLDHMMQVKESLGLLMLFNMLGLVLLMLVNDLMMLFMVMELQSYSLYLMTSSYNSSSNSVKSGLYYFMVGSMGSFMMLDGTVSMYEEMGLTNMELMSLCMNNMNMYDMLMMLFGLFMKVGLAPFYSYSMVIYTLAPTMVTYYMSLTPKLSMLTMMLSLMSMMNLMNTDSSSNIINIISFIILLSMMIGSMGGMNVMKIKTLLAYSSLLNVSYMLLAIISNNGYSYIGYLFYILQYSITHINTFSIILLLPIYSHINLLNNNNNNNNNNNSNNNNNNNYDVMLSKYSPMEYISQFKVLLNKNSYLCMALCMSFFSLMGMPPTSGFYGKLLMLMSALGNGYMFLSMLLVISSSMSAYYYGSVMKELCFDMLNNYNVFNYKDMSMYNKSLLSKMESMYMTLKDNMMSKNGNNMLLMNNNMSYQMSLMTLMIMLIMFEYDNMMRGTFMVTLNMM'
df_pro.at['A0A2U9RB35','sequence'] = 'MFMDFNILDVTLMGYYMLMSTTSFESLGMLVGTMMPVMTELWFNSSQLFGMLYTLMYVGAVVVLFVFMLSVLNKKEEYEPQYSSLLMMLLMLFLDDLDYNDSNDCDMNNMNIYDNFNLLNIHNNSDLTTMGSLMFTEYSLVLILIALILLLGMMGMMVMMKSSSSLKKD'
df_pro.at['A0A2U9RBJ3','sequence'] = 'MNMRSRSNLQLFPFHLVAPSPWPMLMSFALLSFMLTFALTMHGYMLSLWWLLLSASVVLATMFLWTRDMVVEGTYLGDHTMAVRRGLNLGFLLFVGSEVAMFMSLFWAYFHSAMGPTMEMGATWPPVGMTSMKPTELPLLNTMMLLASGATVTWAHHAMLSKDRNGTLLGLLMTTMLIVVFVLCQMCEYTFATFTIADSVYGSVFYAATGLHFMHMVMLAMMLLVCYARMYLYHFTSTHHLNLETTILYLHVLDIIWLFLYIVFYWWGC'
df_pro.at['A0A2U9RBK5','sequence'] = 'MLAYRKRQPYMSLVNSYLMDSPQPSSLNYWYNLGSLLGLCLVMQMASGVFLAMHYSSHMELAFDSVEHMMRDVNLGWLMRYMHANGASFFFGCMYMHVGKALYYGSYRKPRVLVWMMGVMMLMATMATAFMGYCTVYGQMSHWGATVMTNLLSAMPFIGADIVPFMWGGFSVSNPTMQRFFALHYLLPFMLAALVMLHLMALHVHGSSNPMGMTGNMDRLSMHPYFMFKDAMTMFVFLMLFCFFVFFSPNTLGHPDNYIPGNPLVTPASIVPEWYLLPFYAILRSIPDKLGGVLAMFGALLILLILPITDRSLYRGLAFKPISKFFFWLFVFNFILLGNLGQLHVEVPFMQLGVNCTFFYFAYFLIIVPGISTLENILFYLGRIN'
df_pro.at['A0A2U9RBW3','sequence'] = 'MQLVLAGKYMGAGMATMGLTGAGMGMAMVFAALMNGTSRNPGLRNTLFPFAILGFATSEATGLFCLMISFLLLYGV'
df_pro.at['A0A2U9RCC6','sequence'] = 'MIYTLMEWMESMEYLLMAMMFMVGLLLSVAYLTVAERKTMAYMQRRLGPNAVGYYGVLMAMADALKLLVKEMVLVHNGEYMYMMGGPLMTLYSVLMSFGVVPLNTGLGMMEGEYSLMFLLGSGSLGVLGTVMVGWMSNSKYTMLASVRTTAQLMSYELMLTTMVFMMVLLVNSLNMELMMENQVYMNYMMPLLPLCMMYFMAALAETARPPFDNMEAESELVSGHMTELSASPFVIFFLSEYTSMVLMSVLTMTFFFGGYLSWPLDINNLLFLDLLPQNNLYWFIEGLIFSGLFIMKTNFFLFTYIWVRASFPRLRFDLLILFCWYVLLPLLFAFCLFMPCVLFAGNGLALIV'
df_pro.at['A0A2U9RCF9','sequence'] = 'MKHQLTFGDYINRWLFSTNAKDMAVLYFMFSLFCGMLGSFMSLVLRLELAAPGNQLLMGNHQLFNVMVTAHAMLMVFFLMMPVTMGFFGNYLVPLMMGASDMSFARLNNMSFWTLPPALMCLMSSALMENGAGTGWTVYPPLSGMQSHSGPSVDMAMFSLHLTSMSSLLGAMNFMTTMFNMRTLGMTMSKLPLFVWSVLFTAMLLVMTLPVLSAGVTLLLMDRNFNTSFYESAGGGDPMLYQHLFWFFGHPEVYILMLPGFGMVSHVVSTYSKKPVFGQMGMMYAMGSMGFLGLLVWSHHMYVVGLDVDSRAYFTSATMVMAVPTGMKMFSWLATLYGGSLRFTTPMLYGVAFMFMFTVGGLSGVLLSNASLDMAFHDTYYVIGHFHYVLSLGALFAIFAGYYYWSPMVTGLYYNNVLSKIQFWLFFIGANVTFLPMHFLGMNGMPRRIPDYPDAFAGFNAICSFGAVLSIISLLFFGYVIYDQLINGTANKMYSTNSLLKDPDFFESNNTFNSNEVKSESIEFLLNYPPMFHTFNTVAMQS'
df_pro.at['A0A2U9RCH4','sequence'] = 'MLNLMMLNTLGMMSMMGLTKKRMLWMKANKNMMSKSMLLLCMMSLYMVMNELNSMNMNDYGYLRVMNMLNMDYGYDGMGLMFMMLTMYLGVMSMMSNWENMMNKRGTYMMLLMLLVMLMGLNFMCLDLMSFYMFFESTLMPLYLMMGMYGGSNKNKAAYYVLIYTLCSSLFMLLSMMLMYVMYNNVDYVTMNSKMMSMELQSVLLMGLMMGLGVKTPLVPVHTWLPTVHSESPLGGSMLLAGMVLKLAMYGMMRMLLSNMSEMSLLYAPMLYVMGMMTMMYSSLVTLSQTDLKVIMAYSSISHLGVCMLGMFANNLLGIIGSLILCLAHGFVSPALFILMGGVLYDRYHNRVLFYYQGLTQFMPLFAVYLVFFSFSNTGTPLSYNFIGEFLCITGAIISNPTIGIFTALSVLLSACYQMKLTNKLTSGTFSKYLHVTYDLLYRENFILFFLALPTLVLGIVPNFINNILYLSSNALIYSIVI'

df_pro.at['P40011','sequence'] = 'MTKETKKIIQVLSQYATCDMSDALVKHQYKLGGSLPNLHRYSANETKTTMVGKAYTVLFAPKDDPREAISGGYIDKLPEDSVLVIGTTLNLQKKDAPFTKMTNALYGGLMSTRANYLRSRGTVVFGRVRDIDEHLELKRDVFSYGLGTAAHGPVVKMVGINVPLEVLVDGYPEPHVETINPGDIIMGDDNGVVRVADSEELSSIMEYIPKRVQADTLVKEDILNGLECTPAQKKRRAGL'
df_pro.at['P40011','MW (g/mmol)'] = protein_mw_from_sequence([df_pro.at['P40011','sequence']], path='input/PROTEIN_amino_acid_map.txt')[0]
# update MW of io248592 enzyme
df_enz_rxns.at['io248592','MW (g/mmol)'] = df_pro.at['P40011','MW (g/mmol)'] * 3
# remove sheme_c from PTM_stoich_edits of A0A2U9RBK5, since it doesn't currently have a modeled version in the mitochondria
df_pro.at['A0A2U9RBK5','PTM_stoich_edits'] = 'hemeO_m:2,pheme_m:2|hemeA_m:1,hemeO_m:1,pheme_m:2|hemeA_m:2,pheme_m:2'
# save updated df_pro
df_pro.to_csv(df_pro_path, sep='\t')

# no evidence that a homotrimer of JL09_g4214 is involved in this rxn, despite its subunit structure on UniProt indicating it's involved in its own rxn, just to be safe
df_enz_rxns.at['JL09_g3328JL09_g3622JL09_g4214','protein_stoich'] = df_enz_rxns.at['JL09_g3328JL09_g3622JL09_g4214','protein_stoich'].replace('JL09_g4214:3','JL09_g4214:1')
# listed as a heterodimer on UniProt; and despite FRA2's identity being unclear, better to assume 1 copy of JL09_g1478 than 2 copies to reduce overestimation of protein costs
df_enz_rxns.at['JL09_g1478JL09_g4909','protein_stoich'] = df_enz_rxns.at['JL09_g1478JL09_g4909','protein_stoich'].replace('JL09_g1478:2,','JL09_g1478:1,')
df_enz_rxns.at['JL09_g1478JL09_g2712','protein_stoich'] = df_enz_rxns.at['JL09_g1478JL09_g2712','protein_stoich'].replace('JL09_g1478:2,','JL09_g1478:1,').replace('JL09_g2712:2','JL09_g2712:1')
## divide MW by 2
df_enz_rxns.at['JL09_g1478JL09_g2712','MW (g/mmol)'] = df_enz_rxns.at['JL09_g1478JL09_g2712','MW (g/mmol)'] / 2
# for all enzymes containing JL09_g2504, apply "Tetramer of two components I and two components II"; if there are 2 different proteins in their protein_stoich, set both their values to 2 (otherwise, set the single protein's value to 1)
for enz in df_enz_rxns.index:
	if 'JL09_g2504' in df_enz_rxns.at[enz,'protein_stoich']:
		proteins_in_enz = df_enz_rxns.at[enz,'protein_stoich'].split(',')
		if len(proteins_in_enz) == 2:
			new_stoich_list = []
			for pro_stoich in proteins_in_enz:
				pro, stoich = pro_stoich.split(':')
				new_stoich_list.append(f'{pro}:2')
			df_enz_rxns.at[enz,'protein_stoich'] = ','.join(new_stoich_list)
		elif len(proteins_in_enz) == 1:
			df_enz_rxns.at[enz,'protein_stoich'] = f'JL09_g2504:1'
# remove all copies of warning for duplicate allele JL09_g6147, since it corresponds to the same protein as another allele
# remove warnings related to enzymes with extra subunits I couldn't find, or model-database mismatches that aren't important for our purposes
warnings_to_remove = ['JL09_g6147 not found in df_pro',"['JL09_g6147 not found in df_pro','Protein JL09_g6147 not found in df_pro for reaction ATPS_m with compartments {'m', 'c'}.']",
						'Enzyme JL09_g3024 for reaction METAT_c has 1 subunits, while Uniprot subunit structure predicts 4 subunits.',
						'Enzyme JL09_g1785 for reaction UAGPT2_c has 1 subunits, while Uniprot subunit structure predicts 2 subunits.',
						'Enzyme JL09_g5771 for reaction UAGPT2_c has 1 subunits, while Uniprot subunit structure predicts 2 subunits.',
						'Enzyme JL09_g3162 for reaction C3STDH2_c has 1 subunits, while Uniprot subunit structure predicts 4 subunits.',
						'Enzyme JL09_g3162 for reaction C3STDH1_c has 1 subunits, while Uniprot subunit structure predicts 4 subunits.',
						'Enzyme JL09_g2412 for reaction DNMPPA_c has 2 subunits, while Uniprot subunit structure predicts 4 subunits.',
						'Enzyme JL09_g1478JL09_g2712 for reaction GTHO_c has 3 subunits, while Uniprot subunit structure predicts 2 subunits.',
						'Enzyme JL09_g1478JL09_g2712 for reaction GTHO_c has 4 subunits, while Uniprot subunit structure predicts 2 subunits.',
						'Enzyme JL09_g4083 for reaction AMID_c has 1 subunits, while Uniprot subunit structure predicts 3 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3071JL09_g3593JL09_g4155JL09_g6155JL09_g643JL09_g920iomATP8iomATP9 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3593JL09_g4155JL09_g5735JL09_g6155JL09_g643JL09_g920iomATP8iomATP9 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3071JL09_g3593JL09_g6155JL09_g643JL09_g920iomATP6iomATP8iomATP9 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3593JL09_g5735JL09_g6155JL09_g643JL09_g920iomATP6iomATP8iomATP9 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3071JL09_g3593JL09_g6155JL09_g643JL09_g920iomATP8iomATP9iomCOX1 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3593JL09_g5735JL09_g6155JL09_g643JL09_g920iomATP8iomATP9iomCOX1 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme ATPS_mENZ_JL09_g4155JL09_g5735 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme ATPS_mENZ_JL09_g5735iomATP6 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme ATPS_mENZ_JL09_g5735iomCOX1 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme ATPS_mENZ_JL09_g3071JL09_g4155 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme ATPS_mENZ_JL09_g3071iomATP6 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme ATPS_mENZ_JL09_g3071iomCOX1 for reaction ATPS_m has 13 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3071JL09_g3593JL09_g4155JL09_g6155JL09_g643JL09_g920 for reaction ATPS_m has 11 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Enzyme JL09_g1386JL09_g2023JL09_g219JL09_g2239JL09_g2373JL09_g3593JL09_g4155JL09_g5735JL09_g6155JL09_g643JL09_g920 for reaction ATPS_m has 11 subunits, while Uniprot subunit structure predicts 17 subunits.',
						'Charge mismatch for leutrna_c: model has charge 0, database has charges [1]',
						'KEGG Compound mismatch for dextrin_e: model has C01935, ChEBI 28675 has [13000]',
						'KEGG Compound mismatch for dextrin_en: model has C01935, ChEBI 28675 has [13000]',
						## ignore since it's not crucial so long as other options exist
						'Could not find metabolite ID for CHEBI ID 24480 when updating protein stoichiometry',
						'CHEBI ID 23380 is potential match for model metabolite(s) and cofactor(s); consider adding it to the GSM model and rerunning this script.',
						'CHEBI ID 20882 is potential match for model metabolite(s) and cofactor(s); consider adding it to the GSM model and rerunning this script.',
						'CHEBI ID 49550 is potential match for model metabolite(s) and cofactor(s); consider adding it to the GSM model and rerunning this script.',
						'KEGG Compound mismatch for ficytC_m: model has C00125, ChEBI 15719 has [14293]',
						'KEGG Compound mismatch for ficytC_m: model has C00125, ChEBI 16951 has [14322]',
						'KEGG Compound mismatch for focytC_m: model has C00126, ChEBI 15856 has [14363]',
						'Charge mismatch for ca2_e: model has charge 2, database has charges [0]',
						'Charge mismatch for ca2_c: model has charge 2, database has charges [0]',
						'Charge mismatch for cu2_e: model has charge 2, database has charges [0]',
						'Charge mismatch for cu2_c: model has charge 2, database has charges [0]',
						'Charge mismatch for zn2_e: model has charge 2, database has charges [0]',
						'Charge mismatch for zn2_c: model has charge 2, database has charges [0]',
						'Charge mismatch for ca2_rm: model has charge 2, database has charges [0]',
						'Charge mismatch for cu2_mm: model has charge 2.0, database has charges [0]',
						'Charge mismatch for cu2_m: model has charge 2.0, database has charges [0]',
						'Charge mismatch for cobalt2_e: model has charge 2, database has charges [0]',
						'Charge mismatch for cobalt2_c: model has charge 2, database has charges [0]',
						'Putting CHEBI ID 60052 in compartment c when updating protein stoichiometry, since compartment m not found in met_chebi_map',
						'Could not find metabolite ID for CHEBI ID 60052 when updating protein stoichiometry',
						'Metabolite with CHEBI ID 24480 not in c, so compartment m used when updating protein stoichiometry'
					  ]
while any(warning in warnings_summary for warning in warnings_to_remove):
    for warning in warnings_to_remove:
        if warning in warnings_summary:
            warnings_summary.remove(warning)