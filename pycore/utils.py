# RBA utility functions
import re
from dataclasses import dataclass, field
from typing import Dict

max_id_length = 63 # GAMS has a limit of 63 characters for variable names
ptmpro_prefix = 'PTMPRO'
max_id_length_for_ptmpro = max_id_length - len(ptmpro_prefix) - 6 # length available for the base ID and PTM info after adding the PTMPRO prefix and "WASTE" and hyphen

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def metabolites_dict_from_reaction_equation_RBA(eqn, metabolites_to_ignore=[], split=False):
    import re
    rStr, pStr = re.split('-->|->|<--|<-|<=>|<->', eqn)
    rStr = rStr.strip(' ')
    pStr = pStr.strip(' ')
    rs = re.split(r' \+ | \+|\+ |\+', rStr)
    ps = re.split(r' \+ | \+|\+ |\+', pStr)

    r_dict = dict()
    for r in rs:
        if ' ' in r:
            val, met = re.split(r'\s+', r)
            if is_number(val):
                r_dict[met] = -float(val)
            else:
                r_dict[met] = '-' + val
        else:
            r_dict[r] = -1.0

    p_dict = dict()
    for p in ps:
        if ' ' in p:
            val, met = re.split(r'\s+', p)
            if is_number(val):
                p_dict[met] = float(val)
            else:
                p_dict[met] = val
        else:
            p_dict[p] = 1.0

    if metabolites_to_ignore:
        for met in metabolites_to_ignore:
            r_dict.pop(met, None)
            p_dict.pop(met, None)
    if split:
        return r_dict, p_dict
    else:
        return merge_two_dicts(r_dict, p_dict)

def build_reaction_equation_from_metabolites_dict_RBA(met_dict, arrow='<=>', floatdecimal=None):
    """Takes dict (keys=metabolite IDs, values=stoichiometric coefficients) and returns reaction equation string"""
    lhs = []; rhs = [];
    floatdecimal_expr = r'{:.' + str(floatdecimal) + 'f}' if floatdecimal is not None else '{}'
    for k,v in met_dict.items():
        if is_number(v):
            v = float(v)
            if v == -1:
                lhs.append(k)
            elif v == 1:
                rhs.append(k)
            elif v < 0 and v != -1 and v.is_integer():
                lhs.append(' '.join([str(-int(v)), k]))
            elif v > 0 and v != 1 and v.is_integer():
                rhs.append(' '.join([str(int(v)), k]))
            elif v < 0 and v != -1:
                lhs.append(' '.join([(floatdecimal_expr).format(-v), k]))
            elif v > 0 and v != 1:
                rhs.append(' '.join([(floatdecimal_expr).format(v), k]))
        else:
            if v[0] == '-':
                lhs.append(v[1:] + ' ' + k)
            else:
                rhs.append(v + ' ' + k)
                
    return ' '.join([ ' + '.join(lhs), arrow, ' + '.join(rhs)])

@dataclass
class met_ID:
    """Dataclass for metabolite naming scheme with methods to assemble and parse names.
    
    :param tag: The tag for the metabolite (default "MET").
    :param base_id: The base identifier for the metabolite (i.e., the name without any protein isoform, cleavage, PTM, or compartment information).
    :param isoform: The isoform number (if applicable). Mainly useful for splicing isoforms of proteins.
    :param cleavage_ranges: A tuple of tuples indicating cleavage ranges (if applicable). Each inner tuple contains start and end positions of amino acids included.
    :param PTMs: A dictionary of cofactors/modifications (or string in the format of "PTM_stoich_edits" column) with their quantities. Use IDs from the stoichiometric model.
    :param compartment: The compartment identifier (default "c" for cytosol).

    Example usage:
        met = met_ID(tag="PTMPRO", base_id="ProteinX", isoform=2, cleavage_ranges=((10, 50), (70, 100)), PTMs={"CofactorA": 2, "CofactorB": 1}, compartment="m")
        print(str(met))  # Outputs: PTMPRO-ProteinX_I2_AAs10e50s70e100__PTM_CofactorA_x2_CofactorB_x1_m
        
        parsed_met = met_ID.from_name("PTMPRO-ProteinX_I2_AAs10e50s70e100__PTM_CofactorA_x2_CofactorB_x1_m")
        print(parsed_met)  # Outputs: met_ID(tag='PTMPRO', base_id='ProteinX', isoform=2, cleavage_ranges=((10, 50), (70, 100)), PTMs={'CofactorA': 2, 'CofactorB': 1}, compartment='m')
    """
    tag: str = "MET"
    base_id: str = ""
    isoform: int = None
    cleavage_ranges: tuple[tuple[int, int]] = None
    PTMs: Dict[str, int] = field(default_factory=dict)
    compartment: str = "c"

    # important info for proteins: splicing isoforms, cleavage, PTM stoichiometry
    isoform_regex = re.compile(r"_I(\d+)(?=_|$)")
    # cleavage_regex should capture multiple ranges, e.g., _AAs10e50s20e100_ for two ranges, AAs10e50 for 1, AAs1e9s20e30s40e50 for 3, etc.
    cleavage_regex = re.compile(r"_AAs((\d+e\d+)(s\d+e\d+)*)(?=_|$)")
    cleavage_range_regex = re.compile(r"(\d+)e(\d+)")
    # PTM section starts with "_PTM_" followed by each metabolite and its quantity separated by "x", e.g., _PTM_CofactorA_x2_CofactorB_x1
    PTM_prefix = "_PTM_"
    PTM_separator = "_"
    PTM_quantity_separator = "_x"

    def __str__(self) -> str:
        """Assembles the naming scheme."""
        # add isoform string if applicable
        full_ID_elements = {'start': f"{self.tag}-{self.base_id}", 'isoform': f"_I{self.isoform}_" if self.isoform is not None else "", 'cleavage': "", 'PTMs': "", 'compartment': f"_{self.compartment}" if self.compartment else ""}
        
        def full_ID_from_elements(elements):
            return ''.join([elements['start'], elements['isoform'], elements['cleavage'], elements['PTMs'], elements['compartment']])
        # add cleavage ranges if applicable
        if self.cleavage_ranges is not None:
            cleavage_str = "s".join([f"{start}e{end}" for start, end in self.cleavage_ranges])
            full_ID_elements['cleavage'] = f"_AAs{cleavage_str}_" if cleavage_str else ""
        # add cofactors/modifications if applicable
        if self.PTMs:
            if type(self.PTMs) is str:
                # convert from string representation to dict
                self.PTMs = {k:v for k,v in (item.split(':') for item in self.PTMs.split(','))}
            PTMs_str = self.PTM_separator.join([
                f"{ptm}{self.PTM_quantity_separator}{qty}" for ptm, qty in self.PTMs.items()
            ])
            full_ID_elements['PTMs'] = f"{self.PTM_prefix}{PTMs_str}"
            if len(full_ID_from_elements(full_ID_elements)) > max_id_length_for_ptmpro:
                # replace PTMs_str with version where qty of 1 uses "_" instead of "x1"
                PTMs_str = self.PTM_separator.join([
                    f"{ptm}__" if float(qty) == 1 else f"{ptm}{self.PTM_quantity_separator}{qty}" for ptm, qty in self.PTMs.items()
                ])
                full_ID_elements['PTMs'] = f"{self.PTM_prefix}{PTMs_str}"
        return full_ID_from_elements(full_ID_elements)

    @classmethod
    def from_name(cls, name: str):
        """Extracts parts from an existing name string."""
        # Split the main part into each section
        tag, rest = name.split('-', 1)
        # Split by the last underscore to isolate the compartment
        compartment_start = rest.rfind('_')
        base_id = None
        isoform = None
        cleavage_ranges = None
        PTMs = dict()
        # Check for isoform
        isoform_match = cls.isoform_regex.search(name)
        if isoform_match:
            isoform = int(isoform_match.group(1))
            base_id = rest.split(isoform_match.group(0))[0] if not base_id else base_id
            # set compartment_start to the position after the isoform match
            compartment_start = rest.find(isoform_match.group(0)) + len(isoform_match.group(0))
        # Check for cleavage ranges
        cleavage_match = cls.cleavage_regex.search(name)
        if cleavage_match:
            cleavage_ranges = []
            for range_match in cls.cleavage_range_regex.finditer(cleavage_match.group(1)):
                start, end = int(range_match.group(1)), int(range_match.group(2))
                cleavage_ranges.append((start, end))
            base_id = rest.split(cleavage_match.group(0))[0] if not base_id else base_id
            compartment_start = rest.find(cleavage_match.group(0)) + len(cleavage_match.group(0))
        # Check for PTMs (keep in mind that some names may omit PTM info if too long)
        ptm_start = name.find(cls.PTM_prefix)
        if ptm_start != -1:
            # PTM section found (after "_PTM_", before compartment)
            ptm_section_raw = re.search(r'(_PTM_.*)', name)
            ptm_section = re.sub(r'__(?=_|$)', '_x1', ptm_section_raw.group(1)) # replace "__" with "_x1" to indicate quantity of 1 for easier parsing
            ptm_section = re.search(r'_PTM_(.+_x-*\d+\.*\d*)', ptm_section)
            if ptm_section:
                # ptm_items = re.sub(r'__(?=_|$)', '_x1', ptm_section.group(1))
                ptm_items = ptm_section.group(1)
                ptm_quantities = re.split(r'_x(-?\d+\.*\d*)(_|$)', ptm_items)
                ptm_quantities = [item for item in ptm_quantities if item and item != '_']
                # each quantity is preceded by its metabolite ID
                for i in range(0, len(ptm_quantities)-1, 2):
                    ptm_id = ptm_quantities[i]
                    qty = float(ptm_quantities[i+1]) if '.' in ptm_quantities[i+1] else int(ptm_quantities[i+1])
                    PTMs[ptm_id] = qty
                # set base_id to rest without PTM section
                base_id = rest[:ptm_start - len(tag) - 1] if not base_id else base_id
                compartment_start = rest.find(ptm_section_raw.group(0)) + len(ptm_section_raw.group(0))

        if base_id is None:
            base_id = rest
        compartment = rest[compartment_start+1:] if compartment_start != -1 and compartment_start + 1 < len(rest) else None
        if compartment:
            global model
            import sys
            if 'model' not in globals():
                model = sys.modules['__main__'].__dict__.get('model', None)
            if model is None:
                from pycore import RBA_defaults
                # update global model variable if not already set
                model = RBA_defaults.load_SM()
            if compartment not in model.compartments:
                compartment = None
        base_id = base_id[:len(base_id)-len(compartment)-1] if compartment and base_id.endswith('_' + compartment) else base_id
        return cls(tag=tag, base_id=base_id, isoform=isoform, cleavage_ranges=cleavage_ranges, PTMs=PTMs, compartment=compartment)

def extract_details_from_rxnid(rxn_id):
    idsplit = rxn_id.split('-')
    tag = idsplit[0]
    rxn_base_id = idsplit[1]
    enz_id = rxn_id[(len(tag) + len(rxn_base_id) + 2):]
    rxn_dir = rxn_base_id.split('_')[-1]
    rxn_base_id = rxn_base_id[:-len(rxn_dir)-1]
    return(tag,rxn_base_id,rxn_dir,enz_id)

def build_uniprot_query(query: str = None, proteome_id: str = None, uniprot_ids: list[str] = None):
    """
    Builds a Uniprot query URL for the given proteome ID or list of Uniprot IDs.
    
    :param query: A custom query string (optional). If provided, it will be used as the query string, and the other query parameters will be ignored.
    :type query: str
    :param proteome_id: The ID of the proteome to query (e.g., UP000029867).
    :type proteome_id: str
    :param uniprot_ids: List of Uniprot accessions to query (e.g., ['A0A2U9RCF9', 'Q67890']).
    :type uniprot_ids: list[str]
    """
    # field names determined via "Advanced" option for searches on UniProt
    # cc_alternative_products for alternative splicing
    # go_c for cellular component GO terms
    fields = 'accession,reviewed,id,protein_name,gene_names,organism_name,annotation_score,ec,protein_families,sequence,mass,cc_sequence_caution,ft_chain,cc_ptm,cc_cofactor,cc_subunit,encoded_in,cc_subcellular_location,go_c,ft_intramem,ft_transmem,go_id,keyword,cc_function,cc_activity_regulation,ph_dependence,cc_miscellaneous,cc_alternative_products'
    if not query:
        # if both proteome_id and uniprot_ids are provided, use both
        if proteome_id and uniprot_ids:
            query = f'proteome:{proteome_id} OR {" OR ".join([f"accession:{id}" for id in uniprot_ids])}'
        elif proteome_id:
            query = f'proteome:{proteome_id}'
        elif uniprot_ids:
            query = f'{" OR ".join([f"accession:{id}" for id in uniprot_ids])}'
        else:
            raise ValueError('Either query, proteome_id, or uniprot_ids must be provided')
    return f'https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields={fields}&format=tsv&query=%28{query}%29'

def protein_mw_from_sequence(sequences: list[str], path: str = None) -> float:
    """Calculate the molecular weight of a protein from its amino acid sequence."""
    # find PROTEIN_amino_acid_map.txt in the same dir as SM_path
    import os, pandas as pd
    if path is None:
        path = os.path.join(os.path.dirname(SM_path), 'PROTEIN_amino_acid_map.txt')
    df_aamap = pd.read_csv(path, sep='\t')
    df_aamap.index = df_aamap.aa_abbv.to_list()
    aa_weights = df_aamap['MW'].to_dict()
    mws = [sum(aa_weights.get(aa.upper(), 0) / 1000 for aa in seq) for seq in sequences]
    return mws