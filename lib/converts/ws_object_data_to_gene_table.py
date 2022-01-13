import sys
import os
import json
import logging
import pandas as pd


"""
Gene table looks like

            locusId (str):sysName (str):type (int):scaffoldId (str):begin (int):end (int):
                strand (str +/-):name (str):desc (str):GC (float [0,1]):nTA (int)
"""

def obj_data_to_gene_table(obj_data_fp, cdss_method=True,
                            parse_fail_limit=0.1):
    """
    Args:
        obj_data_fp: Path to JSON data generated by KBase
        which contains all the info needed to generate 
        a genes table
    Returns:
        gene_table_df (Pandas DataFrame): A pandas dataframe of all the genes info.
                                        Has the following columns:
            locusId (str):sysName (str):type (int):scaffoldId (str):begin (int):end (int):
                strand (str +/-):name (str):desc (str):GC (float [0,1]):nTA (int)

            Note: Possible to add: AA_seq (str)
    Description:
        
    """

    # We load the JSON object
    full_object_dict = json.loads(open(obj_data_fp).read())
    if not isinstance(full_object_dict, dict) or not "data" in full_object_dict:
        raise Exception("Program expects Genome JSON object dict " + \
                        "to be dict with key 'data'.")
    
    data_list = full_object_dict['data']
    if len(data_list) > 1:
        logging.warning("Expecting data list to have length 1, " + \
                        "instead length: " + str(len(data_list)))
    data_object = data_list[0]
    data_level_3 = data_object['data']

    if cdss_method:
        cdss_list = data_level_3['cdss']
        gene_table_list = []
        nCDS = len(cdss_list)
        if nCDS == 0:
            raise Exception("No CDSs found in this genome object. Cannot " + \
                            "create gene table. Please contact developers " + \
                            "with further questions.")
        logging.info(f"Parsing {nCDS} CDSs")
        nFails = 0
        for CDS_info in cdss_list:
            try:
                gene_table_list.append(parse_CDS_info(CDS_info))
            except Exception as inst:
                nFails += 1
                logging.critical(CDS_info)
                logging.info("Failed to parse above CDS info, due to below:")
                logging.info(inst)
        
        # We know nCDS > 0
        parse_fail_ratio = nFails/nCDS
        if parse_fail_ratio > parse_fail_limit:
            raise Exception("Failed to succesfully parse KBase genome object." + \
                            " Please contact developers.")
        logging.info("Failed to parse the following ratio of CDSs: " + \
                     f"'{nFails/nCDS}'.")

        # Converting gene_table_list into pandas dataframe
        gene_table_d = {x:[] for x in gene_table_list[0].keys()}
        for gt_d in gene_table_list:
            for k in gene_table_d.keys():
                gene_table_d[k].append(gt_d[k])
   
    # Note this dataframe also contains amino acid sequence of gene
    else:
        raise Exception("No other methods known")

    gene_table_df = pd.DataFrame.from_dict(gene_table_d)
    cols = [
            "locusId", "sysName", "type", "scaffoldId", "begin", "end", 
            "strand", "name", "desc", "GC", "nTA"
            ]
    gene_table_df = gene_table_df[cols]




    return gene_table_df





def parse_CDS_info(CDS_info):
    """
    Args:
        CDS_info (python d):
           'aliases' (list<alias_list (multiple)>):
                alias_list 
                    list<'locus_tag', str> AND/OR 
                    list<'old_locus_tag', str> AND/OR
                    list<'protein_id', str>
            'dna_sequence' (str): The actual DNA sequence
            'functions' (list<str>): First object of list is the function
            'location' (list<scaffold (str), bp (int), strand ("+/-"), length (nt)>)
    Returns:
        gene_table_list_d (dict):
            "locusId":str
            "sysName": ?str
            "type": 1
            "scaffoldId": str
            "begin": int
            "end": int
            "strand": str ("+"/"-")
            "name": str (always "unknown" in this case)
            "desc": str
            "GC": float
            "nTA": int
            "AA_seq": Amino Acid sequence of gene
    """


    gene_table_list_d = {}
    
    #Getting locusId
    aliases_l = CDS_info["aliases"]
    locusId_obj = aliases_l[0]
    if locusId_obj[0] != "locus_tag":
        locus_tag_found = False
        for i in range(1, len(aliases_l)):
            if aliases_l[i][0] == "locus_tag":
                locus_tag_found = True
                locusId_obj = aliases_l[i]
                break
            logging.critical(f"Found locus_tag at different loc of list: {i}")
    else:
        locus_tag_found = True

    if not locus_tag_found:
        raise Exception("Expecting locus_tag from genome object, did not find it.")
    else:
        gene_table_list_d["locusId"] = locusId_obj[1]
        gene_table_list_d["sysName"] =  locusId_obj[1]
    
    # Getting scaffold, location, strand
    scaffold, bp_loc, strand, nt_len  = get_location_info(CDS_info["location"][0])
    gene_table_list_d["scaffoldId"] = scaffold
    gene_table_list_d["begin"] = bp_loc
    gene_table_list_d["end"] = bp_loc + nt_len
    gene_table_list_d["strand"] = strand

    # Getting description 
    gene_table_list_d["desc"] = CDS_info["functions"][0] 


    # Getting GC and nTA
    DNA_seq = CDS_info["dna_sequence"].upper() 
    gene_table_list_d["GC"] = (DNA_seq.count("G") + DNA_seq.count("C"))/float(len(DNA_seq))
    gene_table_list_d["nTA"] = DNA_seq.count("TA")

    
    # Undecidable parts (from the data object)
    gene_table_list_d["type"] = 1
    gene_table_list_d["name"] = "unknown"

    # Adding protein sequence
    gene_table_list_d["AA_seq"] = CDS_info["protein_translation"].upper()

    return gene_table_list_d


def get_location_info(loc_list):
    if len(loc_list) != 4:
        print(loc_list)
        raise Exception("loc_list not as expected ^ ")
    else:
        return loc_list
        

def main():
    # Test this against KBase genome data object
    logging.basicConfig(level=logging.DEBUG)
    args = sys.argv
    inp_json = args[1]
    op_df = args[2]
    result_df = obj_data_to_gene_table(inp_json, cdss_method=True)
    result_df.to_csv(op_df, sep="\t", index=False)


if __name__ == "__main__":
    main()


