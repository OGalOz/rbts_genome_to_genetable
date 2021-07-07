#python3
"""
Central place to both download genbank and genome fna
and convert them into a gene table file with all columns,
including the GC content and nTA.
We additionally upload the Gene Table data object.
"""

import os
from Bio import SeqIO
import logging
import datetime
from central.genbank_to_gene_table import genbank_and_genome_fna_to_gene_table



def genome_ref_to_gene_table(genome_ref, gfu, tmp_dir,
                             ws, dfu, gene_table_name,
                             test_bool=False):
    """
    Args:
        genome_ref (str): 
        gfu: GenomeFileUtil Class Object
        dfu: DataFileUtil Class Object
        ws: Workspace Object
        tmp_dir (str): Path to directory where work is done
        gene_table_name (str): Output name of gene_table
        test_bool (bool): Whether we should get the workspace
                          from the genome ref or from the workspace.

    Returns:
        genes_GC_fp (str): Path to location of genes table

    Description:
        We use the GenomeFileUtil Object to download the genome
        information, like the FNA file and the genbank file
        to the location 'tmp_dir', then we create a gene table
        with the function 'genbank_and_genome_fna_to_gene_table'.
        The file is created at the location 'tmp_dir/genes.GC'
        the header columns will always be the same
    """

    # Download genome in GBK format and convert it to fna:
    # gt stands for genome table
    genome_fna_fp, gbk_fp, ws_id = DownloadGenomeToFNA(
            gfu, genome_ref, tmp_dir, test_bool)
    genome_scientific_name = GetGenomeOrganismName(ws, genome_ref)

    gene_table_fp = os.path.join(tmp_dir, "genes.GC")
    # This function creates the gene_table at the location gene_table_fp
    num_lines = genbank_and_genome_fna_to_gene_table(gbk_fp, genome_fna_fp, gene_table_fp)
    
    res = upload_gene_table_object_to_KBase(gene_table_fp, dfu, ws,
                                            num_lines, 
                                            genome_ref, organism_scientific_name,
                                            gene_table_name,
                                            ws_id=ws_id)



# We download Genome Files: gfu is Genome File Util
def DownloadGenomeToFNA(gfu, genome_ref, scratch_dir, test_bool):
    """
    Args: 
        GFU Object, str (A/B/C), str path
        test_bool: Get workspace from genome object.
    Outputs: [fna_fp (str), gbk_fp (str)]
    """

    GenomeToGenbankResult = gfu.genome_to_genbank({'genome_ref': genome_ref})

    logging.info(GenomeToGenbankResult)

    raise Exception("Defined stop - ")

    genbank_fp = GenomeToGenbankResult['genbank_file']['file_path']

    genome_fna_fp = get_fa_from_scratch(scratch_dir)


    return [genome_fna_fp, genbank_fp]


def get_fa_from_scratch(scratch_dir):
    """
    Careful... May not work in the Future
    Inputs:
        scratch_dir: (str) Path to work dir/ tmp etc..
    Outputs:
        FNA fp: (str) Automatic download through GenbankToGenome
    """
    
    scratch_files = os.listdir(scratch_dir)
    fna_paths = []
    for f in scratch_files:
        if f[-2:] == "fa":
            fna_paths.append(os.path.join(scratch_dir, f))
   
    if len(fna_paths) == 0:
        raise Exception("Could not find '.fa' file in scratch directory, files "
                        "in directory are " + ", ".join(scratch_files))
    elif len(fna_paths) > 1:
        raise Exception("Found multiple '.fa' files in scratch directory, files "
                        "in directory are " + ", ".join(scratch_files) + "." + \
                        " the '.fa' files found are: " + ", ".join(fna_paths) + "." + \
                        " Only expecting a single '.fa' file from the genome object. Cannot continue.")
    else:
        fna_fp = fna_paths[0]

    return fna_fp


def GetGenomeOrganismName(ws, genome_ref):
    """
    # Getting the organism name using WorkspaceClient
    ws: workspace client object
    genome_ref: (str) A/B/C
    """
    res = ws.get_objects2(
        {
            "objects": [
                {
                    "ref": genome_ref,
                    "included": ["scientific_name"],
                }
            ]
        }
    )
    scientific_name = res["data"][0]["data"]["scientific_name"]
    return scientific_name




def genbank_and_genome_fna_to_gene_table(gbk_fp, gnm_fp, op_fp):
    """
    Args:
        gbk_fp: (str) Path to genbank file.
        gnm_fp: (str) Path to genome fna file
        op_fp: (str) Path to write genes table to

    We use GenBank Records and Features

    """

    
    id2seq = parseFASTA(gnm_fp)

    out_FH = open(op_fp, 'w')
    # This is the output file start line:
    out_FH.write("locusId\tsysName\ttype\t" + \
                    "scaffoldId\tbegin\tend\tstrand\t" + \
                    "name\tdesc\tGC\tnTA\n")
    
    gb_record_generator = SeqIO.parse(open(gbk_fp, "r"), "genbank")

    # TYPE - Note that we don't like misc_feature or gene
    # May need to skip anything besides values under 10
    types_dict = {"CDS" : 1, "rRNA" : 2, "tRNA" : 5, 
                   "RNA" : 6, "transcript" : 6,
                   "pseudogene": 7, "misc_feature": 20, 
                   "gene": 21}

    for gb_record in gb_record_generator:

        locus_tag = gb_record.name
        #Genome sequence:
        g_seq = str(gb_record.seq)

        scaffold = findSeqInId2Seq(g_seq, id2seq)

        g_len = len(g_seq)

        #Genome features (list of features):
        g_features = gb_record.features
        #DEBUG
        #print(g_features[0])
        g_feat_len = len(g_features)

        """
        scaffoldId_exists= False
        if "scaffold_name" in config_dict:
            scaffold_id = config_dict["scaffold_name"]
            scaffoldId_exists = True
        """

        try:
            for i in range(g_feat_len):


                current_feat = g_features[i]

                if current_feat.type == "source":
                    print("Feature is of type 'source':")
                    print(current_feat)
                    continue

                
                #scaffold already set above
                begin = "null"
                end = "null"
                strand = "null"
                desc = "null"
                typ = "null"
                #locus_tag = "null"
                sysName = "null"
                name = "null"
                GC = "null"
                nTA = "null"


                """
                # Scaffold Id
                if scaffoldId_exists:
                    if scaffold_id in g_features[i].qualifiers:
                        scaffold = g_features[i].qualifiers[scaffold_id]
                    else:
                        logging.debug("Could not find scaffold id "
                                "{} in qualifiers:".format(scaffold_id))
                        logging.debug(g_features[i].qualifiers)
                        scaffold = "1"
                else:
                    scaffold = "1"

                """

                # Begin
                begin = str(current_feat.location.start + 1)
                # End
                end = str(current_feat.location.end + 1)

                # Strand
                if current_feat.strand == 1:
                    strand = "+"
                elif current_feat.strand == -1:
                    strand = "-"
                else:
                    logging.critical("Could not recognize strand type.")
                    raise Exception("Parsing strand failed.")

                # Desc (Description)
                if "product" in current_feat.qualifiers.keys():
                    desc = str(current_feat.qualifiers['product'][0])
                else:
                    desc = current_feat.type
                    logging.critical("Could not find description in current_feat: ")
                    logging.critical(current_feat)
                    #continue

                typ_str = current_feat.type.strip()
                if typ_str in types_dict:
                    typ = str(types_dict[typ_str])
                else:
                    logging.info("Could not recognize type from feature: " \
                            + typ_str)
                    typ = "0"
                if typ == "1":
                    out_FH.write("\t".join([locus_tag, sysName, typ, scaffold,
                            begin, end, strand, name, desc, GC, nTA]) + "\n")
                else:
                    logging.info(f"Did not write annotation for gene with type {typ}")
        except:
            logging.critical("Could not parse all features in genbank file.")
            out_FH.close()
            raise Exception("Parsing genbank file into gene table failed")
    
    out_FH.close()

    return None


def upload_gene_table_object_to_KBase(gene_table_fp, dfu, ws,
                                      num_lines,
                                      genome_ref, organism_scientific_name,
                                      gene_table_name,
                                      ws_id=None):
    """
    Args:
        gene_table_fp (str) Path to gene table
        dfu: DataFileUtil object
        num_lines (int): Number of lines in the file at gene_table_fp
        genome_ref (str): The permanent id of the genome object from which this is taken.
        organism_scientific_name (str): The name of the organism related to the genome.
        gene_table_name (str): Output name
        ws_id (str or None): If ws_id is None, then we proceed as normal, if it is a string then
                            we use that.
    """
    # We create the handle for the object:
    file_to_shock_result = dfu.file_to_shock(
        {"file_path": gene_table_fp, "make_handle": True, "pack": "gzip"}
    )
    # The following var res_handle only created for simplification of code
    res_handle = file_to_shock_result["handle"]

    # We create a better Description by adding date time and username
    date_time = datetime.datetime.utcnow()
    #new_desc = "Uploaded by {} on (UTC) {} using Uploader. User Desc: ".format(
    #        self.params['username'], str(date_time))
    column_header_list = "locusId,sysName,type," + \
                    "scaffoldId,begin,end,strand," + \
                    "name,desc,GC,nTA".split(",")


    # We create the data for the object
    genes_data = {
        "file_type": "KBaseRBTnSeq.RBTS_InputGenesTable",
        "input_genes_table": res_handle["hid"],
        # below should be shock
        "handle_type": res_handle["type"],
        "shock_url": res_handle["url"],
        "shock_node_id": res_handle["id"],
        "compression_type": "gzip",
        "file_name": res_handle["file_name"],
        "utc_created": str(date_time),
        "column_header_list": column_header_list,
        "column_headers_str": ", ".join(column_header_list),
        "num_lines": str(num_lines),
        "related_genome_ref": genome_ref,
        "related_organism_scientific_name": organism_scientific_name
    }

    if ws_id is None:
        ws_info = ws.get_workspace_info({'workspace': params['workspace_name']})
        ws_id = ws_info[0]


    save_object_params = {
        "id": ws_id,
        "objects": [
            {
                "type": "KBaseRBTnSeq.RBTS_InputGenesTable",
                "data": genes_data,
                "name": genetable_name,
            }
        ],
    }
    # save_objects returns a list of object_infos
    dfu_object_info = dfu.save_objects(save_object_params)[0]
    print("dfu_object_info: ")
    print(dfu_object_info)
    return {
        "Name": dfu_object_info[1],
        "Type": dfu_object_info[2],
        "Date": dfu_object_info[3],
    }


def validate_params(params):
    """
    Args:
        params (d):
            genome_ref (str)
            output_name (str)
    """
    for x in ["genome_ref", "output_name"]:
        if x not in params:
            raise Exception(f"Expecting parameter {x} as an input, but not found. " + ", ".join(params.keys()))

    if len(params["genome_ref"].split("/")) != 3:
        raise Exception(f"Expecting genome ref in format 'A/B/C', instead got {params['genome_ref']}")

    if " " in params['output_name']:
        raise Exception(f"Output name cannot contain spaces. Output name: {params['output_name']}")

    test_bool = False
    if 'test_run' in params:
        test_bool=True

    return [params["genome_ref"], params["output_name"], test_bool]




