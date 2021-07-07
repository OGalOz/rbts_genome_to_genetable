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
from converts.genbank_to_gene_table import genbank_and_genome_fna_to_gene_table 



def genome_ref_to_gene_table(genome_ref, gfu, tmp_dir,
                             ws, ws_name,
                             dfu, gene_table_name,
                             test_bool=False,
                             upload_bool=True):
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
        upload_bool (bool): Whether we should upload the object to KBase 

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
    genome_fna_fp, gbk_fp = DownloadGenomeToFNA(gfu, genome_ref, tmp_dir)
    genome_scientific_name, ws_id = GetGenomeOrganismName(ws, genome_ref, test_bool)

    res_dir = os.path.join(tmp_dir, "results")
    os.mkdir(res_dir)
    gene_table_fp = os.path.join(res_dir, "genes.GC")
    # This function creates the gene_table at the location gene_table_fp
    num_lines = genbank_and_genome_fna_to_gene_table(gbk_fp, genome_fna_fp, gene_table_fp)
   
    if upload_bool:
        res = upload_gene_table_object_to_KBase(gene_table_fp, dfu, ws, ws_name,
                                            num_lines, 
                                            genome_ref, genome_scientific_name,
                                            gene_table_name,
                                            ws_id=ws_id)
    else:
        res = {}

    return [res, res_dir]



# We download Genome Files: gfu is Genome File Util
def DownloadGenomeToFNA(gfu, genome_ref, scratch_dir):
    """
    Args: 
        GFU Object, str (A/B/C), str path
    Outputs: [fna_fp (str), gbk_fp (str)]
    """

    GenomeToGenbankResult = gfu.genome_to_genbank({'genome_ref': genome_ref})

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


def GetGenomeOrganismName(ws, genome_ref, test_bool):
    """
    Args:
        ws: workspace client object
        genome_ref: (str) A/B/C
        test_bool: Get workspace from genome object.
    Returns:
        scientific_name (str)
        ws_id (int)
    Description:
        Getting the organism name using WorkspaceClient

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
    ws_id = None
    if test_bool:
        ws_id = res["data"][0]["info"][6]
    return scientific_name, ws_id



def upload_gene_table_object_to_KBase(gene_table_fp, dfu, ws,
                                      ws_name,
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
        ws_id (int or None): If ws_id is None, then we proceed as normal, if it is an int then
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
    column_headers_str = "locusId, sysName, type, scaffoldId, begin, end, strand, name, desc, GC, nTA"
    column_header_list = column_headers_str.split(', ')


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
        "column_headers_str": column_headers_str,
        "num_lines": str(num_lines),
        "related_genome_ref": genome_ref,
        "related_organism_scientific_name": organism_scientific_name
    }

    if ws_id is None:
        ws_info = ws.get_workspace_info({'workspace': ws_name})
        ws_id = ws_info[0]


    save_object_params = {
        "id": ws_id,
        "objects": [
            {
                "type": "KBaseRBTnSeq.RBTS_InputGenesTable",
                "data": genes_data,
                "name": gene_table_name,
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

    upload_bool = True
    if 'upload_bool' in params:
        upload_bool = params['upload_bool']

    return [params["genome_ref"], params["output_name"], test_bool, upload_bool]


