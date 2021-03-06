# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import shutil

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.KBaseReportClient import KBaseReport
from converts.central import genome_ref_to_gene_table, validate_params
#END_HEADER


class rbts_genome_to_genetable:
    '''
    Module Name:
    rbts_genome_to_genetable

    Module Description:
    A KBase module: rbts_genome_to_genetable
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/OGalOz/rbts_genome_to_genetable.git"
    GIT_COMMIT_HASH = "a3b5e2f390d62f264f0457abcca760d0c0b14171"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.ws_url = config['workspace-url']
        #END_CONSTRUCTOR
        pass


    def run_rbts_genome_to_genetable(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_rbts_genome_to_genetable
        res = self.genome_to_genetable(ctx, params)
        output = res[0]
        #END run_rbts_genome_to_genetable

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_rbts_genome_to_genetable return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def genome_to_genetable(self, ctx, params):
        """
        :param params: instance of type "GeneTableParams" (genome_ref is the
           ref of the genome) -> structure: parameter "genome_ref" of String
        :returns: instance of type "GenomeToGeneTableResult" (exit_code (int)
           success is 0, failure is 1 filepath (str) path to where the gene
           table is) -> structure: parameter "exit_code" of Long, parameter
           "filepath" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN genome_to_genetable
        # This is the method that's called by other apps
        logging.info("Beginning conversion from genome_ref to gene table.")


        if "app_test" in params and params["app_test"]:
            x = os.listdir(self.shared_folder)
            if len(x) > 0:
                shutil.rmtree(self.shared_folder)
                logging.info("Cleaning shared folder after multiple tests.")
                os.mkdir(self.shared_folder)

        logging.basicConfig(level=logging.DEBUG)
        dfu = DataFileUtil(self.callback_url)
        gfu = GenomeFileUtil(self.callback_url)
        # We need the workspace object to get info on the workspace the app is running in.
        token = os.environ.get('KB_AUTH_TOKEN', None)
        ws = Workspace(self.ws_url, token=token)
    
        #if not isinstance(params, dict):
        #    raise Exception(f"params must be 'dict', instead {type(params)}.")
        if "genome_ref" not in params:
            raise KeyError("'genome_ref' must be one of params when " + \
                            "calling 'genome_to_genetable'. Current params: "
                            ", ".join(params.keys()))


        # Gene table written to location tmp_dir/g2gt_results/genes.GC
        # We leave many inputs to this function empty
        # since the upload_bool is FALSE.
        # ( ws, dfu, ws_name, gene_table_name empty )
        res, res_dir, gene_table_fp = genome_ref_to_gene_table(params['genome_ref'],
                                                gfu, 
                                                self.shared_folder,
                                                ws, "",
                                                {}, "",
                                                test_bool=False,
                                                upload_bool=False,
                                                local_func=True)
        
        output = {"exit_code": res,
                "filepath": gene_table_fp}
        #END genome_to_genetable

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method genome_to_genetable return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
