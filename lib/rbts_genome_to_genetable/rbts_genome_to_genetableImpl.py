# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

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
    GIT_COMMIT_HASH = "1c2a902daaee74b87f015dc037b8a514cbffa0ab"

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
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_rbts_genome_to_genetable
        logging.basicConfig(level=logging.DEBUG)
        dfu = DataFileUtil(self.callback_url)
        gfu = GenomeFileUtil(self.callback_url)
        # We need the workspace object to get info on the workspace the app is running in.
        token = os.environ.get('KB_AUTH_TOKEN', None)
        ws = Workspace(self.ws_url, token=token)


        genome_ref, output_name, test_bool, upload_bool = validate_params(params)

        # Actual program
        res, res_dir, gene_table_fp = genome_ref_to_gene_table(genome_ref, gfu, self.shared_folder,
                                               ws, params['workspace_name'],
                                               dfu, output_name, test_bool=test_bool,
                                                upload_bool=upload_bool)

        logging.info("Results:")
        # Name, Type, Date
        logging.info(res)


        # Returning file in zipped format:-------------------------------
        file_zip_shock_id = dfu.file_to_shock({'file_path': res_dir,
                                              'pack': 'zip'})['shock_id']

        dir_link = {
                'shock_id': file_zip_shock_id, 
               'name':  'results.zip', 
               'label':'genes_table_output_dir', 
               'description': 'The directory of outputs from running' \
                + ' Genome to Genes Table.'
        }

        report_util = KBaseReport(self.callback_url)


        report_info = report_util.create_extended_report({
                                        'message': "Finished running Genome to Genes Table.",
                                        'file_links': [dir_link],
                                        'workspace_name': params['workspace_name']
                                        })

        logging.info("report_info after creating extended report.")
        logging.info(report_info)

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
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
           success is 0, failure is 1) -> structure: parameter "exit_code" of
           Long
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN genome_to_genetable
        logging.info("Beginning conversion from genome_ref to gene table.")
    
        #if not isinstance(params, dict):
        #    raise Exception(f"params must be 'dict', instead {type(params)}.")
        if "genome_ref" not in params:
            raise Exception("'genome_ref' must be one of params when calling 'genome_to_genetable'. Current params: "
                            ", ".join(params.keys()))

        #logging.basicConfig(level=logging.INFO)
        gfu = GenomeFileUtil(self.callback_url)
        token = os.environ.get('KB_AUTH_TOKEN', None)

        # Gene table written to location tmp_dir/g2gt_results/genes.GC
        # We leave many inputs to this function empty
        # since the upload_bool is FALSE.
        # ( ws, dfu, ws_name, gene_table_name empty )
        res, res_dir, gene_table_fp = genome_ref_to_gene_table(params['genome_ref'],
                                                gfu, 
                                                self.shared_folder,
                                                {}, "",
                                                {}, "",
                                                test_bool=False,
                                                upload_bool=False,
                                                local_func=True)
        
        output = {"exit_code": res}
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
