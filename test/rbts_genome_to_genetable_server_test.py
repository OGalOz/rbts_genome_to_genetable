# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from rbts_genome_to_genetable.rbts_genome_to_genetableImpl import rbts_genome_to_genetable
from rbts_genome_to_genetable.rbts_genome_to_genetableServer import MethodContext
from rbts_genome_to_genetable.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class rbts_genome_to_genetableTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('rbts_genome_to_genetable'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'rbts_genome_to_genetable',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = rbts_genome_to_genetable(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_RBTS_Genome_to_Genes_Table_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa
        cls.only_first = False

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    """
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods


        # We add 'test_run' key to params in order to get the right workspace ID
        test_params = {
                "workspace_name": self.wsName,
                "genome_ref": "62572/2/1",
                "output_name": "Test_Genes_Table2",
                "test_run": 1,
                "upload_bool": True 
        }

        ret = self.serviceImpl.run_rbts_genome_to_genetable(self.ctx, test_params)
    """
    def test_app_function(self):
        test_params = {
                "genome_ref": "66889/5/1",
                "workspace_name": self.wsName,
                "app_test": True
        }

        ret = self.serviceImpl.genome_to_genetable(self.ctx, test_params)

    def test_local_function3(self):
        if not self.only_first:
            test_params = {
                    "genome_ref": "66889/3/1",
                    "workspace_name": self.wsName,
                    "app_test": True
            }

            ret = self.serviceImpl.genome_to_genetable(self.ctx, test_params)

    def test_local_function1(self):
        if not self.only_first:
            test_params = {
                    "genome_ref": "62686/4/1",
                    "workspace_name": self.wsName,
                    "app_test": True
            }

            ret = self.serviceImpl.genome_to_genetable(self.ctx, test_params)
    def test_local_function2(self):
        if not self.only_first:
            test_params = {
                    "genome_ref": "63063/9/1",
                    "workspace_name": self.wsName,
                    "app_test": True
            }
            ret = self.serviceImpl.genome_to_genetable(self.ctx, test_params)
    def test_multiple_genomes(self):
        if not self.only_first:
            genome_refs = ["62686/4/1", "63063/9/1"]

            for i in range(len(genome_refs)):
                gnm_ref = genome_refs[i]
                test_params = {
                        "workspace_name": self.wsName,
                        "genome_ref": gnm_ref,
                        "output_name": "BASE_TEST_" + gnm_ref.replace("/","_"),
                        "app_test": True,
                        "test_num": str(i+1)
                }
                ret = self.serviceImpl.run_rbts_genome_to_genetable(self.ctx, test_params)
