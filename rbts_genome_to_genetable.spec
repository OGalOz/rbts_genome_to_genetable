/*
A KBase module: rbts_genome_to_genetable
*/

module rbts_genome_to_genetable {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_rbts_genome_to_genetable(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

    /*
        This function takes a genome ref and creates a gene table
    */

    /*
        genome_ref is the ref of the genome
    */
    typedef structure {
        string genome_ref;
    } GeneTableParams

    /*
        exit_code (int) success is 0, failure is 1
    */
    typedef structure {
        int exit_code;
    } GenomeToGeneTableResult

    funcdef genome_to_genetable(GeneTableParams params) returns (GenomeToGeneTableResults result) authentication required;

};
