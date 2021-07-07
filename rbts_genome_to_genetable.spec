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

};
