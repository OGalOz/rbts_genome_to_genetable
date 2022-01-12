
export chome=$PWD
export impl=$PWD/lib/rbts_genome_to_genetable/rbts_genome_to_genetableImpl.py
export tst=$PWD/test/rbts_genome_to_genetable_server_test.py
export tmp_dir=$PWD/test_local/workdir/tmp/
export ui_dir=$PWD/ui/narrative/methods/run_rbts_genome_to_genetable/
export uspec=$PWD/ui/narrative/methods/run_rbts_genome_to_genetable/spec.json
export udisp=$PWD/ui/narrative/methods/run_rbts_genome_to_genetable/display.yaml



#Docker fix
docker run -it -v /var/run/docker.sock:/run/docker.sock alpine chmod g+w /run/docker.sock

# clean up
find . -name '.DS_Store' -type f -delete
