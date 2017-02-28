# dowload data
wget http://geneontology.org/gene-associations/goa_human.gaf.gz
wget http://geneontology.org/ontology/go-basic.obo
# go_id name namespace
python3 transformer.py > go_namespace.txt
