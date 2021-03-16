import sys
import pandas as pd
from Bio import Entrez
import argparse
import time
import csv

arguments = sys.argv

parser = argparse.ArgumentParser(description='Get RefSeq ID')
parser.add_argument("-i", "--input", required=True, help="Select input file")
parser.add_argument("-o", "--output", required=True, help="Direct the output to a file")
args = parser.parse_args()

input_file = args.input
output_file = args.output

# ---------------------------------------------------------------
# Load GI ids list
# ---------------------------------------------------------------
acc_df = pd.read_csv(input_file, sep='\t', header=None, names=(['acc_id']))

#   Convert Accession IDs column to a list
list_acc = acc_df['acc_id'].tolist()


print(list_acc)
# ---------------------------------------------------------------
# Get GenBank ID - NCBI
# ---------------------------------------------------------------
def get_assembly_id(acc_id):
    """Get Assembly ID from Accession ID"""
    from Bio import Entrez
    handle = Entrez.esearch(db='assembly', term=acc_id)
    record = Entrez.read(handle)
    handle.close()
    assembly = record["IdList"][0]
    return assembly

def get_assembly_summary(acc_id):
    """Get esummary from an assembly ID
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """
    from Bio import Entrez
    #provide your own mail here
    Entrez.email = "carmen.ausejo@sund.ku.dk" #email
    assembly = get_assembly_id(acc_id)
    esummary_handle = Entrez.esummary(db="assembly", id=assembly, report="full")
    esummary_record = Entrez.read(esummary_handle)
    genbank_id = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
    refseq_id = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
    return(assembly, genbank_id, refseq_id)

# ---------------------------------------------------------------
# Get GenBank ID - NCBI
# ---------------------------------------------------------------
genbank_dic = {}

for idx,acc_id in enumerate(list_acc):
    try:
        genbank_dic[acc_id] = (get_assembly_summary(acc_id))
    except:
        print("Accesion ID not found: " + acc_id)
        genbank_dic[acc_id] = ('-', '-', '-')

    print(acc_id)
    print( str(idx+1) + "/" + str(len(list_acc)))

    time.sleep(4)  # Delays for 4 seconds

# Save GenBank IDs in a Df
entrez_df = pd.DataFrame(genbank_dic.items(), columns=['Accession_id','esummary_record'])
entrez_df[['Assembly_id', 'GenBank_id', 'Refseq_id']]= pd.DataFrame(entrez_df.esummary_record.values.tolist(), index= entrez_df.index)
entrez_df = entrez_df[['Accession_id','Assembly_id', 'GenBank_id', 'Refseq_id']]
entrez_df.to_csv(output_file, sep="\t", quoting=csv.QUOTE_NONE, index=False, header=True)
print("Done!")