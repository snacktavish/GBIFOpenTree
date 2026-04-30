
import json
import csv
import argparse
import dendropy
from opentree import OT

parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
parser.add_argument('-s','--study_id')
#parser.add_argument('-t','--tree_id')


args = parser.parse_args()


def generate_gbif_dwca_from_phylesystem(study_id,
#                                        tree_id,
                                        tip_label='^ot:originalLabel'):
    """
  
    Spaces vs underscores kept being an issue, so all spaces are coerced to underscores when data are read in.
    :param study_id: OpenTree study id of the phylogeny
    :param tree_id: OpenTree tree id of the phylogeny, some studies have several phylogenies.
    :return: Object of class ATT.
    """
    assert(tip_label in ['^ot:originalLabel', 'otu', "^ot:ottTaxonName", "^ot:ottId"])
    try:
        study = OT.get_study(study_id)
        otus = OT.get_otus(study_id)
        trees = [tree['@id'] for tree in study.response_dict['data']['nexml']['trees']['tree']]
        otu_json = otus.response_dict
    except AttributeError:
        sys.stderr.write("Failure getting study {} from OpenTree phylesystem.".format(study_id))
        sys.exit()
    ## this gets the taxa that are in the subtree with all of their info - ott_id, original name,
    assert ( len(otu_json.keys())==1)
    otu_dict = otu_json[list(otu_json.keys())[0]]['otuById']
    original_label_to_otu = {otu_dict[otu_id]['^ot:originalLabel']:otu_id for otu_id in otu_dict} 
    all_otu_keys = []
    for otu_id in otu_dict:
        for subkey in otu_dict[otu_id].keys():
            all_otu_keys.append(subkey)
        otu_dict[otu_id]["dynamicProperties"] = []
        otu_dict[otu_id]["phyloTreeFileName"] = []
    all_otu_keys =  set(all_otu_keys)
    for tree_id in trees:
        tree_reference = str(study_id)+"_"+str(tree_id)
        tree_response = OT.get_tree(study_id, tree_id, tree_format="nexus")
        tree_obj = dendropy.Tree.get(data=tree_response.response_dict['content'].decode('utf-8'), schema="nexus")
        phylo_treefile_name = tree_reference+".nex"
        tree_obj.write(path=phylo_treefile_name, schema="nexus")
        treed_taxa = {}
        for leaf in tree_obj.leaf_node_iter():
            tn = leaf.taxon
            otu_id = original_label_to_otu[tn.label]
            otu_dict[otu_id]["otuId"] = otu_id
            otu_dict[otu_id]["phyloTreeTipLabel"] = otu_dict[otu_id].get('^ot:originalLabel',"")
            otu_dict[otu_id]["phyloTreeFileName"].append(phylo_treefile_name)
            otu_dict[otu_id]["scientificName"] = otu_dict[otu_id].get('^ot:ottTaxonName',"")
            otu_dict[otu_id]["dynamicProperties"].append({"phyloTreeTipLabel":otu_dict[otu_id].get('^ot:originalLabel',""),"phyloTreeFileName":phylo_treefile_name}) 
    header = ["otuId","scientificName","phyloTreeTipLabel","phyloTreeFileName","dynamicProperties"] + list(all_otu_keys)

    with open(study_id+"occurrence.csv", "w") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for otu in otu_dict:
            data = [otu_dict[otu].get(value, "") for value in header]
            writer.writerow(data)


                                         

generate_gbif_dwca_from_phylesystem(args.study_id)






"""
def bulk_tnrs_load(filename):
    "Read in outputs from OpenTree Bulk TNRS,
    translates to a Physcraper OTU dictionary.
    :param filename: input json file
    "
    otu_dict = {}
    with open(filename) as data_file:
        input_dict = json.load(data_file)
    if "names" in input_dict.keys():
        for name in input_dict["names"]:
            i = 1
            otu = "otu" + name['id'].strip('name')
            while otu in otu_dict.keys():
                otu = "{}_{}".format(otu, i)
                i += 1
            otu_dict[otu] = {"^ot:originalLabel":name["originalLabel"]}
            if name.get("ottTaxonName"):
                otu_dict[otu]["^ot:ottTaxonName"] = name["ottTaxonName"]
            if name.get("ottId"):
                otu_dict[otu]["^ot:ottId"] = name["ottId"]
            for source in name.get("taxonomicSources", []):
                #debug(source)
                if source:
                    taxsrc = source.split(":")
                    assert len(taxsrc) == 2, taxsrc
                    otu_dict[otu]["^{}:taxon".format(taxsrc[0])] = taxsrc[1]
        for otu in otu_dict:
            otu_dict[otu]["^physcraper:status"] = "original"
            otu_dict[otu]["^physcraper:last_blasted"] = None
            otu_dict[otu]["^physcraper:ingroup"] = "unknown"
    else:
        for otu in input_dict:
            assert input_dict[otu]["^physcraper:status"], otu_dict[otu]
            otu_dict = input_dict
    return otu_dict
"""