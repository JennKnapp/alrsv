# -*- coding: utf-8 -*-
"""
Created on Thu May 25 2023

@author: Jenn Knapp
email: jknapp@uwaterloo.ca
"""
"""
Purpose: Retrieve clade-definining mutations from any nextstrain phylogenetic tree
    and create clade definition .json files 
Reguires:
phylogenetic_tree.json (downloaded from nextstrain.org)

rsv_a_genome.json
clades/

"""
import json
import os

def process_phylogenetic_tree(tree):
    """
    Process the phylogenetic tree and generate JSON files for each clade.
    """
    print("Processing phylogenetic tree...")
    clades = set()
    traverse_tree(tree["tree"], clades)  # Traverse the tree to collect unique clade names
    
    for clade_name in clades:
        clade_mutations = extract_all_mutations(tree["tree"], clade_name)  # Extract mutations for the clade
        write_json_file(clade_name, clade_mutations)  # Write mutations to a JSON file
        print(f"Processed clade: {clade_name}")
    
    print("Processing completed.")

def traverse_tree(node, clades, target_clade=None, parent_mutations=None):
    """
    Traverse the phylogenetic tree to collect unique clade names.
    If a target clade is provided, only mutations from that clade and its parent nodes are collected.
    """
    if "node_attrs" in node and "clade_membership" in node["node_attrs"]:
        clade_name = node["node_attrs"]["clade_membership"]["value"]
        if clade_name == target_clade or parent_mutations is not None:
            clades.add(clade_name)  # Add clade name to the set

    if "branch_attrs" in node and "mutations" in node["branch_attrs"]:
        branch_mutations = node["branch_attrs"]["mutations"]
        if "nuc" in branch_mutations:
            if parent_mutations is not None:
                parent_mutations.update(branch_mutations["nuc"])  # Update parent mutations with current mutations
            else:
                parent_mutations = set(branch_mutations["nuc"])

    if "children" in node:
        for child in node["children"]:
            traverse_tree(child, clades, target_clade, parent_mutations)  # Recursively traverse child nodes

def extract_all_mutations(node, target_clade=None):
    """
    Extract all mutations from a specific clade and its parent nodes.
    If a target clade is provided, only mutations from that clade are extracted.
    Exclude mutations with positions less than 1.
    """
    mutations = set()

    if target_clade and "node_attrs" in node and "clade_membership" in node["node_attrs"]:
        clade_membership = node["node_attrs"]["clade_membership"]["value"]
        if clade_membership == target_clade:
            if "branch_attrs" in node and "mutations" in node["branch_attrs"]:
                branch_mutations = node["branch_attrs"]["mutations"]
                if "nuc" in branch_mutations:
                    for mutation in branch_mutations["nuc"]:
                       # ref_nt = mutation[0]
                       # position = int(mutation[1:-1]) # - 106 # for ToBRFV
                       # alt_nt = mutation[-1]
                       # if position >= 1:
                       #     new_mutation = f"{ref_nt}{position}{alt_nt}"
                        mutations.add(mutation)
            return list(mutations)  # Return mutations if target clade is found

    if "children" in node:
        for child in node["children"]:
            child_mutations = extract_all_mutations(child, target_clade)  # Extract mutations from child nodes
            mutations.update(child_mutations)  # Add child mutations to the set

    return list(mutations)

def write_json_file(clade_name, mutations):
    """
    Write clade mutations to a JSON file.
    """
    data = {
        "label": clade_name,
        "description": f"{clade_name} defining mutations",
        "sources": [],
        "tags": [clade_name],
        "sites": mutations,
            }

    output_dir = "clades"
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, f"{clade_name}.json")
    try:
        with open(filename, "w") as json_file:
            json.dump(data, json_file, indent=4)
        print(f"Generated JSON file: {filename}")
    except IOError:
        print(f"Error writing JSON file: {filename}")

# Load the input phylogenetic tree from a JSON file
with open("rsv_a_genome.json", "r") as json_file:
    nextstrain_data = json.load(json_file)

# Process the phylogenetic tree
process_phylogenetic_tree(nextstrain_data)

