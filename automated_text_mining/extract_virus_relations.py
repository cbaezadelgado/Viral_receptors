#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 13:01:07 2022

@author: Carlos Baeza-Delgado
"""

# Import modules

import os
import sys
import json
import subprocess
from datetime import date
from time import sleep
import yaml
from Bio import Entrez
from alive_progress import alive_bar

# Run from this script's directory so relative paths (scripts/, BERN2/,
# spacy/, openie/) resolve wherever the repository is installed
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from scripts.preprocess_abstract import parse_abstracts
from scripts.preprocess_abstract import get_batches_abstracts
from scripts.preprocess_abstract import get_annotated_abstracts
from scripts.preprocess_abstract import get_annotated_abstracts_parallel
from scripts.entrez_functions import load_json_database
from scripts.entrez_functions import update_json_database
from scripts.entrez_functions import extract_annotations_and_ids
from scripts.entrez_functions import extract_sp_ids
from scripts.entrez_functions import extract_gene_ids
# from scripts.entrez_functions import extract_sp_and_gene_ids
from scripts.entrez_functions import download_genes
from scripts.entrez_functions import download_species
from scripts.entrez_functions import download_abstracts
# from scripts.process_mentions import parse_all_mentions
from scripts.process_mentions import parse_all_mentions_parallel
from scripts.process_mentions import filter_mentions
from scripts.get_extracted_relations import get_relations


# Load config
with open('config.yaml',
          mode='r',
          encoding='utf-8') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

PATH = config['path']
DB_PATH = config['db_path']
PMIDS_FILE = config['pmids']
output_tag = config['output_tag']
GPU = config['GPU']
N_PROC_L = config['n_processes'][0]
N_PROC_H = config['n_processes'][1]
Entrez.email = config['entrez_mail']
Entrez.api_key = config['entrez_api']
Entrez.max_tries = config['entrez_maxtries']
Entrez.sleep_between_tries = config['entrez_sleep_tries']
relations = config['relations']

if len(relations) == 0:
    print('\n\n WARNING: a least one valid relations must be included in ' +
          'config file\n\n')
    sys.exit()

if not os.path.exists(PATH):
    print(f"\n\n WARNING: directory '{PATH}' NOT FOUND\n\n")
    sys.exit()

DATE = date.today().strftime('%d-%m-%Y')
config['date'] = DATE
stats = []

# Create temp directory
if not os.path.exists(os.path.join(PATH, 'temp')):
    os.mkdir(os.path.join(PATH, 'temp'))


# Check previous files
mentions_to_parse = []
for key in ['parse_species', 'parse_disease', 'parse_genes']:
    if config[key]:
        mentions_to_parse.append(key.split('_')[1])

if config['add_date']:
    suff = f'{output_tag}_{DATE}'
else:
    suff = f'{output_tag}'

config_file_temp = f'config_{suff}.yaml'

if os.path.exists(os.path.join(PATH, config_file_temp)):
    EXISTS = True
    N_FILE = 1
    while EXISTS:
        config_file_temp = f'config_{suff}_{N_FILE}.yaml'
        if os.path.exists(os.path.join(PATH, 'temp', config_file_temp)):
            N_FILE += 1
            config_file_temp = f'config_{suff}_{N_FILE}.yaml'
        else:
            N_FILE = f'_{N_FILE}'
            EXISTS = False
else:
    N_FILE = ''

suffix = f'{suff}{N_FILE}'


# Save config file
config_file = f'config_{suffix}.yaml'
with open(os.path.join(PATH, config_file),
          mode='w',
          encoding='utf-8') as file:
    documents = yaml.dump(config, file)


# Create temp directory
if not os.path.exists(os.path.join(PATH, 'temp')):
    os.mkdir(os.path.join(PATH, 'temp'))


# Load file with pmids
with open(os.path.join(PATH, PMIDS_FILE),
          mode='r',
          encoding='utf-8') as f:
    # pmids = list(set(f.read().split('\n')))
    pmids = [x.strip() for x in f.readlines()]

pmids = list(set(pmids))

print(f'\n- FILE {PMIDS_FILE} with {len(set(pmids))} unique PMIDs ' +
      'LOADED successfully.')


# Load databases
abstracts_raw = load_json_database(DB_PATH,
                                   'abstracts.json')
stored_annotated_abstracts = load_json_database(DB_PATH,
                                                'annotated_abstracts.json')
annotated_parsed_abstracts = load_json_database(DB_PATH,
                                              'annotated_parsed_abstracts.json')
entrez_species = load_json_database(DB_PATH,
                                    'entrez_species.json')
entrez_genes = load_json_database(DB_PATH,
                                  'entrez_genes.json')
species_mentions = load_json_database(DB_PATH,
                                      'species_mentions.json')
genes_mentions = load_json_database(DB_PATH,
                                    'genes_mentions.json')

# Load alias-virus file
with open(os.path.join(DB_PATH, 'alias_virus.json'),
          mode='r',
          encoding='utf-8') as json_file:
    alias_virus = json.load(json_file)


# Download abstracts from NCBI
pmids_to_download = [pmid for pmid in pmids
                     if pmid not in abstracts_raw]


if len(pmids_to_download) > 0:
    print(f'\n- {len(pmids_to_download)} PMIDs not stored in database, ' +
          'proceeding to download them:')
    abstracts_list = download_abstracts(pmids_to_download,
                                        n_processes=N_PROC_H)

    # Remove empty abstracts from list
    for a in abstracts_list[:]:
        if 'abstract' not in a:
            abstracts_list.remove(a)
        elif a['abstract'] == 'None':
            abstracts_list.remove(a)

    abstracts = {}
    for a in abstracts_list:
        abstracts[a['pmid']] = a
        if 'sent_id' in a:
            del abstracts[a['pmid']]['sent_id']

    abstracts_raw.update(abstracts)
    update_json_database(DB_PATH, 'abstracts.json', abstracts_raw)

    # Get pmids not download and save file
    failed_pmids = {pmid for pmid in pmids if pmid not in abstracts_raw}
    if len(failed_pmids) > 0:
        with open(os.path.join(PATH, f'failed_pmids_{suffix}.txt'),
                  mode='w',
                  encoding='utf-8') as f:
            for bad_pmid in failed_pmids:
                f.write(f'{bad_pmid}\n')
    if len(failed_pmids) > 0:
        print(f'\n- {len(failed_pmids)} PMIDs failed.')

    raw_pmids = set(abstracts_raw.keys())
    annotated_pmids = set(stored_annotated_abstracts.keys())
    pmids_to_annotate = raw_pmids.difference(annotated_pmids)
    
    abstracts_to_annotate = [v for k, v in abstracts_raw.items() if k in
                             pmids_to_annotate]

    # Parse abstracts
    if len(abstracts_to_annotate) > 0:
        print(f'\n\n- {len(abstracts_to_annotate)} abstracts not annotated, ' +
                  "let's annotate them all!!")
        parsed_abstracts = parse_abstracts(abstracts_to_annotate)
        batches = get_batches_abstracts(parsed_abstracts,
                                        size=500)
        # Get mentions from BERN2

        print('\n\n- Launching BERN2 app')
        if GPU:
            _ = subprocess.Popen(['bash', 'scripts/launch_bern2_gpu.sh'],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            with alive_bar(100) as bar:
                for i in range(100):
                    sleep(0.25)
                    bar()
        else:
            _ = subprocess.Popen(['bash', 'scripts/launch_bern2_cpu.sh'],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            with alive_bar(100) as bar:
                for i in range(100):
                    sleep(0.25)
                    bar()

        subprocess.call(['bash',
                         'empty_directories.sh'],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
        for n_batch, batch in enumerate(batches):
            try:
                ann_abs_new = get_annotated_abstracts_parallel(batch,
                                                               n_batch+1,
                                                               len(batches),
                                                               n_process=3)
                # ann_abs_new = get_annotated_abstracts(batch,
                #                                       n_batch+1,
                #                                       len(batches))
                # sleep(1)
                subprocess.call(['bash',
                                 'scripts/empty_directories.sh'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            except Exception:
                if GPU:
                    _ = subprocess.Popen(['bash',
                                          'scripts/close_bern2_gpu.sh'],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
                else:
                    _ = subprocess.Popen(['bash',
                                          'scripts/close_bern2_cpu.sh'],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)

            # Update annotated abstracts database
    
            annotated_abstracts_new_dict = {}
            for sent in ann_abs_new:
                if sent['pmid'] in annotated_abstracts_new_dict:
                    annotated_abstracts_new_dict[sent['pmid']].append(sent)
                else:
                    annotated_abstracts_new_dict[sent['pmid']] = [sent]

            stored_annotated_abstracts.update(annotated_abstracts_new_dict)
            update_json_database(DB_PATH,
                                 'annotated_abstracts.json',
                                 stored_annotated_abstracts)
        if GPU:
                _ = subprocess.Popen(['bash', 'scripts/close_bern2_gpu.sh'],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
        else:
            _ = subprocess.Popen(['bash', 'scripts/close_bern2_cpu.sh'],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    else:
        print('\n\n- All PMIDS found in database, retrieving annotated abstracts\n')

else:
    print('\n\n- All PMIDs found in database, retrieving annotated abstracts\n')
    failed_pmids = set()


# annotated_abstracts = [x for pmid, sents in stored_annotated_abstracts.items()
#                        for x in sents if pmid in pmids]

annotated_abstracts = []
for pmid, sents in stored_annotated_abstracts.items():
    if pmid in pmids:
        for sent in sents:
            annotated_abstracts.append(sent)


for i, sent in enumerate(annotated_abstracts):
    sent['sent_id'] = i

u_pmids = {x['pmid'] for x in annotated_abstracts}
stat1 = f'- {len(annotated_abstracts)} sentences annotated from ' +\
        f'{len(u_pmids)} unique PMIDS ({len(failed_pmids)} failed).'
stats.append(stat1)

abstracts_to_parse = []
for pmid, sents in stored_annotated_abstracts.items():
    if pmid not in annotated_parsed_abstracts:
        for sent in sents:
            abstracts_to_parse.append(sent)
# Parse mentions

print()
(sp_ids,
 gene_ids,
 sp_mentions_new,
 gene_mentions_new) = extract_annotations_and_ids(abstracts_to_parse)

sp_ids = sp_ids.difference(set(entrez_species.keys()))
gene_ids = gene_ids.difference(set(entrez_genes.keys()))
sp_mentions_new = sp_mentions_new.difference(set(species_mentions.keys()))
gene_mentions_new = gene_mentions_new.difference(set(genes_mentions.keys()))


if 'species' in mentions_to_parse:
    if len(sp_mentions_new) > 0:
        sp_ids, species_mentions_temp = extract_sp_ids(sp_ids,
                                                       sp_mentions_new,
                                                       alias_virus)
        species_mentions.update(species_mentions_temp)
        update_json_database(DB_PATH, 'species_mentions.json', species_mentions)

    if len(sp_ids) > 0:
        entrez_species_new = download_species(list(sp_ids),
                                              size=200,
                                              n_processes=N_PROC_H)

        entrez_species.update(entrez_species_new)
        update_json_database(DB_PATH, 'entrez_species.json', entrez_species)


if 'genes' in mentions_to_parse:
    if len(gene_mentions_new) > 0:
        # mentions_batches = get_batches_abstracts(list(gene_mentions_new))
        # for i, batch in enumerate(mentions_batches):
            # print(f'\n- BATCH {i}/{len(mentions_batches)}.', end=' ')
        gene_ids, gene_mentions_temp = extract_gene_ids(gene_ids,
                                                        gene_mentions_new)
        genes_mentions.update(gene_mentions_temp)
        update_json_database(DB_PATH, 'genes_mentions.json', genes_mentions)

    if len(gene_ids) > 0:
        entrez_genes_new = download_genes(list(gene_ids),
                                          size=20,
                                          n_processes=N_PROC_H)

        entrez_genes.update(entrez_genes_new)
        update_json_database(DB_PATH, 'entrez_genes.json', entrez_genes)


# Update databases
update_json_database(DB_PATH, 'entrez_species.json', entrez_species)
update_json_database(DB_PATH, 'entrez_genes.json', entrez_genes)
update_json_database(DB_PATH, 'species_mentions.json', species_mentions)
update_json_database(DB_PATH, 'genes_mentions.json', genes_mentions)
sleep(2)

# Parse mentions
print()
batches = get_batches_abstracts(abstracts_to_parse,
                                size=1000)
for n_batch, batch in enumerate(batches):
    annotated_parsed_batch = parse_all_mentions_parallel(batch,
                                                         alias_virus,
                                                         species_mentions,
                                                         genes_mentions,
                                                         entrez_species,
                                                         entrez_genes,
                                                         mentions_to_parse,
                                                         n_batch,
                                                         len(batches),
                                                         n_proc=20)
    annotated_parsed_to_store = {}
    for sent in annotated_parsed_batch:
        if sent['pmid'] in annotated_parsed_to_store:
            annotated_parsed_to_store[sent['pmid']].append(sent)
        else:
            annotated_parsed_to_store[sent['pmid']] = [sent]
    
    annotated_parsed_abstracts.update(annotated_parsed_to_store)
    update_json_database(DB_PATH,
                         'annotated_parsed_abstracts.json',
                         annotated_parsed_abstracts)

annotated_parsed = []
for pmid, sents in annotated_parsed_abstracts.items():
    if pmid in pmids:
        for sent in sents:
            annotated_parsed.append(sent)

for i, sent in enumerate(annotated_parsed):
    sent['sent_id'] = i

# Filter mentions
print()
annotated_mentions, _ = filter_mentions(annotated_parsed,
                                        relations)
pmids_mentions = len({x['pmid'] for x in annotated_mentions})

RELS_PRINT = ['-'.join(rel) for rel in relations]
RELS_PRINT = ' or '.join(RELS_PRINT)


stat2 = f'- {len(annotated_mentions)} sentences with co-occurring ' + \
        f'{RELS_PRINT} mentions, corresponding to {pmids_mentions} abstracts.'
stats.append(stat2)

# Save files

files_dict = {f'abstracts_annotated_{suffix}.json': annotated_abstracts,
              f'abstracts_annotated_parsed_{suffix}.json': annotated_parsed,
              f'mentions_{suffix}.json': annotated_mentions}

for file_name, database in files_dict.items():
    with open(os.path.join(PATH, 'temp', file_name),
              mode='w',
              encoding='utf-8') as f:
        json.dump(database, f)
    sleep(1)


print('\n\n- Extracting relations from SPACY....... ', end='')
mentions_file = f'mentions_{suffix}.json'
if GPU:
    subprocess.call(['bash',
                     'scripts/launch_spacy_gpu.sh',
                     f"{os.path.join(PATH, 'temp', mentions_file)}"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
else:
    subprocess.call(['bash',
                     'scripts/launch_spacy_cpu.sh',
                     f"{os.path.join(PATH, 'temp', mentions_file)}"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)

print('DONE')


print('\n\n- Extracting relations from OpenIE...... ', end='')
mentions_file = f'mentions_{suffix}.json'
subprocess.call(['bash',
                 'scripts/launch_openie.sh',
                 f"{os.path.join(PATH, 'temp', mentions_file)}"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
print('DONE')


spacy_file = f'spacy_relations_{suffix}.json'
openie_file = f'openie_relations_{suffix}.json'

df_rel = get_relations(config, spacy_file, openie_file)
rel_sents = len(df_rel.sent_id.unique())
stat = f'- A total of {rel_sents} sentences with {RELS_PRINT} relations:'
stats.append(stat)

for rel in relations:
    df = df_rel[df_rel.relation_with == rel[1]]
    rel_pmids = len(df.pmid.unique())
    rel_sents = len(df.sent_id.unique())

    stat = f'\t- {rel_sents} sentences with virus-{rel[1]} relation, ' + \
           f'corresponding to {rel_pmids} abstracts.'
    stats.append(stat)

if len(relations) == 1:
    all_file = f'{relations[0][0]}_{relations[0][1]}_relations_{suffix}.json'
    df_rel.drop('relation_with', axis=1, inplace=True)
    df_rel.rename(columns={'other_entity': relations[0][1]}, inplace=True)
    df_rel.to_excel(os.path.join(PATH, all_file.replace('.json', '.xlsx')),
                    index=False)
    stat = f'- Extracted relations saved in {all_file}'
    stats.append(stat)
else:
    types = df_rel.relation_with.unique()
    if len(types) == 1:
        df_rel.drop('relation_with', axis=1, inplace=True)
        df_rel.rename(columns={'other_entity': types[0]}, inplace=True)
    all_file = f'all_relations_{suffix}.json'
    df_rel.to_excel(os.path.join(PATH, all_file.replace('.json', '.xlsx')),
                    index=False)
    stat = f'- All relations saved in {all_file}'
    stats.append(stat)

    if len(types) > 1:
        for t in types:
            df_temp = df_rel[df_rel['relation_with'] == t]
            # df_temp.drop('relation_with', axis=1, inplace=True)
            df_temp = df_temp.drop('relation_with', axis=1)
            df_temp.rename(columns={'other_entity': t}, inplace=True)
            rel_pmids_t = len(df_temp.pmid.unique())
            rel_sents_t = len(df_temp.sent_id.unique())
            temp_file = f'virus_{t}_relations_{suffix}.json'
            stat = f'- Virus-{t} relations saved in {temp_file}'
            stats.append(stat)
            
            df_temp.to_excel(os.path.join(PATH,
                                          temp_file.replace('.json', '.xlsx')),
                             index=False)

stats.append('- Temporary files saved in temp/')

print('\n\n')
print('*'*70)
for stat in stats:
    print(f'\n{stat}')
print('\n')
print('*'*70)

with open(os.path.join(PATH, f'stats_{suffix}.txt'),
          mode='w',
          encoding='utf-8') as f:
    for stat in stats:
        f.write(f'{stat}\n')

print('\n\n\nPROGRAM FINISHED\n')

sys.exit()
