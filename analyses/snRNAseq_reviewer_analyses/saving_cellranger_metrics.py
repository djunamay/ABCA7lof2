from bs4 import BeautifulSoup
import ast
import glob
import pandas as pd
import re

meta = pd.read_csv('./raw_data/metadata/single_cell_individual_metadata.csv')

### Extracts the dictionaries from the HTML file from cellranger 

def extract_dictionaries(text):
    results = []
    stack = []
    # Keep track of all (start, end) pairs that correspond to a dictionary substring.
    for i, char in enumerate(text):
        if char == '{':
            stack.append(i)
        elif char == '}':
            if stack:
                start = stack.pop()
                candidate = text[start:i+1]
                try:
                    # Attempt to convert the substring to a Python object.
                    obj = ast.literal_eval(candidate)
                    if isinstance(obj, dict):
                        results.append(obj)
                except Exception:
                    # If evaluation fails, skip this substring.
                    pass
    return results

def extract_html_data(path):
    with open(path, 'r', encoding='utf-8') as f:
        html = f.read()

    # Parse the HTML
    soup = BeautifulSoup(html, 'html.parser')
    script_tag = soup.find("script", type="text/javascript", text=re.compile("const data ="))
    match = re.search(r"const data = (.*);", script_tag.string, re.DOTALL)
    json_text = match.group(1)

    dictionaries = extract_dictionaries(json_text)

    df = pd.DataFrame(dictionaries[-1]['rows']).T
    df.columns = df.iloc[0]
    df = df[1:].reset_index(drop=True)
    return df


summary_metrics_out = {}
for path in glob.glob('./raw_data/cellranger_counts_out/**/outs/metrics_summary.csv', recursive=True):
    data = pd.read_csv(path)
    summary_metrics_out[path.split('/')[3]] = data
    
  
other_metrics = {}
for path in glob.glob('./raw_data/cellranger_counts_out/**/outs/web_summary.html', recursive=True):
    data = extract_html_data(path)
    other_metrics[path.split('/')[3]] = data
    
all_data = pd.concat([pd.concat([other_metrics[x], summary_metrics_out[x]], axis=1) for x in other_metrics.keys()])
all_data.drop('Reference Path', axis=1, inplace=True)

all_data['study/sequencing batch'] = ['this study' if x in set(meta[meta['seq_batch']=='JBM']['sample_id']) else 'Mathys,2019' for x in all_data['Sample ID']]

all_data.to_csv('snrnaseq_cellranger_metrics.csv', index=False)
