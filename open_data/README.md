# Data availability
<table>
  <tr>
    <th colspan="3">Human postmortem – snRNAseq processing</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li>fastq files</li>
        <li>raw aggregated counts matrix</li>
        <li>QC’ed aggregated counts matrix</li>
        <li>patient metadata</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Raw sequencing reads</li>
        <li>Gene × cell count matrix</li>
        <li>Quality-filtered count matrix</li>
        <li>Sample (patient) metadata</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: NA</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_processing">analyses/snRNAseq_processing</a></li>
        <li>Available: Raw</li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>
      <ul>
        <li>RefCellTypeMarkers.adultBrain.rds</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Reference cell-type marker definitions for adult brain cell types</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: PsychENCODE</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_processing/get_marker_genes.ipynb">get_marker_genes.ipynb</a></li>
        <li>Available: Accessory</li>
      </ul>
    </td>
  </tr>

  <tr>
    <th colspan="3">Human postmortem – snRNAseq statistics</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li>all_paths.csv</li>
        <li>WikiPathways_2019_Human.npy</li>
        <li>DEGs (soon)</li>
        <li>pathway enrichments (soon)</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>List of all analysis file paths</li>
        <li>Pre-compiled pathway gene sets</li>
        <li>Differential expression results (pending)</li>
        <li>Enrichment summaries (pending)</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: WikiPathways_2019_Human</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_stats/compute_stats.ipynb">compute_stats.ipynb</a></li>
        <li>Available: Processed</li>
      </ul>
    </td>
  </tr>

  <tr>
    <th colspan="3">Human postmortem – gene-pathway partitioning</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li>all_scores_0825.csv</li>
        <li>leading_edge_0825Ex.csv</li>
        <li><code>*_loss.npy</code></li>
        <li><code>*_labs.npy</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Per-gene, per-pathway score matrix</li>
        <li>Leading-edge gene lists for each pathway</li>
        <li>NumPy arrays of cluster losses & labels</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/snRNAseq_stats/compute_stats.ipynb">compute_stats.ipynb</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_score_partitioning/projections.ipynb">projections.ipynb</a></li>
        <li>Available: Processed</li>
      </ul>
    </td>
  </tr>

  <tr>
    <th colspan="3">Induced-neurons – LCMS</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li>5041.SUB14737_LipidXData.csv</li>
        <li>1096.SUB12877_lipidXData.csv</li>
        <li>2685.SUB15127_LipidXData.csv</li>
        <li>7689.SUB15127_MetabolomicsData_Targeted.xlsx</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Quantified lipid features (three batches)</li>
        <li>Targeted metabolomics summary spreadsheet</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: Harvard LCMS Core</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/iN_LCMS/lipidomics/SUB14737_lipidomics_choline.ipynb">SUB14737_lipidomics_choline.ipynb</a></li>
        <li>Available: Processed</li>
      </ul>
    </td>
  </tr>
  <tr>
    <td colspan="3">
      <ul>
        <li>(placeholder for raw LCMS instrument files, e.g. vendor <code>.raw</code>)</li>
        <li>Placeholder for raw LCMS instrument output files</li>
        <li>Source: NA</li>
        <li>Used in: NA</li>
        <li>Available: Raw</li>
      </ul>
    </td>
  </tr>

  <tr>
    <th colspan="3">Induced-neurons – bulk RNA-seq</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li>fastq files</li>
        <li>counts.txt</li>
        <li>g2_degs.csv</li>
        <li>y622_degs.csv</li>
        <li>choline_degs.csv</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Raw reads for each condition</li>
        <li>Gene count matrix</li>
        <li>DEG tables for G2, Y622, and choline treatments</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/bulkRNAseq">bulkRNAseq</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/bulkRNAseq">bulkRNAseq</a></li>
        <li>Available: Processed</li>
      </ul>
    </td>
  </tr>

  <tr>
    <th colspan="3">Induced-neurons – O₂ consumption rates</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li>batch*.csv</li>
        <li>batch_*_df_quant.csv</li>
        <li><code>*.xlsx</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Processed oxygen consumption rates (per run)</li>
        <li>Quantified summary tables</li>
        <li>Raw Excel exports from Seahorse analyzer</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/iN_O2_consumption/seahorse_updpated.ipynb">seahorse_updpated.ipynb</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/iN_O2_consumption/seahorse_updpated.ipynb">seahorse_updpated.ipynb</a></li>
        <li>Available: Processed</li>
      </ul>
    </td>
  </tr>
</table>

