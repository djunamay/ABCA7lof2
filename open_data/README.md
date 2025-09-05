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
        <li><code>*.fastq</code></li>
        <li><code>raw_counts.mtx</code></li>
        <li><code>qc_counts.mtx</code></li>
        <li><code>colData, rowData</code></li>
        <li><code>input_stats.rds [ ]</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Raw sequencing reads</li>
        <li>raw aggregated per-cell counts matrix</li>
        <li>QC’ed aggregated per-cell counts matrix</li>
        <li>per-cell and per-gene metadata</li>
        <li>Inputs for DEG analysis</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: NA</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_processing">analyses/snRNAseq_processing</a></li>
        <li>Available: <a href="https://www.synapse.org/#!Synapse:syn53461705">raw and unredacted data</a></li>
        <li>Available: <a href="https://singlecell.broadinstitute.org/single_cell/study/SCP3182/a-single-cell-atlas-of-abca7-loss-of-function-in-human-brain#study-download">redacted data</a></li>
        <li>Available: <a href="https://osf.io/mnysb/">input stats</a></li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>
      <ul>
        <li><code>PanglaoDB_Franzen2019.csv</code></li>
        <li><code>RefCellTypeMarkers.adultBrain.rds</code></li>
        <li><code>RefCellTypeMarkers.adultBrain_Wang2018.csv</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Reference cell-type marker definitions for adult brain cell types and excitatory layer annotations</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: PsychENCODE</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_processing">analyses/snRNAseq_processing</a></li>
        <li>Available: <a href="https://osf.io/5k8v2/files/osfstorage">celltype markers</a></li>
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
        <li><code>supplementary_table_6.csv</code></li>
        <li><code>WikiPathways_2019_Human.npy</code></li>
        <li><code>GO_Biological_Process_2023.npy</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>differentially expressed genes per celltype</li>
        <li>pathway databases</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: WikiPathways and Gene Ontology databases, <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_processing">analyses/snRNAseq_processing</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_stats">analyses/snRNAseq_stats</a></li>
        <li>Available: <a href="https://osf.io/v6y3d/files/osfstorage">pathways</a></li>
        <li>Available: <a href="https://osf.io/5cnfy/wiki/home/</li">DEGs</a>
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
        <li><code>umap_data_figure_1.npz</code></li>
        <li><code>kl_data_figure_2.npz</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>data related to the analysis in Figure 1</li>
        <li>data related to the analysis in Figure 2</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/snRNAseq_stats">analyses/snRNAseq_stats</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/snRNAseq_score_partitioning">analyses/snRNAseq_score_partitioning</a></li>
        <li>Available: <a href="https://osf.io/mnysb/">processed data</a></li>
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
        <li><code>*_LipidXData.csv</code></li>
        <li><code>MetabolomicsData_Targeted.xlsx</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>LCMS lipidomic quantifications for WT, p.Glu50fs*3, p.Tyr622*, p.Tyr622* +/- CDP-choline</li>
        <li>LCMS metabolomic data for WT, p.Tyr622*, p.Tyr622* +/- CDP-choline</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: Harvard LCMS Core</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/iN_LCMS">analyses/iN_LCMS</a></li>
        <li>Available: <a href="https://osf.io/mnysb/">processed data</a></li>
      </ul>
    </td>
  </tr>
<tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li><code>*.raw</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>raw LCMS lipidomic data for WT, p.Glu50fs*3, p.Tyr622*, p.Tyr622* +/- CDP-choline</li>
        <li>raw LCMS metabolomic data for WT, p.Tyr622*, p.Tyr622* +/- CDP-choline</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: Harvard LCMS Core</li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/iN_LCMS">analyses/iN_LCMS</a></li>
        <li>Available: <a href="https://osf.io/u68k3/">raw data</a></li>
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
        <li><code>*.fastq</code></li>
        <li><code>counts.txt</code></li>
        <li><code>supplementary_table_9.csv</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Raw reads for each condition</li>
        <li>Gene count matrix</li>
        <li>DEG tables for WT, p.Glu50fs*3, p.Tyr622*, p.Tyr622* +/- CDP-choline</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/bulkRNAseq">analyses/bulkRNAseq</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/tree/main/analyses/bulkRNAseq">analyses/bulkRNAseq</a></li>
        <li>Available: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299277">raw data</a></li>
        <li>Available: <a href="https://osf.io/5cnfy/wiki/home/">DEGs</a></li>
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
        <li><code>*.xlsx</code></li>
        <li><code>seahorse_quantifications.npz</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Raw Excel exports from Seahorse analyzer</li>
        <li>Processed and quantified data</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/iN_O2_consumption">analyses/iN_O2_consumption</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/iN_O2_consumption">analyses/iN_O2_consumption</a></li>
        <li>Available: <a href="https://osf.io/kptbw/">raw seahorse data</a></li>
        <li>Available: <a href="https://osf.io/mnysb/">processed seahorse data</a></li>
      </ul>
    </td>
  </tr>

  <tr>
    <th colspan="3">Induced-neurons / cortical organoids - Imaging Data</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li><code>*.czi</code></li>
        <li><code>image_quantifications.npz</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Raw images</li>
        <li>Image quantification</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/imaging">analyses/imaging</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/imaging">analyses/imaging</a></li>
        <li>Available: <a href="https://osf.io/htb23/">Raw Images</a></li>
        <li>Available: <a href="https://osf.io/mnysb/">Image quantifications</a></li>
      </ul>
    </td>
  </tr>

  <tr>
    <th colspan="3">Induced-neurons / cortical organoids – Amyloid and Ephys Data</th>
  </tr>
  <tr>
    <th>File name</th>
    <th>Description</th>
    <th>Details</th>
  </tr>
  <tr>
    <td>
      <ul>
        <li><code>*.xlsx</code></li>
        <li><code>*.csv</code></li>
        <li><code>*.xlsx</code></li>
      </ul>
    </td>
    <td>
      <ul>
        <li>amyloid ELISA data with analysis</li>
        <li>iN ephys data</li>
        <li>Raw ephys data</li>
      </ul>
    </td>
    <td>
      <ul>
        <li>Source: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/amyloid_ephys">analyses/amyloid_ephys</a></li>
        <li>Used in: <a href="https://github.com/djunamay/ABCA7lof2/blob/main/analyses/amyloid_ephys">analyses/amyloid_ephys</a></li>
        <li>Available: <a href="https://osf.io/jy3fe/">iN ephys data</a></li>
        <li>Available: <a href="https://osf.io/em92c/">cortical organoid ephys data</a></li>
        <li>Available: <a href="https://osf.io/szuvt/">ELISA data</a></li>
      </ul>
    </td>
  </tr>
</table>

