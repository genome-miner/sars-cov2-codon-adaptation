# =======================================================================
# Project Title:
# AI-Assisted Computational Framework for mRNA Vaccine Sequence
# Analysis and Optimization
#
# Research Question:
# Can codon usage bias and adaptation of SARS-CoV-2 variant
# coding sequences to the human host be quantitatively characterized
# using compositional and codon-level features, and can these
# features predict the Codon Adaptation Index (CAI), a key
# determinant of translational efficiency relevant to mRNA
# vaccine antigen design?
# =======================================================================


# =======================================================================
# Import Required Libraries
# =======================================================================

import os
import json
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
from scipy.stats import kruskal
import matplotlib.pyplot as plt
from Bio.Data import CodonTable
from collections import Counter
from scipy.spatial.distance import jensenshannon
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error


# ======================================================================
# Base Directory
# ======================================================================
BASE_DIR = Path.cwd()
DATA_FILES = BASE_DIR / "Data"
FIGURES_DIR = BASE_DIR / "Figures"

# ======================================================================
# Reproducibility
# ======================================================================

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)


# ======================================================================
# Phase 1:
# Viral GenBank Dataset Processing
# ======================================================================

dictionary = {
    "COVID_19": [
        DATA_FILES / "COVID_19_Complete_Sequence.gb",
        DATA_FILES / "COVID_19_Extra1_OV093088.1.gb",
        DATA_FILES / "COVID_19_Extra2_OU161735.1.gb"
    ],
    "Omicron_strain": [
       DATA_FILES / "Omicron_Strain_Complete_Sequence.gb",
       DATA_FILES / "Omicron_Strain_Extra1_OR529199.1.gb",
       DATA_FILES / "Omicron_Strain_Extra2_ON115270.1.gb",
       DATA_FILES / "Omicron_Strain_Extra3_ON115272.1.gb"
    ],
    "Beta_strain": [
        DATA_FILES / "Beta_Strain_Complete_Sequence.gb",
        DATA_FILES / "Beta_Strain_Extra1_OR936719.gb"
    ],
    "Delta_strain": [
        DATA_FILES / "Delta_Strain_Complete_Sequence.gb",
        DATA_FILES / "Delta_Strain_Extra1_OR936720.1.gb"
    ]
}

data_record = []

allowed_bases = {"A", "T", "G", "C"}
allowed_bases_bytes = np.array([b"A", b"T", b"G", b"C"])

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
sense_codons = set(standard_table.forward_table.keys())
stop_codons = set(standard_table.stop_codons)

ordered_codons = np.array(sorted(sense_codons))
n_codons_total = len(ordered_codons)
stop_codons_array = np.array(sorted(stop_codons))
aa_lookup_array = np.array([standard_table.forward_table[c] for c in ordered_codons])

global_unique_sequences = set()


# ======================================================================
# Convert NumPy Objects to Native Python Objects
# ======================================================================

def convert_to_python(value):
    if isinstance(value, dict):
        return {str(convert_to_python(key)): convert_to_python(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [convert_to_python(item) for item in value]
    if isinstance(value, np.ndarray):
        return [convert_to_python(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def dictionary_to_json(value):
    return json.dumps(convert_to_python(value), separators=(",", ":"), ensure_ascii=False)


# ======================================================================
# Vectorized Codon Extraction
# ======================================================================

def seq_to_codons_np(seq):
    n = len(seq)
    usable_len = (n // 3) * 3
    if usable_len == 0:
        return np.array([], dtype="<U3")
    byte_arr = np.frombuffer(seq[:usable_len].encode("ascii"), dtype="S1")
    codon_bytes = byte_arr.reshape(-1, 3).view("S3").ravel()
    return codon_bytes.astype("<U3")


# ======================================================================
# CDS Validation
# ======================================================================

def validate_cds(cds_sequence):
    validation = {"Valid": False, "Reason": None, "Start_codon": False, "Stop_codon": False, "Internal_stop": False, "Protein_sequence": None}

    if len(cds_sequence) == 0:
        validation["Reason"] = "Empty sequence"
        return validation

    seq_arr = np.frombuffer(cds_sequence.encode("ascii"), dtype="S1")
    invalid_mask = ~np.isin(seq_arr, allowed_bases_bytes)

    if invalid_mask.any():
        invalid_bases = sorted(set(seq_arr[invalid_mask].astype("U1").tolist()))
        validation["Reason"] = f"Invalid bases: {invalid_bases[0]}"
        return validation

    if len(cds_sequence) % 3 != 0:
        validation["Reason"] = "Length not divisible by 3"
        return validation

    codons = seq_to_codons_np(cds_sequence)

    validation["Start_codon"] = bool(codons[0] == "ATG")
    validation["Stop_codon"] = bool(codons[-1] in stop_codons)

    if not validation["Start_codon"]:
        validation["Reason"] = "Missing start codon"
        return validation

    if not validation["Stop_codon"]:
        validation["Reason"] = "Missing terminal stop codon"
        return validation

    coding_codons = codons[:-1]

    is_stop = np.isin(coding_codons, stop_codons_array)
    validation["Internal_stop"] = bool(is_stop.any())

    if validation["Internal_stop"]:
        validation["Reason"] = "Internal stop codon detected"
        return validation

    is_sense = np.isin(coding_codons, ordered_codons)

    if not is_sense.all():
        bad_codon = coding_codons[~is_sense][0]
        validation["Reason"] = f"Invalid sense codon: {bad_codon}"
        return validation

    idx = np.searchsorted(ordered_codons, coding_codons)
    validation["Protein_sequence"] = "".join(aa_lookup_array[idx].tolist())

    validation["Valid"] = True
    validation["Reason"] = "Valid"

    return validation


# ======================================================================
# Process Each GenBank Dataset
# ======================================================================

for key, file_list in dictionary.items():

    unique_sequences = set()

    print("\n==============================")
    print(f"Processing dataset: {key}")
    print("==============================")

    for value in file_list:

        if not os.path.exists(value):
            raise FileNotFoundError(f"GenBank file not found: {value}")

        for record in SeqIO.parse(value, "genbank"):

            print(f"Record ID: {record.id}")

            for feature in record.features:

                if feature.type != "CDS":
                    continue

                try:
                    cds_sequence = str(feature.extract(record.seq)).upper()
                except Exception as error:
                    print("Skipping CDS because extraction failed:")
                    print(error)
                    continue

                validation = validate_cds(cds_sequence)

                if not validation["Valid"]:
                    print("Skipping CDS:")
                    print(validation["Reason"])
                    continue

                cds_length = len(cds_sequence)

                if cds_sequence in unique_sequences:
                    continue

                unique_sequences.add(cds_sequence)

                global_duplicate = cds_sequence in global_unique_sequences
                global_unique_sequences.add(cds_sequence)

                gene_name = str(feature.qualifiers.get("gene", ["Unknown"])[0])
                product_name = str(feature.qualifiers.get("product", ["Unknown"])[0])
                protein_id = str(feature.qualifiers.get("protein_id", ["Unknown"])[0])

                data_record.append({
                    "Data_name": key,
                    "ID": record.id,
                    "Gene": gene_name,
                    "Product": product_name,
                    "Protein_ID": protein_id,
                    "Sequence": cds_sequence,
                    "Protein_sequence": validation["Protein_sequence"],
                    "Length": cds_length,
                    "Protein_length": len(validation["Protein_sequence"]),
                    "Start_codon": bool(validation["Start_codon"]),
                    "Stop_codon": bool(validation["Stop_codon"]),
                    "Internal_stop": bool(validation["Internal_stop"]),
                    "Global_duplicate": bool(global_duplicate)
                    })

# ======================================================================
# Create CDS Metadata DataFrame
# ======================================================================

data = pd.DataFrame(data_record)

if data.empty:
    raise ValueError("No valid CDS sequences were extracted.")

global_duplicate_count = int(data["Global_duplicate"].sum())

print("\n==============================")
print("Duplicate Summary")
print("==============================")
print(f"Cross-dataset duplicate records: {global_duplicate_count}")

data = data[~data["Sequence"].duplicated(keep="first")].reset_index(drop=True)

print("\n==============================")
print("CDS Dataset Overview")
print("==============================")
print("Number of unique CDS sequences:", len(data))
print("\nDataset distribution:")
print(data["Data_name"].value_counts())

data.to_csv("Metadata.csv", index=False)


# ======================================================================
# Phase 2:
# Human Reference Dataset Processing
# ======================================================================

human_table = pd.read_csv(DATA_FILES / "Homo_sapiens_codon_usage.csv")

required_columns = {"Codon", "Amino Acid", "Fraction"}
missing_columns = required_columns - set(human_table.columns)

if missing_columns:
    raise ValueError(f"Missing columns in human codon usage table: {missing_columns}")

human_table["Codon"] = human_table["Codon"].astype(str).str.upper().str.strip()
human_table["Amino Acid"] = human_table["Amino Acid"].astype(str).str.strip()
human_table["Fraction"] = pd.to_numeric(human_table["Fraction"], errors="coerce")

human_table = human_table.dropna(subset=["Codon", "Amino Acid", "Fraction"])
human_table = human_table[human_table["Codon"].isin(sense_codons)].copy()

highest_fraction = human_table.groupby("Amino Acid")["Fraction"].transform("max")

if (highest_fraction <= 0).any():
    raise ValueError("One or more amino-acid groups have no positive reference codon usage.")

human_table["Weight"] = human_table["Fraction"] / highest_fraction

codon_weight = {str(codon): float(weight) for codon, weight in zip(human_table["Codon"], human_table["Weight"])}

missing_weights = sense_codons - set(codon_weight.keys())

if missing_weights:
    raise ValueError(f"Human codon table is missing weights for: {sorted(missing_weights)}")

CAI_PSEUDOCOUNT = 1e-6

weight_array = np.array([codon_weight[c] for c in ordered_codons], dtype=np.float64)
weight_array = np.maximum(weight_array, CAI_PSEUDOCOUNT)
log_weight_array = np.log(weight_array)


# ======================================================================
# Human CDS Codon Frequency Calculation
# ======================================================================

human_seq_record = SeqIO.parse(DATA_FILES / "cds_from_genomic.fna", "fasta")
human_counts_total = np.zeros(n_codons_total, dtype=np.int64)

for record in human_seq_record:

    human_cds = str(record.seq).upper()

    if len(human_cds) == 0 or len(human_cds) % 3 != 0:
        continue

    seq_arr = np.frombuffer(human_cds.encode("ascii"), dtype="S1")

    if not np.isin(seq_arr, allowed_bases_bytes).all():
        continue

    codons = seq_to_codons_np(human_cds)

    if codons.size == 0:
        continue

    is_sense = np.isin(codons, ordered_codons)
    sense_codons_this = codons[is_sense]

    if sense_codons_this.size == 0:
        continue

    idx = np.searchsorted(ordered_codons, sense_codons_this)
    human_counts_total += np.bincount(idx, minlength=n_codons_total)

total_human_codons = int(human_counts_total.sum())

if total_human_codons == 0:
    raise ValueError("No valid human sense codons were extracted.")

human_freq_array = human_counts_total.astype(np.float64) / total_human_codons
human_codon_freq = {str(codon): float(freq) for codon, freq in zip(ordered_codons.tolist(), human_freq_array.tolist())}

print("\n==============================")
print("Human Codon Reference")
print("==============================")
print("Number of human sense codons:", total_human_codons)
print("Number of codon weights:", len(codon_weight))
print("Number of codon frequencies:", len(human_codon_freq))


# ======================================================================
# Phase 3:
# Viral Sequence Feature Extraction
# ======================================================================

data_features = []
human_distribution = human_freq_array

for index, row in data.iterrows():

    sequence = str(row["Sequence"]).upper()
    sequence_length = len(sequence)

    # ---- Nucleotide composition ----

    seq_arr = np.frombuffer(sequence.encode("ascii"), dtype="S1")
    bases_unique, bases_counts = np.unique(seq_arr, return_counts=True)
    base_count_map = {str(base): int(count) for base, count in zip(bases_unique.astype("U1").tolist(), bases_counts.tolist())}

    G_content = base_count_map.get("G", 0)
    C_content = base_count_map.get("C", 0)
    A_content = base_count_map.get("A", 0)
    T_content = base_count_map.get("T", 0)

    GC_content = float(((G_content + C_content) / sequence_length) * 100)
    AT_content = float(((A_content + T_content) / sequence_length) * 100)

    GC_skew = float((G_content - C_content) / (G_content + C_content)) if (G_content + C_content) > 0 else 0.0
    AT_skew = float((A_content - T_content) / (A_content + T_content)) if (A_content + T_content) > 0 else 0.0

    # ---- Codon extraction ----

    codons = seq_to_codons_np(sequence)
    is_sense = np.isin(codons, ordered_codons)
    virus_sense_codons = codons[is_sense]

    if virus_sense_codons.size == 0:
        print(f"Skipping {row['ID']}: No sense codons available.")
        continue

    total_viral_sense_codons = int(virus_sense_codons.size)
    sense_idx = np.searchsorted(ordered_codons, virus_sense_codons)

    counts = np.bincount(sense_idx, minlength=n_codons_total)
    virus_distribution = counts.astype(np.float64) / total_viral_sense_codons
    virus_frequency = {str(codon): float(freq) for codon, freq in zip(ordered_codons.tolist(), virus_distribution.tolist())}

    # ---- Amino acid translation ----

    amino_acids_arr = aa_lookup_array[sense_idx]
    translated_sequence = "".join(amino_acids_arr.tolist())

    if translated_sequence != row["Protein_sequence"]:
        raise ValueError(f"Translation mismatch for record: {row['ID']}")

    # ---- CAI ----

    CAI = float(np.exp(np.mean(log_weight_array[sense_idx])))

    # ---- Virus vs Human codon usage ----

    Codon_usage_JSD = float(jensenshannon(virus_distribution, human_distribution, base=2))
    Codon_usage_similarity = float(1 - Codon_usage_JSD)

    # ---- Additional codon features ----

    Unique_codon_fraction = float(np.count_nonzero(counts) / n_codons_total)
    Most_frequent_codon_fraction = float(counts.max() / total_viral_sense_codons)

    nonzero_counts = counts[counts > 0]
    freqs = nonzero_counts.astype(np.float64) / total_viral_sense_codons
    Codon_usage_entropy = float(-np.sum(freqs * np.log2(freqs)))

    # ---- Amino acid composition ----

    aa_unique, aa_counts = np.unique(amino_acids_arr, return_counts=True)
    total_amino_acids = int(amino_acids_arr.size)
    amino_acid_composition = {str(aa): float((count / total_amino_acids) * 100) for aa, count in zip(aa_unique.tolist(), aa_counts.tolist())}
    amino_acids_list = "".join(str(aa) for aa in amino_acids_arr.tolist())

    # ---- Store values (dicts kept as native Python dicts, not JSON strings) ----

    data_features.append({
        "Data_name": str(row["Data_name"]),
        "ID": str(row["ID"]),
        "Gene": str(row["Gene"]),
        "Product": str(row["Product"]),
        "Length": int(row["Length"]),
        "Protein_length": int(row["Protein_length"]),
        "GC_content": GC_content,
        "AT_content": AT_content,
        "GC_skew": GC_skew,
        "AT_skew": AT_skew,
        "CAI": CAI,
        "Codon_usage_JSD": Codon_usage_JSD,
        "Codon_usage_similarity": Codon_usage_similarity,
        "Unique_codon_fraction": Unique_codon_fraction,
        "Most_frequent_codon_fraction": Most_frequent_codon_fraction,
        "Codon_usage_entropy": Codon_usage_entropy,
        "Virus_codon_freq_dict": virus_frequency,             
        "Amino_acid_composition_dict": amino_acid_composition, 
        "Amino_acids": amino_acids_list                        
    })


# ======================================================================
# Phase 4:
# Final Feature Matrix (single file, codon/AA dicts expanded
# into real numeric columns instead of JSON-string columns)
# ======================================================================


features_df = pd.DataFrame(data_features)

if features_df.empty:
    raise ValueError("Feature extraction produced an empty dataset.")

# ---- Expand Virus_codon_freq dict into 61 numeric columns ----

codon_freq_expanded = pd.DataFrame(features_df["Virus_codon_freq_dict"].tolist())
codon_freq_expanded = codon_freq_expanded.reindex(columns=ordered_codons.tolist(), fill_value=0.0)
codon_freq_expanded = codon_freq_expanded.add_prefix("Freq_")

# ---- Expand Amino_acid_composition dict into numeric columns ----

all_amino_acids = sorted(set(aa_lookup_array.tolist()))
aa_comp_expanded = pd.DataFrame(features_df["Amino_acid_composition_dict"].tolist())
aa_comp_expanded = aa_comp_expanded.reindex(columns=all_amino_acids, fill_value=0.0)
aa_comp_expanded = aa_comp_expanded.fillna(0.0)
aa_comp_expanded = aa_comp_expanded.add_prefix("AAcomp_")

# ---- Assemble final single feature matrix ----

scalar_columns = [
    "Data_name",
    "ID", "Gene",
    "Product",
    "Length",
    "Protein_length",
    "GC_content",
    "AT_content",
    "GC_skew",
    "AT_skew",
    "CAI",
    "Codon_usage_JSD",
    "Codon_usage_similarity",
    "Unique_codon_fraction",
    "Most_frequent_codon_fraction",
    "Codon_usage_entropy"
    ]

final_feature_matrix = pd.concat(
    [
        features_df[scalar_columns].reset_index(drop=True),
        codon_freq_expanded.reset_index(drop=True),
        aa_comp_expanded.reset_index(drop=True),
        features_df[["Amino_acids"]].reset_index(drop=True),
    ],
    axis=1
)

print("\n==============================")
print("Final Feature Matrix Summary")
print("==============================")

summary_columns = [
    "Length",
    "Protein_length",
    "GC_content",
    "AT_content",
    "GC_skew",
    "AT_skew",
    "CAI",
    "Codon_usage_JSD",
    "Codon_usage_similarity",
    "Unique_codon_fraction",
    "Most_frequent_codon_fraction",
    "Codon_usage_entropy"
    ]

print(final_feature_matrix[summary_columns].describe())

final_feature_matrix.to_csv("Feature_Matrix_Final.csv", index=False)

print("\nFinal feature matrix shape:", final_feature_matrix.shape)
print("\nFinal feature matrix preview:")
print(final_feature_matrix.head())

print("\n[SUCCESS] Single final feature matrix saved as 'Feature_Matrix_Final.csv'.")


# ======================================================================
# Phase 6: Exploratory Data Analysis (EDA)
# ======================================================================

os.makedirs("Figures", exist_ok=True)

freq_columns = [f"Freq_{c}" for c in ordered_codons]
aacomp_columns = [f"AAcomp_{aa}" for aa in sorted(set(aa_lookup_array.tolist()))]

# ----------------------------------------------------------------------
# CAI Boxplot
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
sns.boxplot(data=final_feature_matrix, x="Data_name", y="CAI")
plt.title("CAI Among Viral Variants")
plt.xlabel("Variant")
plt.ylabel("Codon Adaptation Index (CAI)")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "CAI_Boxplot.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# GC Content Boxplot
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
sns.boxplot(data=final_feature_matrix, x="Data_name", y="GC_content", color="seagreen")
plt.title("GC Content Among Viral Variants")
plt.xlabel("Variant")
plt.ylabel("GC Content")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "GC_Content_Boxplot.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# AT Content Boxplot
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
sns.boxplot(data=final_feature_matrix, x="Data_name", y="AT_content", color="orange")
plt.title("AT Content Among Viral Variants")
plt.xlabel("Variant")
plt.ylabel("AT Content")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "AT_Content_Boxplot.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# Codon Bias Heatmap (virus vs human, per codon, per variant)
# ----------------------------------------------------------------------
codon_bias_df = final_feature_matrix[["Data_name"] + freq_columns].copy()

for codon, freq_col in zip(ordered_codons, freq_columns):
    codon_bias_df[freq_col] = codon_bias_df[freq_col] - human_codon_freq[codon]

codon_bias_avg = codon_bias_df.groupby("Data_name")[freq_columns].mean()
codon_bias_avg.columns = ordered_codons  # relabel back to plain codon names for readability

plt.figure(figsize=(14, 6))
sns.heatmap(codon_bias_avg, center=0, cmap="coolwarm")
plt.title("Mean Codon Usage Bias vs Human Reference (per Variant)")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Codon_Bias_Heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# Feature Correlation Heatmap (scalar summary features)
# ----------------------------------------------------------------------
numeric_features = final_feature_matrix[
    ["GC_content", "AT_content", "Length", "GC_skew", "AT_skew", "CAI",
     "Codon_usage_similarity", "Codon_usage_entropy", "Unique_codon_fraction"]
]

correlation = numeric_features.corr()
print(correlation)

plt.figure(figsize=(8, 6))
sns.heatmap(correlation, annot=True, linewidth=0.05, cmap="coolwarm")
plt.title("Feature Correlation Heatmap")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Feature_Correlation_Heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# Scale scalar features first (needed for BOTH clustering and PCA)
# ----------------------------------------------------------------------
scaler = StandardScaler()
scaled_features = scaler.fit_transform(numeric_features)

# ----------------------------------------------------------------------
# Hierarchical Clustering (scalar summary features, SCALED)
# ----------------------------------------------------------------------
hierarchy = linkage(scaled_features, method="ward")

plt.figure(figsize=(10, 6))
dendrogram(hierarchy, labels=final_feature_matrix["Data_name"].tolist())
plt.title("Hierarchical Clustering of Viral Sequence Features")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Hierarchical_Clustering.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# PCA plot (scalar summary features)
# ----------------------------------------------------------------------
pca = PCA(n_components=2)
principal_components = pca.fit_transform(scaled_features)

pca_df = pd.DataFrame(principal_components, columns=["PC1", "PC2"])
pca_df["Data_name"] = final_feature_matrix["Data_name"]

print("Explained variance (scalar features):", pca.explained_variance_ratio_)

plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Data_name", s=100)
plt.title("PCA Plot of Viral Sequence Features (Scalar Summary)")
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "PCA_Plot.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# PCA and clustering on FULL codon usage space (61-dim)
# ----------------------------------------------------------------------
codon_freq_matrix = final_feature_matrix[freq_columns]

scaled_codon_freq = StandardScaler().fit_transform(codon_freq_matrix)

pca_codon = PCA(n_components=2)
codon_pcs = pca_codon.fit_transform(scaled_codon_freq)

codon_pca_df = pd.DataFrame(codon_pcs, columns=["PC1", "PC2"])
codon_pca_df["Data_name"] = final_feature_matrix["Data_name"]

print("Explained variance (codon usage space):", pca_codon.explained_variance_ratio_)

plt.figure(figsize=(8, 6))
sns.scatterplot(data=codon_pca_df, x="PC1", y="PC2", hue="Data_name", s=100)
plt.title("PCA Plot of Full Codon Usage Frequencies (61 Codons)")
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "PCA_Plot_Codon_Usage.png", dpi=300, bbox_inches="tight")
plt.close()

codon_hierarchy = linkage(scaled_codon_freq, method="ward")

plt.figure(figsize=(10, 6))
dendrogram(codon_hierarchy, labels=final_feature_matrix["Data_name"].tolist())
plt.title("Hierarchical Clustering Based on Full Codon Usage Profile")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Hierarchical_Clustering_Codon_Usage.png", dpi=300, bbox_inches="tight")
plt.close()

# ======================================================================
# Amino Acid Composition Visualization
# ======================================================================
protein_df = final_feature_matrix.melt(
    id_vars=["Data_name"],
    value_vars=aacomp_columns,
    var_name="Amino_Acid",
    value_name="Percentage"
)
protein_df["Amino_Acid"] = protein_df["Amino_Acid"].str.replace("AAcomp_", "", regex=False)

plt.figure(figsize=(12, 6))
sns.barplot(data=protein_df, x="Amino_Acid", y="Percentage", hue="Data_name")
plt.title("Amino Acid Composition Among Viral Variants")
plt.xlabel("Amino Acid")
plt.ylabel("Percentage")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Amino_Acid_Composition.png", dpi=300, bbox_inches="tight")
plt.close()

print("[SUCCESS] All EDA figures saved to /Figures")


# ======================================================================
# Phase 6b: Distribution Check and Statistical Significance
# ======================================================================

# ----------------------------------------------------------------------
# CAI Distribution Histogram (per variant, overlaid)
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
for variant in final_feature_matrix["Data_name"].unique():
    subset = final_feature_matrix.loc[final_feature_matrix["Data_name"] == variant, "CAI"]
    sns.histplot(subset, label=variant, kde=True, stat="density", alpha=0.4)

plt.title("Distribution of CAI Across Viral Variants")
plt.xlabel("Codon Adaptation Index (CAI)")
plt.ylabel("Density")
plt.legend(title="Variant")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "CAI_Distribution_Histogram.png", dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------------
# Kruskal-Wallis Test: is CAI significantly different across variants?
# ----------------------------------------------------------------------
cai_groups = [
    group["CAI"].values
    for _, group in final_feature_matrix.groupby("Data_name")
]

kw_stat, kw_pvalue = kruskal(*cai_groups)

print("\n==============================")
print("Kruskal-Wallis Test: CAI Across Variants")
print("==============================")
print(f"H-statistic: {kw_stat:.4f}")
print(f"p-value: {kw_pvalue:.6f}")

if kw_pvalue < 0.05:
    print("Result: Statistically significant difference in CAI across variants (p < 0.05).")
else:
    print("Result: No statistically significant difference in CAI across variants (p >= 0.05).")

print("\n[SUCCESS] Distribution and statistical test completed.")


# ======================================================================
# Phase 7: Machine Learning
# ======================================================================
scalar_input_cols = ["GC_content", "GC_skew", "AT_skew", "Length"]
X = final_feature_matrix[scalar_input_cols + aacomp_columns]
Y = final_feature_matrix["CAI"]
preprocessor = ColumnTransformer(transformers=[
    ("scalar", StandardScaler(), scalar_input_cols),
    ("aacomp", Pipeline([("scale", StandardScaler()), ("pca", PCA(n_components=3, random_state=RANDOM_SEED))]), aacomp_columns)
])

model = Pipeline([("preprocess", preprocessor), ("ridge", Ridge(alpha=0.1))])
loo = LeaveOneOut()
y_predicted = cross_val_predict(model, X, Y, cv=loo)

r2 = r2_score(Y, y_predicted)
mae = mean_absolute_error(Y, y_predicted)
rmse = np.sqrt(mean_squared_error(Y, y_predicted))

print(f"LOOCV R²: {r2:.4f}")
print(f"LOOCV MAE: {mae:.4f}")
print(f"LOOCV RMSE: {rmse:.4f}")

# Add Data_name (variant) as a feature and compare
X_with_variant = final_feature_matrix[scalar_input_cols + aacomp_columns + ["Data_name"]]

preprocessor_with_variant = ColumnTransformer(transformers=[
    ("scalar", StandardScaler(), scalar_input_cols),
    ("aacomp", Pipeline([("scale", StandardScaler()), ("pca", PCA(n_components=3, random_state=RANDOM_SEED))]), aacomp_columns),
    ("variant", OneHotEncoder(drop="first"), ["Data_name"])
])

model_with_variant = Pipeline([("preprocess", preprocessor_with_variant), ("ridge", Ridge(alpha=0.1))])
y_with_variant = cross_val_predict(model_with_variant, X_with_variant, Y, cv=loo)

r2_variant = r2_score(Y, y_with_variant)
mae_variant = mean_absolute_error(Y, y_with_variant)
rmse_variant = np.sqrt(mean_squared_error(Y, y_with_variant))

print(f"LOOCV R² (with variant): {r2_variant:.4f}")
print(f"LOOCV MAE (with variant): {mae_variant:.4f}")
print(f"LOOCV RMSE (with variant): {rmse_variant:.4f}")


# ======================================================================
# Step 8: Final Model Fit and Feature Coefficients
# ======================================================================

final_model = Pipeline([
    ("preprocess", preprocessor),
    ("ridge", Ridge(alpha=0.1))
])

final_model.fit(X, Y)

feature_names_out = final_model.named_steps["preprocess"].get_feature_names_out()
coefficients = final_model.named_steps["ridge"].coef_

coef_df = pd.DataFrame({
    "Feature": feature_names_out,
    "Coefficient": coefficients
}).sort_values(by="Coefficient", key=abs, ascending=False)

print("\n==============================")
print("Ridge Regression Feature Coefficients (sorted by magnitude)")
print("==============================")
print(coef_df.to_string(index=False))

coef_df.to_csv("Ridge_Feature_Coefficients.csv", index=False)
print("\n[SUCCESS] Feature coefficients saved as 'Ridge_Feature_Coefficients.csv'.")


# ======================================================================
# Step 9: Residual Diagnostics
# ======================================================================

y_pred_final = cross_val_predict(final_model, X, Y, cv=loo)
residuals = Y - y_pred_final

# Predicted vs Actual CAI value
plt.figure(figsize=(6,7))
plt.scatter(Y, y_pred_final, alpha=0.5)
plt.plot([Y.min(), Y.max()], [Y.min(), Y.max()], "r--", label="Perfect Prediction")
plt.xlabel("Actual CAI")
plt.ylabel("Predicted CAI")
plt.legend()
plt.savefig(FIGURES_DIR / "Predicted_vs_Actual_CAI.png", dpi=300, bbox_inches="tight")
plt.show()

# Residuals vs Predicted CAI
plt.figure(figsize=(6,7))
plt.scatter(y_pred_final, residuals, alpha=0.5)
plt.axhline(0, color="red", linestyle="--")
plt.xlabel("Predicted CAI")
plt.ylabel("Residual (Actual - Predicted)")
plt.title("Residuals vs Predicted CAI")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Residuals_vs_Predicted.png", dpi=300, bbox_inches="tight")
plt.show()

print("\n[SUCCESS] Residual diagnostic plots saved to /Figures")