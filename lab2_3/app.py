import io
import numpy as np
import pandas as pd
import streamlit as st
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
import altair as alt

st.set_page_config(page_title="fasta sliding window frequencies", layout="wide")

def read_fasta_to_sequence(text: str) -> str:
    lines = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('>'):
            continue
        lines.append(line)
    seq = "".join(lines).upper().replace(" ", "")
    return seq

def sliding_window_frequencies(seq: str, window_size: int) -> pd.DataFrame:
    n = len(seq)
    if window_size <= 0 or window_size > n:
        return pd.DataFrame()

    alphabet = sorted(set(seq))
    total_windows = n - window_size + 1

    records = []
    for i in range(total_windows):
        window = seq[i:i+window_size]
        counts = {sym: 0 for sym in alphabet}
        for ch in window:
            if ch in counts:
                counts[ch] += 1
        freqs = {sym: counts[sym] / float(window_size) for sym in alphabet}
        freqs["window"] = i
        records.append(freqs)

    df = pd.DataFrame(records).set_index("window")
    df = df[alphabet]
    return df

def run_kmeans(df: pd.DataFrame, n_clusters: int, random_state: int = 42):
    """run k-means on frequency vectors and provide labels, silhouette score, and 2d pca projection"""
    if df.empty or df.shape[0] < n_clusters:
        return None, None, None

    X = df.values
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=random_state)
    labels = km.fit_predict(X)

    sil = None
    if X.shape[0] > n_clusters and len(np.unique(labels)) > 1:
        try:
            sil = float(silhouette_score(X, labels))
        except Exception:
            sil = None

    pca = PCA(n_components=2, random_state=random_state)
    X2 = pca.fit_transform(X)
    pca_df = pd.DataFrame({
        "pc1": X2[:, 0],
        "pc2": X2[:, 1],
        "window": df.index.values,
        "cluster": labels.astype(int),
    })

    return labels, sil, pca_df


st.title("sliding window frequencies from fasta")
st.write("upload a fasta file. i compute relative symbol frequencies per sliding window and plot them.")

with st.sidebar:
    st.header("parameters")
    window_size = st.number_input("window size", min_value=2, max_value=10000, value=30, step=1)
    use_ai = st.checkbox("enable ai analysis (k-means)", value=False)
    if use_ai:
        n_clusters = st.slider("number of clusters (k)", min_value=2, max_value=8, value=3, step=1)

uploaded = st.file_uploader("choose a fasta file", type=["fa", "fasta", "txt"])

if uploaded is not None:
    try:
        raw = uploaded.read().decode("utf-8", errors="ignore")
    except Exception:
        st.error("cannot decode file. please use utf-8.")
        st.stop()

    seq = read_fasta_to_sequence(raw)
    st.write(f"sequence length: **{len(seq)}**")
    if len(seq) == 0:
        st.error("file does not contain a valid sequence.")
        st.stop()

    if len(seq) < window_size:
        st.error(f"sequence has {len(seq)} characters. increase sequence length or reduce window size {window_size}.")
        st.stop()

    df = sliding_window_frequencies(seq, window_size)

    st.subheader("relative frequencies per window")
    st.caption("index is the start position of each window in the sequence.")
    st.dataframe(df, use_container_width=True)

    st.subheader("line chart")
    st.line_chart(df, x=None, y=list(df.columns))

    csv_buf = io.StringIO()
    df.to_csv(csv_buf)
    st.download_button("download csv", data=csv_buf.getvalue(), file_name="frequencies.csv", mime="text/csv")

    if use_ai:
        st.subheader("ai analysis: k-means clustering of windows")
        labels, sil, pca_df = run_kmeans(df, n_clusters)
        if labels is None:
            st.warning("cannot run k-means with current settings.")
        else:
            counts = pd.Series(labels).value_counts().sort_index()
            dist_df = pd.DataFrame({"cluster": counts.index.astype(int), "count": counts.values})
            st.write("cluster distribution:")
            st.bar_chart(dist_df.set_index("cluster"))

            if sil is not None:
                st.write(f"silhouette score: **{sil:.4f}**")

            lab_df = pd.DataFrame({"window": df.index, "cluster": labels})
            st.write("cluster labels per window:")
            st.dataframe(lab_df, use_container_width=True)

            st.write("2d pca projection colored by cluster:")
            chart = (
                alt.Chart(pca_df)
                .mark_circle(size=60)
                .encode(
                    x="pc1:Q",
                    y="pc2:Q",
                    color=alt.Color("cluster:N", legend=alt.Legend(title="cluster")),
                    tooltip=["window:Q", "cluster:N", "pc1:Q", "pc2:Q"],
                )
                .interactive()
            )
            st.altair_chart(chart, use_container_width=True)

            csv_labels = io.StringIO()
            lab_df.to_csv(csv_labels, index=False)
            st.download_button("download cluster labels", data=csv_labels.getvalue(),
                               file_name="cluster_labels.csv", mime="text/csv")
else:
    st.info("upload a fasta file to begin.")
