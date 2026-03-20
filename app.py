import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import pandas as pd

# Configuração e Estilo
st.set_page_config(page_title="AmazonNatural Designer", layout="wide")
st.markdown("""
    <style>
    .main { background-color: #f0f7f0; }
    .stMetric { background-color: #ffffff; padding: 10px; border-radius: 10px; border: 1px solid #e0e0e0; }
    h1 { color: #2e7d32; }
    </style>
    """, unsafe_allow_html=True)

st.title("🌿 AmazonNatural Designer v2.0")
st.caption("Ferramenta de Otimização Molecular - Doutorado em Parasitologia")

# Banco de Dados
data = {
    "Planta": ["Coentro", "Cebolinha", "Copaíba", "Açaí", "Andiroba"],
    "Composto": ["Linalol", "Limoneno", "Beta-cariofileno", "Quercetina", "Gedunina"],
    "SMILES": [
        "CC(C)=CCCC(C)(O)C=C", "CC1=CCC(CC1)C(=C)C", 
        "CC1=C[C@@H]2CC[C@H](C1)C(=C)[C@H]2CC1", 
        "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "CC12CCC(=O)C=C1C(C3C(C2)C45C(O3)C(C(C4(C)C)OC(=O)C)C6=COC=C6)C"
    ]
}
df = pd.DataFrame(data)

# Sidebar
with st.sidebar:
    st.header("⚙️ Configurações")
    modo = st.radio("Fonte:", ["Banco Amazônico", "SMILES Manual"])
    if modo == "Banco Amazônico":
        selecao = st.selectbox("Composto:", df["Composto"])
        smiles_base = df[df["Composto"] == selecao]["SMILES"].values[0]
    else:
        smiles_base = st.text_input("SMILES:")

if smiles_base:
    mol = Chem.MolFromSmiles(smiles_base)
    if mol:
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Análise Original")
            st.image(Draw.MolToImage(mol, size=(400, 400)))
            
            m1, m2 = st.columns(2)
            m1.metric("Peso Mol.", f"{Descriptors.MolWt(mol):.1f}")
            m2.metric("LogP", f"{Descriptors.MolLogP(mol):.1f}")

        with col2:
            st.subheader("Otimização Lead")
            regra = st.selectbox("Estratégia:", ["Nenhuma", "Fluorar (Alifático)", "Metilar Hidroxila", "Trocar OH por F"])
            
            new_mol = mol
            if regra == "Fluorar (Alifático)":
                rxn = AllChem.ReactionFromSmarts('[C:1][H]>>[C:1]F')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]
            elif regra == "Metilar Hidroxila":
                rxn = AllChem.ReactionFromSmarts('[O:1][H]>>[O:1]C')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]
            elif regra == "Trocar OH por F":
                rxn = AllChem.ReactionFromSmarts('[O:1][H]>>[F:1]')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]

            st.image(Draw.MolToImage(new_mol, size=(400, 400)))
            st.metric("Novo LogP", f"{Descriptors.MolLogP(new_mol):.1f}")
            
            sdf = Chem.MolToMolBlock(new_mol)
            st.download_button("💾 Baixar SDF para Docking", sdf, file_name="lead.sdf")

st.divider()
st.info("💡 Dica: Para o Linalol, use 'Fluorar (Alifático)' para ver a modificação na cadeia.")
