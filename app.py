import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import pandas as pd

# Configuração da Página
st.set_page_config(page_title="AmazonNatural Designer PRO", layout="wide")

st.markdown("""
    <style>
    .main { background-color: #f8fbf8; }
    .stButton>button { width: 100%; background-color: #2e7d32; color: white; border-radius: 10px; }
    </style>
    """, unsafe_allow_html=True)

st.title("🌿 AmazonNatural Designer | v2.0")
st.subheader("Otimização de Leads da Biodiversidade Amazônica")

# Banco de Dados Expandido
data = {
    "Planta": ["Coentro (Cheiro-Verde)", "Cebolinha", "Copaíba", "Açaí", "Andiroba"],
    "Composto": ["Linalol", "Limoneno", "Beta-cariofileno", "Quercetina", "Gedunina"],
    "SMILES": [
        "CC(C)=CCCC(C)(O)C=C", 
        "CC1=CCC(CC1)C(=C)C", 
        "CC1=C[C@@H]2CC[C@H](C1)C(=C)[C@H]2CC1", 
        "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "CC12CCC(=O)C=C1C(C3C(C2)C45C(O3)C(C(C4(C)C)OC(=O)C)C6=COC=C6)C"
    ]
}
df = pd.DataFrame(data)

# Interface Lateral
with st.sidebar:
    st.header("🧬 Seleção de Scaffold")
    modo = st.radio("Origem:", ["Banco Amazônico", "Custom SMILES"])
    if modo == "Banco Amazônico":
        selecao = st.selectbox("Escolha o Composto:", df["Composto"])
        smiles = df[df["Composto"] == selecao]["SMILES"].values[0]
    else:
        smiles = st.text_input("Cole o SMILES:")

# Lógica Principal
if smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.info("📌 Molécula Original")
            st.image(Draw.MolToImage(mol, size=(400, 400)))
            
            # Tabela de Propriedades
            props = {
                "Propriedade": ["Peso Molecular", "LogP", "H-Donors", "H-Acceptors"],
                "Valor": [
                    f"{Descriptors.MolWt(mol):.2f}",
                    f"{Descriptors.MolLogP(mol):.2f}",
                    Descriptors.NumHDonors(mol),
                    Descriptors.NumHAcceptors(mol)
                ]
            }
            st.table(pd.DataFrame(props))

        with col2:
            st.success("🛠️ Laboratório de Modificação")
            regra = st.selectbox("Estratégia Bioisostérica:", [
                "Nenhuma", 
                "Fluorar (Carbono Alifático - Ideal para Linalol)",
                "Adicionar Metila (-CH3)",
                "Trocar -OH por -F",
                "Metilar Hidroxila (-OCH3)"
            ])
            
            new_mol = mol
            # Novas regras mais inteligentes:
            if regra == "Fluorar (Carbono Alifático - Ideal para Linalol)":
                rxn = AllChem.ReactionFromSmarts('[C:1][H]>>[C:1]F')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]
            elif regra == "Adicionar Metila (-CH3)":
                rxn = AllChem.ReactionFromSmarts('[C,N,O:1][H]>>[C,N,O:1]C')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]
            elif regra == "Trocar -OH por -F":
                rxn = AllChem.ReactionFromSmarts('[O:1][H]>>[F:1]')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]
            elif regra == "Metilar Hidroxila (-OCH3)":
                rxn = AllChem.ReactionFromSmarts('[O:1][H]>>[O:1]C')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]

            st.image(Draw.MolToImage(new_mol, size=(400, 400)))
            
            # Badge de Lipinski
            if Descriptors.MolWt(new_mol) <= 500 and Descriptors.MolLogP(new_mol) <= 5:
                st.write("✅ **Regra de Lipinski: Aprovada**")
            else:
                st.write("⚠️ **Regra de Lipinski: Alerta**")

            # Download
            sdf = Chem.MolToMolBlock(new_mol)
            st.download_button("💾 Baixar SDF para Docking", sdf, file_name="amazon_lead.sdf")

st.divider()
st.caption("AmazonNatural Designer - Desenvolvido para Pesquisa Acadêmica em Parasitologia.")
