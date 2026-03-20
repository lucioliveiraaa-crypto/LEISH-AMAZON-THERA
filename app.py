import streamlit as st
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
import pandas as pd
import pandas as pd

# METADADOS DE AUTORIA (INPI)
# Autora: Luciana - Doutorado em Parasitologia
# Título: AmazonNatural Designer v1.0

st.set_page_config(page_title="AmazonNatural Designer", layout="wide")

st.title("🌿 AmazonNatural Designer")
st.markdown("Plataforma de Otimização Molecular para Produtos Naturais da Amazônia")

# Banco de Dados de Exemplo (Cheiro-Verde e outros)
data = {
    "Planta": ["Coentro", "Cebolinha", "Copaíba", "Açaí"],
    "Composto": ["Linalol", "Limoneno", "Beta-cariofileno", "Quercetina"],
    "SMILES": ["CC(C)=CCCC(C)(O)C=C", "CC1=CCC(CC1)C(=C)C", "CC1=C[C@@H]2CC[C@H](C1)C(=C)[C@H]2CC1", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"]
}
df = pd.DataFrame(data)

# Sidebar
with st.sidebar:
    st.header("Entrada de Molécula")
    modo = st.selectbox("Escolha a origem:", ["Banco Amazônico", "SMILES Manual"])
    if modo == "Banco Amazônico":
        selecao = st.selectbox("Composto:", df["Composto"])
        smiles = df[df["Composto"] == selecao]["SMILES"].values[0]
    else:
        smiles = st.text_input("Insira o SMILES:")

# Processamento e Visualização
if smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        col1, col2 = st.columns(2)
        with col1:
            st.info("Original")
            st.image(Draw.MolToImage(mol))
            st.write(f"LogP: {round(Descriptors.MolLogP(mol), 2)}")
        with col2:
            st.success("Modificação")
            regra = st.radio("Alteração Bioisostérica:", ["Nenhuma", "Adicionar Flúor (Aromático)", "Metilar -OH"])
            
            new_mol = mol
            if regra == "Adicionar Flúor (Aromático)":
                rxn = AllChem.ReactionFromSmarts('[c:1][H]>>[c:1]F')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]
            elif regra == "Metilar -OH":
                rxn = AllChem.ReactionFromSmarts('[O:1][H]>>[O:1]C')
                p = rxn.RunReactants((mol,))
                if p: new_mol = p[0][0]
                
            st.image(Draw.MolToImage(new_mol))
            st.write(f"Novo LogP: {round(Descriptors.MolLogP(new_mol), 2)}")
            
            # Botão de Download para Docking
            sdf = Chem.MolToMolBlock(new_mol)
            st.download_button("Baixar .SDF para Docking", sdf, file_name="otimizada.sdf")
