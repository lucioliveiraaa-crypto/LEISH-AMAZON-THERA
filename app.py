import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import pandas as pd
from st_mol import show_mol
import py3Dmol

# Configuração da Página
st.set_page_config(page_title="AmazonNatural Designer PRO", layout="wide", page_icon="🌿")

# CSS Personalizado para um visual moderno
st.markdown("""
    <style>
    .main { background-color: #fafafa; }
    .stMetric { background-color: #ffffff; padding: 15px; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
    .stButton>button { border-radius: 20px; height: 3em; background-color: #1b5e20; color: white; font-weight: bold; }
    .stSelectbox label, .stRadio label { font-weight: bold; color: #2e7d32; }
    h1, h2, h3 { color: #1b5e20; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
    </style>
    """, unsafe_allow_html=True)

# Função para Visualização 3D
def render_3d(smiles):
    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m)
    AllChem.MMFFOptimizeMolecule(m)
    m_block = Chem.MolToMolBlock(m)
    view = py3Dmol.view(width=400, height=400)
    view.addModel(m_block, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    return view

# Sidebar Elegante
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/2572/2572101.png", width=100)
    st.title("Configurações")
    st.divider()
    
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
    
    modo = st.radio("Fonte da Molécula:", ["🌿 Banco Amazônico", "⌨️ SMILES Manual"])
    if modo == "🌿 Banco Amazônico":
        selecao = st.selectbox("Selecione o Composto:", df["Composto"])
        smiles_base = df[df["Composto"] == selecao]["SMILES"].values[0]
    else:
        smiles_base = st.text_input("Cole o SMILES:")

# Cabeçalho
st.title("🌿 AmazonNatural Designer PRO")
st.markdown("---")

if smiles_base:
    mol = Chem.MolFromSmiles(smiles_base)
    if mol:
        # Colunas de Conteúdo
        tab1, tab2 = st.tabs(["📊 Análise e Propriedades", "🧬 Laboratório de Modificação"])

        with tab1:
            col_img, col_metrics = st.columns([1, 1])
            with col_img:
                st.subheader("Estrutura 2D")
                st.image(Draw.MolToImage(mol, size=(400, 400)), use_column_width=True)
            
            with col_metrics:
                st.subheader("Propriedades Físico-Químicas")
                c1, c2 = st.columns(2)
                c1.metric("Peso Molecular", f"{Descriptors.MolWt(mol):.2f}")
                c2.metric("LogP", f"{Descriptors.MolLogP(mol):.2f}")
                c1.metric("H-Donors", Descriptors.NumHDonors(mol))
                c2.metric("H-Acceptors", Descriptors.NumHAcceptors(mol))
                
                st.divider()
                st.subheader("Visualização 3D Interativa")
                show_mol(render_3d(smiles_base), height=400, width=400)

        with tab2:
            col_mod_ui, col_mod_view = st.columns([1, 1])
            
            with col_mod_ui:
                st.subheader("Ações de Bioisosterismo")
                regra = st.selectbox("Escolha a estratégia:", [
                    "Nenhuma", "Fluorar Carbono", "Metilar Hidroxila", "Trocar -OH por -F"
                ])
                
                new_mol = mol
                if regra == "Fluorar Carbono":
                    rxn = AllChem.ReactionFromSmarts('[C:1][H]>>[C:1]F')
                    p = rxn.RunReactants((mol,))
                    if p: new_mol = p[0][0]
                elif regra == "Metilar Hidroxila":
                    rxn = AllChem.ReactionFromSmarts('[O:1][H]>>[O:1]C')
                    p = rxn.RunReactants((mol,))
                    if p: new_mol = p[0][0]
                elif regra == "Trocar -OH por -F":
                    rxn = AllChem.ReactionFromSmarts('[O:1][H]>>[F:1]')
                    p = rxn.RunReactants((mol,))
                    if p: new_mol = p[0][0]
                
                st.divider()
                # Lipinski Check
                lipinski = Descriptors.MolWt(new_mol) <= 500 and Descriptors.MolLogP(new_mol) <= 5
                if lipinski:
                    st.success("✅ Atende à Regra de Lipinski")
                else:
                    st.warning("⚠️ Fora dos parâmetros de Lipinski")
                
                sdf = Chem.MolToMolBlock(new_mol)
                st.download_button("💾 Exportar SDF para Docking", sdf, file_name="lead_otimizado.sdf")

            with col_mod_view:
                st.subheader("Resultado da Otimização")
                st.image(Draw.MolToImage(new_mol, size=(400, 400)), use_column_width=True)
                st.write(f"**Novo LogP:** {Descriptors.MolLogP(new_mol):.2f}")

st.divider()
st.caption("Desenvolvido por Luciana | Doutorado em Parasitologia | UFPA/UEPA")
