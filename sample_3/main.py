import streamlit as st
from streamlit_ketcher import st_ketcher
import py3Dmol
from streamlit.components.v1 import html
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import pathlib

# ページ設定
st.set_page_config(layout="wide")

# カスタムCSS
st.markdown("""
    <style>
        .stApp > header {
            background-color: transparent;
        }
        .main > div {
            padding: 0rem !important;
        }
        .element-container {
            margin-bottom: 0.5rem;
        }
        [data-testid="stHorizontalBlock"] {
            align-items: stretch;
        }
        div[data-testid="stDataFrame"] td {
            cursor: pointer;
        }
        .block-container {
            padding-top: 1rem !important;
            padding-bottom: 0rem !important;
        }
    </style>
""", unsafe_allow_html=True)

# アミノメチルシクロヘキサンのSMILES
AMINOMETHYL_CYCLOHEXANE = "NCC1CCCCC1"

# シミュレーション設定のモーダルダイアログ
@st.dialog("Simulation Settings")
def show_simulation_settings():
    st.write("Configure advanced docking parameters")
    
    # タブで設定を整理
    tab1, tab2 = st.tabs(["Basic", "Advanced"])
    
    with tab1:
        st.subheader("Grid Settings")
        grid_size = st.slider("Grid Size (Å)", 16, 32, 24)
        center_x = st.number_input("Center X", value=0.0)
        center_y = st.number_input("Center Y", value=0.0)
        center_z = st.number_input("Center Z", value=0.0)
        
    with tab2:
        st.subheader("Advanced Parameters")
        exhaustiveness = st.slider("Exhaustiveness", 1, 32, 8)
        num_modes = st.slider("Number of Binding Modes", 1, 20, 9)
        energy_range = st.slider("Energy Range (kcal/mol)", 1, 10, 3)
        st.checkbox("Use GPU acceleration", value=True)
        st.checkbox("Consider water molecules", value=False)
        scoring_function = st.selectbox(
            "Scoring Function",
            ["Vina", "Vinardo", "AD4"]
        )
    
    if st.button("Apply Settings"):
        st.session_state.simulation_params = {
            "grid_size": grid_size,
            "center": (center_x, center_y, center_z),
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "energy_range": energy_range,
            "scoring": scoring_function
        }
        st.rerun()

def main():
    st.title("Reverse Docking Tool")

    if 'simulation_params' not in st.session_state:
        st.session_state.simulation_params = None

    # 上部グリッド
    top_left, top_right = st.columns(2, gap="medium")

    with top_left:
        # Ketcherエディタ
        molecule = st.text_input("Input SMILES", AMINOMETHYL_CYCLOHEXANE, label_visibility="collapsed")
        smile_code = st_ketcher(molecule, height=450)
        st.markdown(f"**SMILES:** `{smile_code}`")

    with top_right:
        st.subheader("Search Results")
        # サンプルデータ
        data = {
            'PDB ID': ['9BHQ', '1ABC', '2XYZ', '3DEF', '4GHI'],
            'Similarity': [95, 88, 82, 75, 70],
            'Resolution': [1.8, 2.1, 1.9, 2.3, 2.0],
            'Species': ['H. sapiens', 'H. sapiens', 'M. musculus', 'H. sapiens', 'E. coli'],
            'Method': ['X-ray', 'X-ray', 'Cryo-EM', 'X-ray', 'X-ray']
        }
        df = pd.DataFrame(data)
        
        st.data_editor(
            df,
            use_container_width=True,
            hide_index=True,
            column_config={
                "PDB ID": st.column_config.TextColumn(
                    "PDB ID",
                    help="Click to view structure",
                    width="small",
                ),
                "Similarity": st.column_config.ProgressColumn(
                    "Similarity (%)",
                    format="%d",
                    min_value=0,
                    max_value=100,
                ),
                "Resolution": st.column_config.NumberColumn(
                    "Resolution (Å)",
                    format="%.1f",
                ),
            },
            key='protein_table'
        )

    # 下部グリッド
    bottom_left, bottom_right = st.columns(2, gap="medium")

    with bottom_left:
        st.subheader("3D Ligand View")
        if smile_code:
            try:
                mol = Chem.MolFromSmiles(smile_code)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)
                pdb = Chem.MolToPDBBlock(mol)
                
                view = py3Dmol.view(width=400, height=400)
                view.addModel(pdb, "pdb")
                view.setStyle({'stick':{}})
                view.zoomTo()
                html(view._make_html(), height=400)
            except Exception as e:
                st.error(f"Error generating 3D structure: {e}")
        
        # シミュレーションコントロール
        st.subheader("Simulation Control")
        control_col1, control_col2 = st.columns(2)
        with control_col1:
            if st.button("Run Docking", use_container_width=True, type="primary"):
                st.info("Starting docking simulation...")
        with control_col2:
            if st.button("⚙️ Settings", use_container_width=True):
                show_simulation_settings()
        
        # 設定表示
        st.divider()
        if st.session_state.simulation_params:
            params = st.session_state.simulation_params
            st.markdown("#### Current Settings")
            st.markdown(f"""
            - **Grid size:** {params['grid_size']}Å
            - **Exhaustiveness:** {params['exhaustiveness']}
            - **Modes:** {params['num_modes']}
            - **Scoring:** {params['scoring']}
            """)
        else:
            st.markdown("#### Default Settings")
            st.markdown("""
            - **Grid size:** 24Å
            - **Exhaustiveness:** 8
            - **Modes:** 9
            - **Scoring:** Vina
            """)

    with bottom_right:
        # タンパク質ビューア
        st.subheader("Protein Structure")
        pdb_path = pathlib.Path('./data/9bhq.pdb')
        if pdb_path.exists():
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
                view = py3Dmol.view(width=600, height=400)
                view.addModel(pdb_content, "pdb")
                view.setStyle({'cartoon': {'color': 'spectrum'}})
                view.zoomTo()
                html(view._make_html(), height=400)
        else:
            st.info("Select a protein from the results table to view structure")

if __name__ == "__main__":
    main()