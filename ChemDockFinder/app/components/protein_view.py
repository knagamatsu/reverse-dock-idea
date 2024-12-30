import streamlit as st
import py3Dmol
from streamlit.components.v1 import html
import pathlib
from typing import Optional, Dict
import logging

logger = logging.getLogger(__name__)

def load_pdb_content(pdb_file: pathlib.Path) -> Optional[str]:
    """
    PDBファイルの読み込み
    
    Args:
        pdb_file: PDBファイルのパス
    
    Returns:
        str: PDBファイルの内容、失敗時はNone
    """
    try:
        if not pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
            
        with open(pdb_file, 'r') as f:
            return f.read()
    
    except Exception as e:
        logger.error(f"Error loading PDB file: {e}")
        return None

def display_protein_structure(pdb_content: str, 
                            style: Dict = None,
                            width: int = 600,
                            height: int = 400):
    """
    タンパク質構造の表示
    
    Args:
        pdb_content: PDBファイルの内容
        style: 表示スタイル設定
        width: 表示幅
        height: 表示高さ
    """
    try:
        view = py3Dmol.view(width=width, height=height)
        view.addModel(pdb_content, "pdb")
        
        if style is None:
            style = {'cartoon': {'color': 'spectrum'}}
        
        view.setStyle(style)
        view.zoomTo()
        html(view._make_html(), height=height)
        
    except Exception as e:
        st.error(f"Error displaying protein structure: {e}")

def protein_view_section(pdb_file: Optional[pathlib.Path] = None,
                        docking_result: Optional[Dict] = None):
    """
    タンパク質表示セクションの表示
    
    Args:
        pdb_file: PDBファイルのパス
        docking_result: ドッキング結果（オプション）
    """
    st.subheader("Protein Structure")
    
    if pdb_file and pdb_file.exists():
        # 表示スタイルの選択
        style = st.selectbox(
            "Display Style",
            ["Cartoon", "Ribbon", "Line", "Surface"],
            key="protein_style"
        )
        
        # 色設定の選択
        color_scheme = st.selectbox(
            "Color Scheme",
            ["Spectrum", "Chain", "Secondary Structure", "Element"],
            key="color_scheme"
        )
        
        # PDB内容の読み込み
        pdb_content = load_pdb_content(pdb_file)
        if pdb_content:
            view = py3Dmol.view(width=600, height=400)
            view.addModel(pdb_content, "pdb")
            
            # スタイルの適用
            style_settings = {}
            if style == "Cartoon":
                style_settings = {'cartoon': {}}
            elif style == "Ribbon":
                style_settings = {'ribbon': {}}
            elif style == "Line":
                style_settings = {'line': {}}
            else:  # Surface
                style_settings = {'surface': {}}
            
            # 色の設定
            if color_scheme == "Spectrum":
                style_settings[list(style_settings.keys())[0]]['color'] = 'spectrum'
            elif color_scheme == "Chain":
                style_settings[list(style_settings.keys())[0]]['colorscheme'] = {'prop': 'chain'}
            elif color_scheme == "Secondary Structure":
                style_settings[list(style_settings.keys())[0]]['colorscheme'] = {'prop': 'ss'}
            else:  # Element
                style_settings[list(style_settings.keys())[0]]['colorscheme'] = {'prop': 'element'}
            
            # ドッキング結果がある場合
            if docking_result:
                # 結合ポケットの表示
                pocket_residues = docking_result.get('pocket_residues', [])
                if pocket_residues:
                    view.addStyle({'resi': pocket_residues}, 
                                {'stick': {'color': 'orange'}})
                
                # リガンドの表示
                ligand_pose = docking_result.get('best_pose')
                if ligand_pose:
                    view.addModel(ligand_pose, "pdb")
                    view.setStyle({'model': -1}, {'stick': {'color': 'green'}})
            
            view.setStyle({}, style_settings)
            view.zoomTo()
            html(view._make_html(), height=400)
            
            # タンパク質情報の表示
            st.markdown("""
            **Structure Information:**
            - File: {}
            - Chains: {}
            """.format(pdb_file.name, "A, B, C"))  # チェーン情報は実際のデータから取得する必要あり
            
        else:
            st.error("Failed to load protein structure")
    else:
        st.info("Select a protein from the results table to view structure")