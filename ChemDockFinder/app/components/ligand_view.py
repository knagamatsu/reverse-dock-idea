import streamlit as st
import py3Dmol
from streamlit.components.v1 import html
from rdkit import Chem
from rdkit.Chem import AllChem
import logging
from typing import Optional

logger = logging.getLogger(__name__)

def generate_3d_structure(smiles: str) -> Optional[str]:
    """
    SMILES文字列から3D構造を生成
    
    Args:
        smiles: SMILES文字列
    
    Returns:
        str: PDB形式の構造データ、失敗時はNone
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
            
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        return Chem.MolToPDBBlock(mol)
    
    except Exception as e:
        logger.error(f"Error generating 3D structure: {e}")
        return None

def display_ligand_view(smiles: str, width: int = 400, height: int = 400):
    """
    リガンドの3D表示
    
    Args:
        smiles: SMILES文字列
        width: 表示幅
        height: 表示高さ
    """
    try:
        pdb_data = generate_3d_structure(smiles)
        if pdb_data:
            view = py3Dmol.view(width=width, height=height)
            view.addModel(pdb_data, "pdb")
            view.setStyle({'stick':{}})
            view.zoomTo()
            html(view._make_html(), height=height)
        else:
            st.error("Failed to generate 3D structure")
            
    except Exception as e:
        st.error(f"Error displaying structure: {e}")

def ligand_view_section(smiles: Optional[str]):
    """
    リガンド表示セクションの表示
    
    Args:
        smiles: SMILES文字列
    """
    st.subheader("3D Ligand View")
    
    if smiles:
        # 表示スタイルの選択
        style = st.selectbox(
            "Display Style",
            ["Stick", "Ball and Stick", "Space Fill"],
            key="ligand_style"
        )
        
        # 3D表示
        try:
            pdb_data = generate_3d_structure(smiles)
            if pdb_data:
                view = py3Dmol.view(width=400, height=400)
                view.addModel(pdb_data, "pdb")
                
                # スタイルの適用
                if style == "Stick":
                    view.setStyle({'stick':{}})
                elif style == "Ball and Stick":
                    view.setStyle({'stick':{}, 'sphere':{'radius':0.5}})
                else:  # Space Fill
                    view.setStyle({'sphere':{}})
                
                view.zoomTo()
                html(view._make_html(), height=400)
                
                # 構造情報の表示
                mol = Chem.MolFromSmiles(smiles)
                st.markdown(f"""
                **Molecule Information:**
                - Atoms: {mol.GetNumAtoms()}
                - Bonds: {mol.GetNumBonds()}
                - Molecular Weight: {round(Chem.Descriptors.ExactMolWt(mol), 2)}
                """)
                
            else:
                st.error("Failed to generate 3D structure")
                
        except Exception as e:
            st.error(f"Error displaying structure: {e}")
    else:
        st.info("Enter a chemical structure to view 3D model")