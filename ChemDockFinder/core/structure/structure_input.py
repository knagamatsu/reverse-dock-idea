import streamlit as st
from streamlit_ketcher import st_ketcher
from typing import Tuple, Optional
import pandas as pd
from rdkit import Chem

def validate_smiles(smiles: str) -> bool:
    """SMILES形式の検証"""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def structure_input_section() -> Tuple[str, str]:
    """
    構造入力セクションの表示
    
    Returns:
        Tuple[str, str]: (SMILES文字列, 入力方法)
    """
    # 入力方法の選択
    input_method = st.radio(
        "Input Method",
        ["Draw Structure", "SMILES Input", "File Upload"],
        horizontal=True
    )
    
    smiles = None
    method = None
    
    if input_method == "Draw Structure":
        # Ketcher エディタによる構造描画
        default_smiles = "NCC1CCCCC1"  # アミノメチルシクロヘキサン
        smiles = st_ketcher(default_smiles, height=450)
        method = "drawn"
        
    elif input_method == "SMILES Input":
        # SMILES入力
        smiles = st.text_input(
            "Enter SMILES",
            value="NCC1CCCCC1",
            help="Enter a valid SMILES string"
        )
        method = "smiles"
        
    else:  # File Upload
        # ファイルアップロード
        uploaded_file = st.file_uploader(
            "Upload Structure File",
            type=['smi', 'mol', 'sdf'],
            help="Upload a file containing chemical structure"
        )
        
        if uploaded_file:
            try:
                # ファイルの読み込みと処理
                content = uploaded_file.read().decode('utf-8')
                if uploaded_file.name.endswith('.smi'):
                    smiles = content.strip()
                else:
                    # MOLファイルからSMILESへの変換
                    mol = Chem.MolFromMolBlock(content)
                    if mol:
                        smiles = Chem.MolToSmiles(mol)
                method = "file"
            except Exception as e:
                st.error(f"Error reading file: {e}")
                return None, None
    
    # SMILES の検証
    if smiles:
        if validate_smiles(smiles):
            st.success("Valid structure")
            st.markdown(f"**SMILES:** `{smiles}`")
            return smiles, method
        else:
            st.error("Invalid structure")
            return None, None
    
    return None, None