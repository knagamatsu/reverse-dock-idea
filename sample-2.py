import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
import pandas as pd
import plotly.graph_objects as go
import py3Dmol
from dockstring import load_target

# ページ設定
st.set_page_config(page_title="リガンドドッキングシミュレーター", layout="wide")
st.title("🧬 リガンド-タンパク質ドッキングシミュレーター")

# デフォルトのリガンドデータ
DEFAULT_LIGANDS = {
    "Anagliptin": "CC1=NN2C=C(C=NC2=C1)C(=O)NCC(C)(C)NCC(=O)N3CCC[C@H]3C#N",
    "Omarigriptin": "CS(=O)(=O)n1cc2c(n1)CN([C@H]1CO[C@H](c3cc(F)ccc3F)[C@@H](N)C1)C2",
    "Voglibose": "C1[C@@H]([C@@H]([C@H]([C@@H]([C@]1(CO)O)O)O)O)NC(CO)CO"
}

def calculate_properties(mol):
    """分子の物性を計算"""
    if mol is not None:
        return {
            "分子量": Descriptors.ExactMolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "水素結合ドナー": Descriptors.NumHDonors(mol),
            "水素結合アクセプター": Descriptors.NumHAcceptors(mol),
            "TPSA": Descriptors.TPSA(mol),
            "回転可能結合数": Descriptors.NumRotatableBonds(mol)
        }
    return None

# サイドバーに入力方法を配置
st.sidebar.header("入力設定")
input_method = st.sidebar.radio(
    "入力方法を選択",
    ["既存リガンドから選択", "SMILES式を入力", "構造式エディタで描画"]
)

# 入力セクション
smiles = None
if input_method == "既存リガンドから選択":
    ligand_name = st.sidebar.selectbox("リガンドを選択", list(DEFAULT_LIGANDS.keys()))
    smiles = DEFAULT_LIGANDS[ligand_name]
    st.sidebar.info(f"SMILES: {smiles}")

elif input_method == "SMILES式を入力":
    smiles = st.sidebar.text_input("SMILES式を入力", "")

else:  # 構造式エディタ
    # サイドバーは狭いので、メイン領域にエディタを配置
    st.subheader("構造式エディタ")
    smiles = st_ketcher("")
    if smiles:
        st.sidebar.info(f"SMILES: {smiles}")

# メイン領域をカラムで分割
if smiles:
    col1, col2 = st.columns([3, 2])
    
    with col1:
        # 構造表示と物性値
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            subcol1, subcol2 = st.columns([1, 1])
            with subcol1:
                st.subheader("2D構造")
                img = Draw.MolToImage(mol)
                st.image(img)

            with subcol2:
                st.subheader("物性値")
                props = calculate_properties(mol)
                for name, value in props.items():
                    st.metric(name, f"{value:.2f}" if isinstance(value, float) else value)

            # 3D構造の生成と表示
            st.subheader("3D構造")
            try:
                mol_3d = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_3d, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol_3d)
                
                # mol_blockを生成
                mol_block = Chem.MolToMolBlock(mol_3d)
                
                # py3Dmolビューアの設定
                viewer = py3Dmol.view(width=600, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({'stick':{}})
                viewer.zoomTo()
                
                # HTMLとしてビューアを表示
                html = viewer._repr_html_()
                st.components.v1.html(html, height=400)
            except Exception as e:
                st.error(f"3D構造の生成に失敗しました: {str(e)}")
    
    with col2:
        # ドッキングシミュレーション
        st.subheader("ドッキングシミュレーション")
        if st.button("ドッキングを実行", type="primary"):
            try:
                with st.spinner("ドッキング計算を実行中..."):
                    # ドッキング計算
                    target = load_target('DPP4')
                    score, aux = target.dock(smiles)
                    
                    # 結果の表示
                    score_color = "green" if score < -7.0 else "orange" if score < -5.0 else "red"
                    st.metric("ドッキングスコア", f"{score:.2f}")
                    
                    # スコアの解釈
                    interpretation = {
                        "強い結合": score < -7.0,
                        "中程度の結合": -7.0 <= score < -5.0,
                        "弱い結合": score >= -5.0
                    }
                    
                    for desc, condition in interpretation.items():
                        if condition:
                            st.markdown(f"**結果**: :{score_color}[{desc}が予測されます]")
                    
                    # スコアの詳細説明
                    st.markdown("""
                    **スコアの解釈**:
                    - -7.0以下: 強い結合
                    - -7.0～-5.0: 中程度の結合
                    - -5.0以上: 弱い結合
                    """)
                    
                    # ドッキングポーズの保存と表示
                    st.subheader("ドッキングポーズ")
                    # ここにドッキング結果の3D表示コードを追加

            except Exception as e:
                st.error(f"ドッキング計算でエラーが発生しました: {str(e)}")

else:
    st.info("リガンドを入力してください")

# フッター
st.sidebar.markdown("---")
st.sidebar.markdown("""
### 注意事項
- このツールは教育・研究目的で使用してください
- 予測結果は参考値であり、実際の結合活性とは異なる場合があります
""")
