import streamlit as st
import pathlib
import sys
from pathlib import Path

# プロジェクトルートをPythonパスに追加
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

from app.components.structure_input import structure_input_section
from app.components.results_view import results_view_section
from app.components.ligand_view import ligand_view_section
from app.components.protein_view import protein_view_section
from app.components.settings import settings_section, initialize_settings
from app.dashboard import DashboardManager, mini_dashboard
from core.database.protein_db import ProteinDatabase
from core.database.results_db import ResultsDatabase
from core.docking.calculator import DockingCalculator

# アプリケーションの設定
st.set_page_config(
    page_title="ChemDockFinder",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

def initialize_session():
    """セッション状態の初期化"""
    if 'protein_db' not in st.session_state:
        st.session_state.protein_db = ProteinDatabase()
    if 'results_db' not in st.session_state:
        st.session_state.results_db = ResultsDatabase()
    if 'calculator' not in st.session_state:
        st.session_state.calculator = DockingCalculator()
    if 'selected_protein' not in st.session_state:
        st.session_state.selected_protein = None
    if 'current_calculation' not in st.session_state:
        st.session_state.current_calculation = None

def sidebar():
    """サイドバーの表示"""
    with st.sidebar:
        st.title("ChemDockFinder")
        
        # ナビゲーション
        page = st.radio(
            "Navigation",
            ["Docking", "Dashboard", "Settings"]
        )
        
        st.divider()
        
        # ミニダッシュボード
        if page != "Dashboard":
            mini_dashboard()
        
        return page

def run_docking(smiles: str, protein_file: pathlib.Path):
    """ドッキング計算の実行"""
    try:
        with st.spinner("Preparing calculation..."):
            # 計算の準備
            calc_id = st.session_state.calculator.prepare_calculation(
                smiles=smiles,
                protein_files=[protein_file],
                params=st.session_state.docking_settings
            )
            st.session_state.current_calculation = calc_id
        
        with st.spinner("Running docking simulation..."):
            # 計算の実行
            results = st.session_state.calculator.run_docking(
                calc_id,
                lambda progress: st.progress(progress / 100.0)
            )
            
            # 結果の保存
            st.session_state.results_db.add_calculation_result(calc_id, results)
            
            return results
            
    except Exception as e:
        st.error(f"Error in docking calculation: {e}")
        return None

def docking_page():
    """ドッキングページの表示"""
    # 上部グリッド
    top_left, top_right = st.columns(2, gap="medium")
    
    with top_left:
        # 構造入力エリア
        smiles, input_method = structure_input_section()
        
    with top_right:
        # 検索結果エリア
        results_df = st.session_state.protein_db.get_all_proteins()
        results_view_section(results_df)
    
    # 下部グリッド
    bottom_left, bottom_right = st.columns(2, gap="medium")
    
    with bottom_left:
        # リガンド表示エリア
        if smiles:
            ligand_view_section(smiles)
            
            # 計算コントロール
            st.subheader("Simulation Control")
            control_col1, control_col2 = st.columns(2)
            with control_col1:
                if st.button("Run Docking", use_container_width=True, type="primary"):
                    if st.session_state.selected_protein:
                        results = run_docking(smiles, st.session_state.selected_protein)
                        if results:
                            st.success("Docking completed successfully!")
                    else:
                        st.warning("Please select a protein from the results table")
            
            with control_col2:
                if st.button("⚙️ Settings", use_container_width=True):
                    settings_section()
    
    with bottom_right:
        # タンパク質構造表示エリア
        if st.session_state.selected_protein:
            protein_view_section(
                st.session_state.selected_protein,
                st.session_state.current_calculation
            )

def main():
    """メイン関数"""
    # 初期化
    initialize_session()
    initialize_settings()
    
    # ダッシュボードマネージャーの設定
    dashboard_manager = DashboardManager(st.session_state.results_db)
    
    # サイドバーとページ選択
    page = sidebar()
    
    if page == "Docking":
        docking_page()
    
    elif page == "Dashboard":
        dashboard_manager.display_dashboard()
    
    elif page == "Settings":
        st.title("Settings")
        settings_section()

if __name__ == "__main__":
    main()